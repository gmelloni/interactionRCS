#' Restricted cubic spline interaction OR
#'
#' Generate OR values in a logistic model for a 1 unit increase in a variable at
#' specified points of another interacting variable splined with rcs(df = 3)
#'
#' @param var2values numeric vector of var2 points to estimate
#' @param model model of class lrm. If data is NULL, the function expects to find the data in model$x.
#' @param data data used in the model. If absent, we will attempt to recover the data from the model. Only used for bootstrap
#' @param var1 variable that increases by 1 unit from 0
#' @param var2 variable to spline. var2values belong to var2
#' @param ci calculate 95% CI?
#' @param conf confidence level. Default 0.95
#' @param ci.method confidence interval method. "delta" performs delta method. "bootstrap" performs bootstrapped CI (slower)
#' @param ci.boot.method one of the available bootstrap CI methods from \code{\link[boot]{boot.ci}}. Default percentile
#' @param R number of bootstrap samples if ci.method = "bootstrap". Default 100
#' @param parallel can take values "no", "multicore", "snow" if ci.method = "bootstrap". Default multicore
#' @param ... other parameters for boot
#' @examples
#' library(rms)
#' library(mlbench)
#' data(PimaIndiansDiabetes)
#' # Set age on a 5-year scale
#' PimaIndiansDiabetes$age <- PimaIndiansDiabetes$age/5
#' # Recode diabetes as 0/1
#' PimaIndiansDiabetes$diabetes <- ifelse(PimaIndiansDiabetes$diabetes=="pos" , 1 , 0)
#' myformula <- diabetes ~ mass + age * rcs(glucose, 3)
#' model <- lrm(myformula , data = PimaIndiansDiabetes )
#' rcsOR( var2values = 20:80
#'        , model = model , data = PimaIndiansDiabetes , var1 ="age", var2="glucose"
#'        , ci=TRUE , conf = 0.95 , ci.method = "delta")
#' @return if ci = FALSE, a dataframe with initial values and OR
#' , if ci = TRUE a dataframe with 5 columns, initial values, OR, lower CI, upper CI and SE
#' @importFrom rms lrm rcs
#' @importFrom pryr modify_call
#' @importFrom msm deltamethod
#' @importFrom boot boot boot.ci
#' @importFrom stats vcov coef as.formula qnorm sd
#' @export
rcsOR <- function(var2values , model , data=NULL , var1 , var2
                  , ci=TRUE , conf = 0.95 , ci.method = "delta"
                  , ci.boot.method = "perc" , R = 100 , parallel = "multicore" , ...) {
  # argg <- c(as.list(environment()), list(...))
  # Check correct class for model
  if( !all( c("lrm","rms","glm") %in% class(model) ) ){
    stop("Cubic spline Logistic model must be run with rms::lrm")
  }
  if(!is.numeric(var2values)){
    stop("var2values must be a numeric vector")
  }
  x <- var2values
  if(missing(data)){
    if(is.null(model$x)){
      stop("Missing data")
    } else {
      data <- model$x
    }
  }
  if(!all(c(var1,var2) %in% colnames(data) )){
    stop("var1 or var2 not present in the data")
  }
  # Check that var1 is a 0/1, if not check if the mean is 0
  # if(!all(data[[var1]] %in% c(0,1,NA))){
  #   if(!isTRUE(all.equal(mean(data[[var1]] , na.rm =TRUE),0))){
  #     warning("var1 is not centered on 0 nor a 0/1 variable, results are always reported for a 0 to 1 change in var1.")
  #   }
  # }

  coefMod <- coef(model)
  # k <- attributes(rms::rcs(data[[var2]], 3))$parms
  k <- model$Design$parms[[var2]]
  k2k1 <- (k[2] - k[1])/(k[3] - k[2])
  k3k1 <- (k[3] - k[1])/(k[3] - k[2])
  # Extract parameters
  # beta1 <- coefMod["score.sd"]
  # alph <- coefMod[c("age" , "age'")]
  # lamb <- coefMod[c("score.sd * age" , "score.sd * age'")]
  # Could be rcs(var2)*var1 or var1*rcs(var2). We search for both versions
  myvars <- intersect( c(var1 , var2
                         , paste0(var2 , "'")
                         , paste(var1 , var2 , sep = " * ")
                         , paste(var2 , var1 , sep = " * ")
                         , paste(var1 , paste0(var2 , "'") , sep = " * ")
                         , paste(paste0(var2 , "'"),var1 , sep = " * "))
                       , names(coefMod))
  if(length(myvars)!=5) stop("either var1 or var2 is not in the interaction!")
  mycoef <- coefMod[  myvars ]
  mycoefWhich <- sapply( myvars , function(v) which( names(coefMod) %in% v ))
  a <- mycoef[ c(var2 , paste0(var2,"'"))]
  b <- mycoef[ var1 ]
  l <- mycoef[ setdiff(myvars , c(var2 , paste0(var2,"'") , var1)) ]
  numer <- vapply(x , function(i) {
    max(i - k[1],0)^3  - (max(i - k[2],0)^3)*k3k1 + (max(i - k[3],0)^3)*k2k1
  } , numeric(1))
  denom <- (k[3] - k[1])^2
  numDem <- numer/denom
  sp1 <- vapply(x , function(i) a[1]*i , numeric(1)) + a[2]*numDem
  sp2 <- vapply(x , function(i) l[1]*i , numeric(1)) + l[2]*numDem
  OR <- unname(exp( b + sp2))
  # OR <- unname(exp( b + sp1 + sp2)/exp(sp1))
  # OR <- unname(exp( b*x1 + sp1 + sp2*x1)/exp(b*x2 + sp1 + sp2*x2))

  if(ci){
    alpha <- qnorm( 1 - (1-conf)/2)
    if(ci.method == "delta"){
      # This creates a vector like x1 , x2 , x3 , x7 , x8
      # that tells you the position of the regressor as it appears in the model
      xNum <- paste0("x" , mycoefWhich)
      vcovMod <- vcov(model)
      ORci <- t(vapply( seq_len(length(x)) , function(i) {
        x_i <- x[i]
        numDem_i <- numDem[i]
        # myform <- paste0("~(", xNum[1] , " + " , xNum[2] , "*(" , x_i
        #                , ") + " , xNum[3] , "*(" , numDem_i
        #                , ") + " , xNum[4] , "*(" , x_i
        #                , ") + " , xNum[5] , "*(" , numDem_i
        #                , "))/(" , xNum[2] , "*(" , x_i
        #                , ") + " , xNum[3] , "*(" , numDem_i , "))")
        myform <- paste0("~(", xNum[1] , " + ", xNum[4] , "*(" , x_i
                         , ") + " , xNum[5] , "*(" , numDem_i
                         , "))")
        SE <- NULL
        try(SE<-msm::deltamethod(as.formula(myform), coefMod, vcovMod) , silent = TRUE)
        if(is.null(SE)){
          return(c(OR[i] , NA , NA , NA))
        }
        up<-exp(log(OR[i])+alpha*SE)
        lo<-exp(log(OR[i])-alpha*SE)
        c(OR[i] , lo , up , SE)
      } , numeric(4)))
      ORci <- cbind( Value = x , ORci)
      rownames(ORci) <- x
      colnames(ORci) <- c("Value" , "OR" , "CI_L" , "CI_U" , "SE")
      ORci <- as.data.frame(ORci)
      class(ORci) <- c("OR" , class(ORci))
      return(ORci)
    } else if(ci.method == "bootstrap"){
      if(missing(parallel)){
        parallel <- "multicore"
      }
      if(missing(R)){
        R <- 100
      }
      myBoot <- boot::boot(data = data, statistic = .bootrcsHR, x = x , model = model
                           , R = R , parallel = parallel, var1 = var1 , var2 = var2 )
      SE <- apply(myBoot$t , 2 , sd)
      ORci <- t(vapply( seq_len(length(x)) , function(idx) {
        bci <- boot::boot.ci(boot.out = myBoot,  index = idx , type = ci.boot.method , conf = conf)
        if(ci.boot.method %in% "norm"){
          c(bci$t0 , bci$normal[2] , bci$normal[3] , SE = SE[idx])
        } else {
          c(bci$t0 , bci[[4]][4] , bci[[4]][5] , SE = SE[idx])
        }
      } , numeric(4)))
      ORci <- cbind(x , ORci)
      colnames(ORci) <- c("Value" , "OR" , "CI_L" , "CI_U" , "SE")
      rownames(ORci) <- x
      ORci <- as.data.frame(ORci)
      # class(ORci) <- c("OR" , class(ORci))
      return(ORci)
    } else {
      stop("Only delta and bootstrap are valid CI methods")
    }
  } else {
    OR <- data.frame(Value = x , OR = OR)
    # class(OR) <- c("OR",class(HR))
    return(OR)
  }
}
