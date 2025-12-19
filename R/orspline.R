.bootintOR <- function(data, idx , x , model , var1 , var2){
  df <- data[idx, ]
  mycall <- model$call
  # mycall <- pryr::modify_call(mycall , list(data=quote(df)))
  mycall <- rlang::call_modify(mycall , data=rlang::expr(df))
  # myformula <- model$sformula
  mymodel <- eval(mycall)
  # mymodel <- cph(model$sformula, data=df)
  coefMod <- coef(mymodel)
  if("lrm" %in% class(model)){
    separator <- " * "
  } else {
    separator <- ":"
  }
  myvars <- intersect( c(var1 , var2
                         , paste(var1 , var2 , sep = separator)
                         , paste(var2 , var1 , sep = separator))
                       , names(coefMod))
  mycoef <- coefMod[  myvars ]
  mycoefWhich <- sapply( myvars , function(v) which( names(coefMod) %in% v ))
  intTerm <- mycoef[ length(myvars) ]
  var1Term <- mycoef[ var1 ]
  var2Term <- mycoef[ var2 ]
  OR <- unname(exp( var1Term + x*intTerm))
  OR
}


#' Linear interaction OR
#'
#' Generate OR values for a 1 unit increase in a variable at
#' specified points of another interacting variable in a simple logistic interaction model
#'
#' @param var2values numeric vector of var2 points to estimate
#' @param model model of class lrm Glm or glm. If data is NULL, the function expects to find the data in model$x
#' @param data data used in the model. If absent, it will attempt to recover the data from the model object. Only used for bootstrap CI
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
#' myformula <- diabetes ~ mass + age * glucose
#' model <- glm(myformula , data = PimaIndiansDiabetes , family = binomial())
#' loglinOR( var2values = 20:80
#'        , model = model , data = PimaIndiansDiabetes , var1 ="age", var2="glucose"
#'        , ci=TRUE , conf = 0.95 , ci.method = "delta")
#' @return if ci = FALSE, a vector of estimate of length(var2values),
#' if ci = TRUE a dataframe with 5 columns, initial values, OR, lower CI, upper CI and SE
#' @importFrom rms cph
#' @importFrom rlang call_modify expr
#' @importFrom msm deltamethod
#' @importFrom boot boot boot.ci
#' @importFrom stats vcov coef as.formula qnorm sd glm
#' @export
loglinOR <- function(var2values , model , data , var1 , var2
                     , ci=TRUE , conf = 0.95 , ci.method = "delta"
                     , ci.boot.method = "perc" , R = 100 , parallel = "multicore" , ...) {
  # Check correct class for model
  if( !any( c("lrm","glm","Glm") %in% class(model) ) ){
    stop("Interaction Cox model must be run with rms::lrm, rms::Glm or stats::glm")
  }
  # Check correct family
  if(!"lrm" %in% class(model)){
    if(!"binomial" %in% model$family$family){
      stop("model of class glm but not family binomial")
    }
  }
  if(!is.numeric(var2values)){
    stop("var2values must be a numeric vector")
  }
  x <- var2values
  if(missing(data)){
    if(is.null(model$x) && ci && ci.method == "bootstrap"){
      # We need the data only for bootstrap CI
      stop("Missing data for bootstrap CI")
    } else {
      data <- model$x
    }
  }
  if(!all(c(var1,var2) %in% colnames(data) )){
    stop("var1 or var2 not present in the data")
  }
  coefMod <- coef(model)
  # Different styles to name the interaction term between lrm and glm
  if(any(c("lrm","Glm") %in% class(model))){
    separator <- " * "
  } else {
    separator <- ":"
  }
  myvars <- intersect( c(var1 , var2
                         , paste(var1 , var2 , sep = separator)
                         , paste(var2 , var1 , sep = separator))
                       , names(coefMod))
  if(length(myvars)!=3) stop("either var1 or var2 is not in the interaction!")
  mycoef <- coefMod[  myvars ]
  mycoefWhich <- vapply( myvars , function(v) which( names(coefMod) %in% v ) , numeric(1))
  intTerm <- mycoef[ length(myvars) ]
  var1Term <- mycoef[ var1 ]
  var2Term <- mycoef[ var2 ]
  OR <- unname(exp( var1Term + x*intTerm))
  if(ci){
    alpha <- qnorm( 1 - (1-conf)/2)
    if(ci.method=="normal"){
      stop("not implemented yet")
      SEgeneral <- sqrt( sum( c( diag(model$var)[c(mycoefWhich[var1] ,mycoefWhich[ mycoefWhich[ myvars[3]] ])]
                                 , 2*model$var[mycoefWhich[var1] ,mycoefWhich[ myvars[3]]]) ) )
      SEvar1 <- sqrt(diag(model$var))[mycoefWhich[var1]]
      # wrong
      OR_L <- unname(exp( var1Term - alpha*SEvar1 + x*(intTerm- alpha*SEgeneral) ))
      OR_U <- unname(exp( var1Term + alpha*SEvar1 + x*(intTerm+ alpha*SEgeneral) ))
      # also wrong
      # OR_L <- unname(exp( var1Term + x*(intTerm) - alpha*SEgeneral ))
      # OR_U <- unname(exp( var1Term + x*(intTerm) + alpha*SEgeneral ))
      # also wrong
      # OR_L <- unname(exp( var1Term + x*(intTerm- alpha*SEgeneral) ))
      # OR_U <- unname(exp( var1Term + x*(intTerm+ alpha*SEgeneral) ))
      ORci <- data.frame(Value = x
                         , OR = OR
                         , CI_L = OR_L
                         , CI_U = OR_U
                         , SE = SEgeneral)
      colnames(ORci) <- c("Value" , "OR" , "CI_L" , "CI_U" , "SE")
      ORci <- as.data.frame(ORci)
      # class(ORci) <- c("OR" , class(ORci))
      return(ORci)
    } else if(ci.method == "delta"){
      # This creates a vector like x1 , x2 , x3 , x7 , x8
      # that tells you the position of the regressor as it appears in the model
      xNum <- paste0("x" , mycoefWhich)
      vcovMod <- vcov(model)
      ORci <- t(vapply( seq_len(length(x)) , function(i) {
        x_i <- x[i]
        myform <- paste0("~(", xNum[1] , " + " , xNum[3] , "*(" , x_i , "))")
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
      # class(ORci) <- c("OR" , class(ORci))
      return(ORci)
    } else if(ci.method == "bootstrap"){
      if(missing(parallel)){
        parallel <- "multicore"
      }
      if(missing(R)){
        R <- 100
      }
      myBoot <- boot::boot(data = data, statistic = .bootintOR
                           , x = x , model = model
                           , R = R , parallel = parallel
                           , var1 = var1 , var2 = var2 )
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
    # class(OR) <- c("OR",class(OR))
    return(OR)
  }
}
