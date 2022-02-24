.bootintLIN <- function(data, idx , x , model , var1 , var2){
  df <- data[idx, ]
  mycall <- model$call
  mycall <- pryr::modify_call(mycall , list(data=quote(df)))
  # myformula <- model$sformula
  mymodel <- eval(mycall)
  # mymodel <- cph(model$sformula, data=df)
  coefMod <- coef(mymodel)
  if("Glm" %in% class(model)){
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
  LIN <- unname( var1Term + x*intTerm )
  LIN
}


#' Linear regression interaction estimates
#'
#' Generate linear estimates for a 1 unit increase in a variable at
#' specified points of another interacting variable in a linear interaction model
#'
#' @param var2values numeric vector of var2 points to estimate
#' @param model model of class rms::Glm or stats::glm family gaussian. If data is NULL, the function expects to find the data in model$x
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
#' # Recode diabetes as 0/1
#' PimaIndiansDiabetes$diabetes <- ifelse(PimaIndiansDiabetes$diabetes=="pos" , 1 , 0)
#' myformula <- glucose ~ mass + diabetes * age
#' model <- glm(myformula , data = PimaIndiansDiabetes ,family=gaussian)
#' # Show the effect on glucose of being diabetic at age 20 to 80
#' linLIN( var2values = 20:80
#'        , model = model , data = PimaIndiansDiabetes , var1 ="diabetes", var2="age"
#'        , ci=TRUE , conf = 0.95 , ci.method = "delta")
#' @return if ci = FALSE, a vector of estimate of length(var2values),
#' if ci = TRUE a dataframe with 5 columns, initial values, linear estimates, lower CI, upper CI and SE
#' @importFrom rms Glm
#' @importFrom pryr modify_call
#' @importFrom msm deltamethod
#' @importFrom boot boot boot.ci
#' @importFrom stats vcov coef as.formula qnorm sd glm
#' @export
linLIN <- function(var2values , model , data , var1 , var2
                     , ci=TRUE , conf = 0.95 , ci.method = "delta"
                     , ci.boot.method = "perc" , R = 100 , parallel = "multicore" , ...) {
  # Check correct class for model
  if( !any( c("Glm","glm") %in% class(model) ) ){
    stop("Interaction linear model must be run with rms::Glm or stats::glm")
  }
  # Check correct family
  if(!"gaussian" %in% model$family$family){
    stop("model of class glm but not family gaussian")
  } else {
    if("glm" %in% class(model) && !"Glm" %in% class(model)){
      modelClass <- "glm"
    } else {
      modelClass <- "Glm"
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
  if("Glm" %in% class(model)){
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
  LIN <- unname( var1Term + x*intTerm)
  if(ci){
    alpha <- qnorm( 1 - (1-conf)/2)
    if(ci.method=="normal"){
      stop("not implemented yet")
      SEgeneral <- sqrt( sum( c( diag(model$var)[c(mycoefWhich[var1] ,mycoefWhich[ mycoefWhich[ myvars[3]] ])]
                                 , 2*model$var[mycoefWhich[var1] ,mycoefWhich[ myvars[3]]]) ) )
      SEvar1 <- sqrt(diag(model$var))[mycoefWhich[var1]]
      # wrong
      LIN_L <- unname(exp( var1Term - alpha*SEvar1 + x*(intTerm- alpha*SEgeneral) ))
      LIN_U <- unname(exp( var1Term + alpha*SEvar1 + x*(intTerm+ alpha*SEgeneral) ))
      # also wrong
      # LIN_L <- unname(exp( var1Term + x*(intTerm) - alpha*SEgeneral ))
      # LIN_U <- unname(exp( var1Term + x*(intTerm) + alpha*SEgeneral ))
      # also wrong
      # LIN_L <- unname(exp( var1Term + x*(intTerm- alpha*SEgeneral) ))
      # LIN_U <- unname(exp( var1Term + x*(intTerm+ alpha*SEgeneral) ))
      LINci <- data.frame(Value = x
                         , LIN = LIN
                         , CI_L = LIN_L
                         , CI_U = LIN_U
                         , SE = SEgeneral)
      colnames(LINci) <- c("Value" , "LIN" , "CI_L" , "CI_U" , "SE")
      LINci <- as.data.frame(LINci)
      # class(LINci) <- c("LIN" , class(LINci))
      return(LINci)
    } else if(ci.method == "delta"){
      # This creates a vector like x1 , x2 , x3 , x7 , x8
      # that tells you the position of the regressor as it appears in the model
      xNum <- paste0("x" , mycoefWhich)
      vcovMod <- vcov(model)
      LINci <- t(vapply( seq_len(length(x)) , function(i) {
        x_i <- x[i]
        myform <- paste0("~(", xNum[1] , " + " , xNum[3] , "*(" , x_i , "))")
        SE <- NULL
        try(SE<-msm::deltamethod(as.formula(myform), coefMod, vcovMod) , silent = TRUE)
        if(is.null(SE)){
          return(c(LIN[i] , NA , NA , NA))
        }
        up<-LIN[i]+alpha*SE
        lo<-LIN[i]-alpha*SE
        c(LIN[i] , lo , up , SE)
      } , numeric(4)))
      LINci <- cbind( Value = x , LINci)
      rownames(LINci) <- x
      colnames(LINci) <- c("Value" , "LIN" , "CI_L" , "CI_U" , "SE")
      LINci <- as.data.frame(LINci)
      # class(LINci) <- c("LIN" , class(LINci))
      return(LINci)
    } else if(ci.method == "bootstrap"){
      if(missing(parallel)){
        parallel <- "multicore"
      }
      if(missing(R)){
        R <- 100
      }
      myBoot <- boot::boot(data = data, statistic = .bootintLIN
                           , x = x , model = model
                           , R = R , parallel = parallel
                           , var1 = var1 , var2 = var2 )
      SE <- apply(myBoot$t , 2 , sd)
      LINci <- t(vapply( seq_len(length(x)) , function(idx) {
        bci <- boot::boot.ci(boot.out = myBoot,  index = idx , type = ci.boot.method , conf = conf)
        if(ci.boot.method %in% "norm"){
          c(bci$t0 , bci$normal[2] , bci$normal[3] , SE = SE[idx])
        } else {
          c(bci$t0 , bci[[4]][4] , bci[[4]][5] , SE = SE[idx])
        }
      } , numeric(4)))
      LINci <- cbind(x , LINci)
      colnames(LINci) <- c("Value" , "LIN" , "CI_L" , "CI_U" , "SE")
      rownames(LINci) <- x
      LINci <- as.data.frame(LINci)
      # class(LINci) <- c("LIN" , class(LINci))
      return(LINci)
    } else {
      stop("Only delta and bootstrap are valid CI methods")
    }
  } else {
    LIN <- data.frame(Value = x , LIN = LIN)
    # class(LIN) <- c("LIN",class(LIN))
    return(LIN)
  }
}
