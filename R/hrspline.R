.bootHRSpline <- function(data, idx , x , model , var1 , var2){
  df <- data[idx, ]
  mycall <- model$call
  mycall <- pryr::modify_call(mycall , list(data=df))
  mymodel <- eval(mycall)
  coefMod <- coef(mymodel)
  if("cph" %in% class(model)){
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
  HR <- unname(exp( var1Term + x*intTerm))
  HR
}


#' Generate HR values for a 1 unit increase in a variable at specified points of another variable
#'
#' This function does bla bla bla
#'
#' @param x numeric vector of var2 points to estimate
#' @param model model of class coxph. If data is NULL, the package expects to find the data in the model run with x = TRUE
#' @param data data used in the model. If absent, we will attempt to recover the data from the model object if x = TRUE
#' @param var1 variable that increases by 1 unit. If continuous, 1-SD increase is used
#' @param var2 variable to spline. x values belong to var2
#' @param ci calculate 95% CI?
#' @param ci.method confidence interval method. "delta" performs delta method. "bootstrap" performs bootstrapped CI (slower)
#' @param R number of bootstrap samples if ci.method = "bootstrap". Default 100
#' @param parallel can take values "no", "multicore", "snow" if ci.method = "bootstrap"
#' @param ... other parameters for boot
#' @return if ci = FALSE, a vector of estimate of length(x), if ci = TRUE a dataframe with 5 columns, initial values, HR, lower CI, upper CI and SE
#' @export
HRSpline <- function(x , model , data , var1 , var2
                        , ci=TRUE , conf = 0.95 , ci.method = "delta"
                        , ci.boot.method = "perc" , R = 100 , parallel = "multicore" , ...) {
  # Check correct class for model
  if( !any( c("coxph","cph") %in% class(model) ) ){
    stop("Interaction Cox model must be run with survival::coxph or rms::cph")
  }
  if(missing(data)){
    if(is.null(model$x) && ci && ci.method == "bootstrap"){
      # We need the data only for bootstrap CI
      stop("Missing data for bootstrap CI")
    } else {
      data <- model$x
    }
  }
  coefMod <- coef(model)
  # Different styles to name the interaction term between coxph and cph
  if("cph" %in% class(model)){
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
  HR <- unname(exp( var1Term + x*intTerm))
  if(ci){
    alpha <- qnorm( 1 - (1-conf)/2)
    if(ci.method=="exact"){
      stop("not implemented yet")
        SEgeneral <- sqrt( sum( c( diag(model$var)[c(mycoefWhich[var1] ,mycoefWhich[ myvars[3]])]
                                   , 2*model$var[mycoefWhich[var1] ,mycoefWhich[ myvars[3]]]) ) )
        SEvar1 <- sqrt(diag(model$var))
        # wrong
        # HR_L <- unname(exp( var1Term - alpha*SEvar1 + x*(intTerm- alpha*SEgeneral) ))
        # HR_U <- unname(exp( var1Term + alpha*SEvar1 + x*(intTerm+ alpha*SEgeneral) ))
        # also wrong
        # HR_L <- unname(exp( var1Term + x*(intTerm) - alpha*SEgeneral ))
        # HR_U <- unname(exp( var1Term + x*(intTerm) + alpha*SEgeneral ))
        # also wrong
        # HR_L <- unname(exp( var1Term + x*(intTerm- alpha*SEgeneral) ))
        # HR_U <- unname(exp( var1Term + x*(intTerm+ alpha*SEgeneral) ))
        HRci <- data.frame(Value = x
                           , HR = HR
                           , CI_L = HR_L
                           , CI_U = HR_U
                           , SE = SEgeneral)
        colnames(HRci) <- c("Value" , "HR" , "CI_L" , "CI_U" , "SE")
        return(HRci)
    } else if(ci.method == "delta"){
      # This creates a vector like x1 , x2 , x3 , x7 , x8
      # that tells you the position of the regressor as it appears in the model
      xNum <- paste0("x" , mycoefWhich)
      vcovMod <- vcov(model)
      HRci <- t(vapply( seq_len(length(x)) , function(i) {
        x_i <- x[i]
        myform <- paste0("~(", xNum[1] , " + " , xNum[3] , "*(" , x_i , "))")
        SE <- NULL
        try(SE<-msm::deltamethod(as.formula(myform), coefMod, vcovMod) , silent = TRUE)
        if(is.null(SE)){
          return(c(HR[i] , NA , NA , NA))
        }
        up<-exp(log(HR[i])+alpha*SE)
        lo<-exp(log(HR[i])-alpha*SE)
        c(HR[i] , lo , up , SE)
      } , numeric(4)))
      HRci <- cbind( Value = x , HRci)
      rownames(HRci) <- x
      colnames(HRci) <- c("Value" , "HR" , "CI_L" , "CI_U" , "SE")
      return(as.data.frame(HRci))
    } else if(ci.method == "bootstrap"){
      if(missing(parallel)){
        parallel <- "multicore"
      }
      if(missing(R)){
        R <- 100
      }
      myBoot <- boot::boot(data = data, statistic = .bootHRSpline
                           , x = x , model = model
                           , R = R , parallel = parallel
                           , var1 = var1 , var2 = var2 )
      HRci <- t(vapply( seq_len(length(x)) , function(idx) {
        bci <- boot::boot.ci(boot.out = myBoot,  index = idx
                             , type = ci.boot.method , conf = conf)
        if(ci.boot.method %in% "norm"){
          c(bci$t0 , bci$normal[2] , bci$normal[3])
        } else {
          c(bci$t0 , bci[[4]][4] , bci[[4]][5])
        }
      } , numeric(3)))
      HRci <- cbind(x , HRci)
      colnames(HRci) <- c("Value" , "HR" , "CI_L" , "CI_U")
      rownames(HRci) <- x
      return(as.data.frame(HRci))
    } else {
      stop("Only delta and bootstrap are valid CI methods")
    }
  } else {
    return(as.data.frame(HR))
  }
}
