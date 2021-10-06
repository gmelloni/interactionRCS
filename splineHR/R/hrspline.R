# .bootHRSpline <- function(data, idx , x , model , var1 , var2){
#   df <- data[idx, ]
#   # Grab the call and recreate it with new data
#   mycall <- model$call
#   mycall2 <- pryr::modify_call(mycall , list(data=df))
#   # mymodel <- cph(model$sformula, data=df)
#   mymodel <- eval(mycall2)
#   coefMod <- coef(mymodel)
#   ##### SHOULD THE NODES BE GLOBAL OR RUN SPECIFIC? HERE SPECIFIC
#   # if(missing(k)){
#   k <- attributes(rcs(df[[var2]], 3))$parms
#   # }
#   # k <- attributes(rcs(df[[var2]], 3))$parms # like this knots are unique for every run
#   myvars <- intersect( c(var1 , var2
#                          , paste0(var2 , "'")
#                          , paste(var1 , var2 , sep = " * ")
#                          , paste(var2 , var1 , sep = " * ")
#                          , paste(var1 , paste0(var2 , "'") , sep = " * ")
#                          , paste(paste0(var2 , "'"),var1 , sep = " * "))
#                        , names(coefMod))
#   mycoef <- coefMod[  myvars ]
#   a <- mycoef[ c(var2 , paste0(var2,"'"))]
#   b <- mycoef[ var1 ]
#   l <- mycoef[ setdiff(myvars , c(var2 , paste0(var2,"'") , var1)) ]
#   numer <- vapply(x , function(i) {
#     max(i - k[1],0)^3  - (max(i - k[2],0)^3)*k3k1 + (max(i - k[3],0)^3)*k2k1
#   } , numeric(1))
#   denom <- (k[3] - k[1])^2
#   numDem <- numer/denom
#   sp1 <- vapply(x , function(i) a[1]*i , numeric(1)) + a[2]*numDem
#   sp2 <- vapply(x , function(i) l[1]*i , numeric(1)) + l[2]*numDem
#   HR <- unname(exp( b + sp1 + sp2)/exp(sp1))
# }

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
#' @param ci.method confidence interval method. "delta" performs delta method (fast but inaccurate). "boot" performs bootstrapped CI (slower)
#' @param R number of bootstrap samples if ci.method = "bootstrap". Default 100
#' @param parallel can take values "no" , "multicore" , "snow" if ci.method = "bootstrap"
#' @param ... other parameters for boot
#' @return if ci = FALSE, a vector of estimate of length(x), if ci = TRUE a matrix with 4 columns, initial values, HR, lower CI and upper CI
#' @export
# HRSpline <- function(x , model , data , var1 , var2
#                         , ci=TRUE , conf = 0.95 , ci.method = "delta"
#                         , ci.boot.method = "perc" , R = 100 , parallel = "multicore" , ...) {
#   # Check correct class for model
#   if( !all( c("coxph") %in% class(model) ) ){
#     stop("Interaction Cox model must be run with survival::coxph")
#   }
#
#
#   if(missing(data)){
#     if(is.null(model$x)){
#       stop("Missing data")
#     } else {
#       data <- model$x
#     }
#   }
#   coefMod <- coef(model)
#   k <- attributes(rms::rcs(data[[var2]], 3))$parms
#   k2k1 <- (k[2] - k[1])/(k[3] - k[2])
#   k3k1 <- (k[3] - k[1])/(k[3] - k[2])
#   # Extract parameters
#   # beta1 <- coefMod["score.sd[1]"]
#   # alph <- coefMod[c("age" , "age'")]
#   # lamb <- coefMod[c("score.sd[1] * age" , "score.sd[1] * age'")]
#   myvars <- intersect( c(var1 , var2
#                          , paste0(var2 , "'")
#                          , paste(var1 , var2 , sep = " * ")
#                          , paste(var2 , var1 , sep = " * ")
#                          , paste(var1 , paste0(var2 , "'") , sep = " * ")
#                          , paste(paste0(var2 , "'"),var1 , sep = " * "))
#                        , names(coefMod))
#   mycoef <- coefMod[  myvars ]
#   a <- mycoef[ c(var2 , paste0(var2,"'"))]
#   b <- mycoef[ var1 ]
#   l <- mycoef[ setdiff(myvars , c(var2 , paste0(var2,"'") , var1)) ]
#   numer <- vapply(x , function(i) {
#     max(i - k[1],0)^3  - (max(i - k[2],0)^3)*k3k1 + (max(i - k[3],0)^3)*k2k1
#   } , numeric(1))
#   denom <- (k[3] - k[1])^2
#   numDem <- numer/denom
#   sp1 <- vapply(x , function(i) a[1]*i , numeric(1)) + a[2]*numDem
#   sp2 <- vapply(x , function(i) l[1]*i , numeric(1)) + l[2]*numDem
#   HR <- unname(exp( b + sp1 + sp2)/exp(sp1))
#   if(ci){
#     if(ci.method == "delta"){
#       vcovMod <- vcov(model)
#       HRci <- t(vapply( seq_len(x) , function(i) {
#         x_i <- x[i]
#         numDem_i <- numDem[i]
#         myform <- paste0("~(x1 + x2*" , x_i
#                          , " + x3*" , numDem_i
#                          , " + x10*" , x_i
#                          , " + x11*" , numDem_i
#                          , ")/(x2*" , x_i
#                          , " + x3*" , numDem_i , ")")
#         SE<-msm::deltamethod(as.formula(myform), coefMod, vcovMod)
#         up<-exp(log(HR[i])+qnorm( 1 - (1-conf)/2)*SE)
#         lo<-exp(log(HR[i])-qnorm( 1 - (1-conf)/2)*SE)
#         c(HR[i] , lo , up , SE)
#       } , numeric(4)))
#       HRci <- cbind( Value = x , HRci)
#       rownames(HRci) <- x
#       colnames(HRci) <- c("Value" , "HR" , "CI_L" , "CI_U" , "SE")
#       return(as.data.frame(HRci))
#     } else if(ci.method == "bootstrap"){
#       if(missing(parallel)){
#         parallel <- "multicore"
#       }
#       if(missing(R)){
#         R <- 100
#       }
#       myBoot <- boot::boot(data = data, statistic = .bootHRcubSpline, x = x , model = model
#                            , R = R , parallel = parallel, var1 = var1 , var2 = var2 )
#       HRci <- t(vapply( seq_len(length(x)) , function(idx) {
#         bci <- boot::boot.ci(boot.out = myBoot,  index = idx , type = ci.boot.method , conf = conf)
#         if(ci.boot.method %in% "norm"){
#           c(bci$t0 , bci$normal[2] , bci$normal[3])
#         } else {
#           c(bci$t0 , bci[[4]][4] , bci[[4]][5])
#         }
#       } , numeric(3)))
#       HRci <- cbind(x , HRci)
#       colnames(HRci) <- c("Value" , "HR" , "CI_L" , "CI_U")
#       rownames(HRci) <- x
#       return(as.data.frame(HRci))
#     } else {
#       stop("Only delta and bootstrap are valid CI methods")
#     }
#   } else {
#     return(as.data.frame(HR))
#   }
# }
