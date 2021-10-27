# Add generic
plot <- function (x, ...) {
  UseMethod("plot", x)
}

#' Plot the result of an object of class HR
#'
#' Create a spline var2 by HR of var1
#'
#' @param df data.frame calculated using HRcubSpline
#' @param xlab xlab name
#' @param main plot title
#' @param log if TRUE, plot HR in log scale
#' @param line1 if TRUE, plot horizontal line on HR = 1 or log(HR) = 0
#' @return simple splined plot of your var2 values by HR
#' @examples
#' library(survival)
#' data(cancer)
#' myformula <- Surv(time, status) ~ ph.karno + ph.ecog + rcs(age, 3)*sex
#' model <- cph(myformula , data = lung )
#' myHR <- HRcubSpline( var2values = 40:80
#'                      , model = model , data = lung2 , var1 ="sex", var2="age" , units=1
#'                      , center = 0, ci=TRUE , conf = 0.95 , ci.method = "delta")
#' plot(myHR , ylab = "HR of male VS female" , xlab = "Age")
#' @export
plot.HR <- function(x , xlab = "" , main = "" , log = FALSE
                    , ylab = if(log) "log(HR)" else "HR" ,line1 = TRUE){
  plot( x$Value , x$HR
        , type="n" , xlab = xlab
        , ylab = if(log) "HR (log scale)" else "HR"
        , log = if(log) "y" else ""
        , ylim = c(min(x$CI_L , na.rm = TRUE),max(x$CI_U , na.rm = TRUE))
        , main = main)
  if(line1) { if(log) abline(h=0 , lty = 3 , lwd = 2 , col = "gray") else abline(h=1 , lty = 3 , lwd = 2 , col = "gray")}
  lines( pspline::sm.spline(x$Value , x$HR) , col = "dodgerblue" , lty = 1 , lwd = 3 )
  if("CI_L" %in% colnames(x) && "CI_U" %in% colnames(x)){
    lines( pspline::sm.spline(x$Value , x$CI_L) , col = "black" , lty = 2 , lwd = 2 )
    lines( pspline::sm.spline(x$Value , x$CI_U) , col = "black" , lty = 2 , lwd = 2 )
  }
}
