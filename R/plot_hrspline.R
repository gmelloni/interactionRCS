# # Add generic
# plot <- function (x, ...) {
#   UseMethod("plot", x)
# }

#' Plot the result of an object of class HR
#'
#' Create a spline var2 by HR of var1
#'
#' @param x data.frame calculated using HRcubSpline
#' @param xlab xlab name
#' @param ylab ylab name. Default HR if log=FALSE otherwise HR(log scale)
#' @param main plot title
#' @param log if TRUE, plot HR in log scale
#' @param line1 if TRUE, plot horizontal line on HR = 1 or log(HR) = 0
#' @param color HR line color. Default dodgerblue
#' @return simple splined plot of your var2 values by HR
#' @examples
#' library(rms)
#' library(survival)
#' data(cancer)
#' myformula <- Surv(time, status) ~ ph.karno + ph.ecog + rcs(age, 3)*sex
#' model <- cph(myformula , data = lung )
#' myHR <- rcsHR( var2values = 40:80
#'                      , model = model , data = lung , var1 ="sex", var2="age"
#'                      , ci=TRUE , conf = 0.95 , ci.method = "delta")
#' plotHR(myHR , ylab = "HR of male VS female" , xlab = "Age")
#' @importFrom pspline sm.spline
#' @importFrom graphics par lines abline plot
#' @export
plotHR <- function(x , xlab = "" , main = "" , log = FALSE
                  , ylab = if(log) "log(HR)" else "HR" ,line1 = TRUE
                  , color="dodgerblue"){
  plot( x$Value , x$HR
        , type="n" , xlab = xlab
        , ylab = if(log) "HR (log scale)" else "HR"
        , log = if(log) "y" else ""
        , ylim = c(min(x$CI_L , na.rm = TRUE),max(x$CI_U , na.rm = TRUE))
        , main = main)
  if(line1) {
    if(log)
      abline(h=0 , lty = 3 , lwd = 1 , col = "black")
    else
      abline(h=1 , lty = 3 , lwd = 1 , col = "black")
  }
  lines( pspline::sm.spline(x$Value , x$HR) , col = color , lty = 1 , lwd = 3 )
  if("CI_L" %in% colnames(x) && "CI_U" %in% colnames(x)){
    lines( pspline::sm.spline(x$Value , x$CI_L) , col = "darkgray" , lty = 2 , lwd = 2 )
    lines( pspline::sm.spline(x$Value , x$CI_U) , col = "darkgray" , lty = 2 , lwd = 2 )
  }
}
