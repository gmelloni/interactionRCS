# # Add generic
# plot <- function (x, ...) {
#   UseMethod("plot", x)
# }

#' Plot the result of HR, OR or linear estimates
#'
#' Create a spline var2 by 1 unit increase of var1
#'
#' @param x data.frame calculated using any of the function of this package
#' @param xlab xlab name
#' @param ylab ylab name. Default is the estimate column name if log=FALSE otherwise Estimate(log scale)
#' @param main plot title
#' @param log if TRUE, plot the estimate in log scale
#' @param line1 if TRUE, plot horizontal line on 1 or 0 (if log=TRUE)
#' @param color line color. Default dodgerblue
#' @return simple splined plot of estimates of var1 at var2 values
#' @examples
#' library(rms)
#' library(survival)
#' data(cancer)
#' myformula <- Surv(time, status) ~ ph.karno + ph.ecog + rcs(age, 3)*sex
#' model <- cph(myformula , data = lung )
#' myHR <- rcsHR( var2values = 40:80
#'                , model = model , data = lung , var1 ="sex", var2="age"
#'                , ci=TRUE , conf = 0.95 , ci.method = "delta")
#' plotINT(myHR , ylab = "HR of male VS female" , xlab = "Age")
#' @importFrom pspline sm.spline
#' @importFrom graphics par lines abline plot
#' @export
plotINT <- function(x , xlab = "" , main = "" , log = FALSE
                  , ylab = NULL ,line1 = TRUE
                  , color="dodgerblue"){
  if(is.null(ylab)){
    ylab <- colnames(x)[2]
    if(ylab == "LIN"){
      ylab <- "Linear Estimate"
    }
    if(log){
      ylab <- paste0(ylab , " (log scale)")
    }
  }
  if("CI_L" %in% colnames(x) && "CI_U" %in% colnames(x)){
    ylim <- c(min(x$CI_L , na.rm = TRUE),max(x$CI_U , na.rm = TRUE))
  } else {
    ylim <- c(min(x[ , 2] , na.rm = TRUE),max(x[ , 2] , na.rm = TRUE))
  }
  plot( x$Value , x[ , 2]
        , type="n" , xlab = xlab
        , ylab = ylab
        , log = if(log) "y" else ""
        , ylim = ylim
        , main = main)
  if(line1) {
    if(log)
      abline(h=0 , lty = 3 , lwd = 1 , col = "black")
    else
      abline(h=1 , lty = 3 , lwd = 1 , col = "black")
  }
  lines( pspline::sm.spline(x$Value , x[ , 2]) , col = color , lty = 1 , lwd = 3 )
  if("CI_L" %in% colnames(x) && "CI_U" %in% colnames(x)){
    lines( pspline::sm.spline(x$Value , x$CI_L) , col = "darkgray" , lty = 2 , lwd = 2 )
    lines( pspline::sm.spline(x$Value , x$CI_U) , col = "darkgray" , lty = 2 , lwd = 2 )
  }
}
