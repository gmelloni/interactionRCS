#' Plot the result of HRSpline functions
#'
#' This function does bla bla bla
#'
#' @param df data.frame calculated using HRcubSpline
#' @param log if TRUE, plot HR in log scale
#' @param xlab xlab name
#' @param main plot title
#' @return simple splined plot of your var2 values by HR of
#' @export
plot.HRSpline <- function(df , xlab = "" , main = "" , log = FALSE){
  plot( df$Value , df$HR
        , type="n" , xlab = xlab
        , ylab = if(log) "HR (log scale)" else "HR"
        , log = if(log) "y" else ""
        , ylim = c(min(df$CI_L , na.rm = TRUE),max(df$CI_U , na.rm = TRUE))
        , main = "")
  lines( pspline::sm.spline(df$Value , df$HR) , col = "dodgerblue" , lty = 1 , lwd = 3 )
  lines( pspline::sm.spline(df$Value , df$CI_L) , col = "grey" , lty = 2 , lwd = 2 )
  lines( pspline::sm.spline(df$Value , df$CI_U) , col = "grey" , lty = 2 , lwd = 2 )
}
