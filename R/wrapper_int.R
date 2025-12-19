#' Returns the estimates of for an unspecified interaction model
#'
#' This function is a dispatcher that generate OR, HR or linear estimates values for a simple or restricted cubic spline
#' interaction model from a logistic, Cox or linear regression
#'
#' @param var2values numeric vector of var2 points to estimate
#' @param model model of class cph, coxph, lrm, glm or Glm. If data is NULL, the function expects to find the data in model$x.
#' @param data data used in the model. If absent, we will attempt to recover the data from the model. Only used for bootstrap and glm class models
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
#' # Logistic model predicting diabetes over BMI, age and glucose
#' myformula <- diabetes ~ mass + age * rcs( glucose , 3 )
#' model <- lrm(myformula , data = PimaIndiansDiabetes )
#' intEST( var2values = 20:80
#'        , model = model , data = PimaIndiansDiabetes , var1 ="age", var2="glucose"
#'        , ci=TRUE , conf = 0.95 , ci.method = "delta")
#' # Linear model predicting BMI over diabetes, age and glucose
#' myformula2 <- mass ~ diabetes + age * rcs( glucose , 3 )
#' model2 <- glm(myformula2 , data = PimaIndiansDiabetes , family = "gaussian")
#' intEST( var2values = 20:80
#'        , model = model2 , data = PimaIndiansDiabetes , var1 ="age", var2="glucose"
#'        , ci=TRUE , conf = 0.95 , ci.method = "delta")

#' @return if ci = FALSE, a dataframe with initial values and OR/HR/linear estimates
#' , if ci = TRUE a dataframe with 5 columns, initial values, OR/HR/linear estimates, lower CI, upper CI and SE
#' @importFrom rms lrm rcs cph
#' @importFrom survival coxph
#' @importFrom rlang call_modify
#' @importFrom msm deltamethod
#' @importFrom boot boot boot.ci
#' @importFrom stats vcov coef as.formula qnorm sd glm
#' @export
intEST <- function(var2values , model , data , var1 , var2
                   , ci=TRUE , conf = 0.95 , ci.method = "delta"
                   , ci.boot.method = "perc" , R = 100 , parallel = "multicore" , ...){
  # argg <- c(as.list(environment()), list(...))
  argg <- match.call(expand.dots = TRUE)
  modelClass <- class(model)
  isRCS <- any( grepl(paste0("rcs(",var2) , attributes(model$terms)$term.labels , fixed = TRUE))
  if(any(c("cph" , "coxph") %in% modelClass)){
    if(isRCS){
      argg[[1]] <- as.name("rcsHR")
    } else {
      argg[[1]] <- as.name("loglinHR")
    }
  } else if("lrm" %in% modelClass){
    if(isRCS){
      argg[[1]] <- as.name("rcsOR")
    } else {
      argg[[1]] <- as.name("loglinOR")
    }
  } else if(any(c("Glm" , "glm") %in% modelClass)){
    family <- model$family$family
    if(family %in% c("gaussian","quasi")){
      if(isRCS){
        argg[[1]] <- as.name("rcsLIN")
      } else {
        argg[[1]] <- as.name("linLIN")
      }
    } else if(family == "binomial"){
      if(isRCS){
        argg[[1]] <- as.name("rcsOR")
      } else {
        argg[[1]] <- as.name("loglinOR")
      }
    } else {
      stop("glm or Glm object not of family gaussian nor binomial")
    }

  } else {
    stop("Object of unrecognized class. coxph, cph, glm or Glm are all accepted")
  }
  return(eval.parent(argg))
}
