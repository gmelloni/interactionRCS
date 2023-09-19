<!-- badges: start -->
[![R-CMD-check](https://github.com/gmelloni/interactionRCS/workflows/R-CMD-check/badge.svg)](https://github.com/gmelloni/interactionRCS/actions)
<!-- badges: end -->

# interactionRCS
#### A tool to calculate and plot Hazard Ratios, Odds Ratios or linear estimates in a simple or restricted cubic splined interaction model
#### Version 1.1 (February 25, 2022)
---

### Description
`interactionRCS` facilitates interpretation and presentation of results from a regression model (linear, logistic, Cox) where an interaction between the main predictor of interest X (binary or continuous) and another continuous covariate Z has been specified. In particular, `interactionRCS` allows for
  basic interaction assessment (i.e. log-linear/linear interaction models where a product term between the two predictors is included) 
  as well as settings where the second covariate is flexibly modeled with restricted cubic splines. Confidence intervals for 
  the predicted effect measures (beta, OR, HR) can be calculated with either bootstrap or the delta method. Lastly, `interactionRCS`
  produces a plot of the effect measure over levels of the other covariate.

### Installation
To install the latest version of `interactionRCS`, type the following lines in a web-aware R environment.

```
if(!"devtools" %in% rownames(installed.packages())){
  install.packages("devtools")
}
devtools::install_github("https://github.com/gmelloni/interactionRCS.git")
# or alternative devtools::install_git("https://github.com/gmelloni/interactionRCS.git")
library(interactionRCS)
```

### Usage 
After estimating a regression model (linear, logistic, Cox) such as `model<-glm(y~ ...)` estimate and plot interactions with:

```
int<-estINT(model=model, ...)

plotINT(int, ...)
```

For a detailed introduction to `interactionRCS` and code examples please refer to this [vignette](https://raw.githack.com/gmelloni/interactionRCS/main/inst/extdata/vignette.html)

### Authors
Giorgio Melloni, Hong Xiong, Andrea Bellavia

TIMI study group, Department of Cardiovascular Medicine, Brigham and Womens Hospital / Harvard Medical School
