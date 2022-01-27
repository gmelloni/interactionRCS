<!-- badges: start -->
[![R-CMD-check](https://github.com/gmelloni/interactionRCS/workflows/R-CMD-check/badge.svg)](https://github.com/gmelloni/interactionRCS/actions)
<!-- badges: end -->

# interactionRCS
#### A tool to calculate and plot Hazard Ratios after a Cox model in which an interaction between the main predictor and a continuous covariate has been specified.
#### Version 0.1 (October 27, 2021)
---

### Description
`interactionRCS` estimates the Hazard Ratio for one covariate of interest (binary or continuous) over levels of a second continuous 
  covariate with which an interaction has been specified, after ftting a Cox model with `cph`. In particular, `interactionRCS` allows for
  basic interaction assessment (i.e. log-linear interaction model where a product term between the two predictors is included) 
  as well as settings where the second covariate is flexibly modeled with restricted cubic splines. Confidence intervals for 
  the predicted Hazard Ratios can be calculated with either bootstrap or the delta method. Lastly, `interactionRCS`
  produces a plot of the hazard ratio over levels of the other covariate.

### Installation
To install the latest version of `interactionRCS`, type the following lines in a web-aware R environment.

```
if(!"devtools" %in% rownames(installed.packages())){
  install.packages("devtools")
}
devtools::install_github("https://github.com/gmelloni/interactionRCS.git")
library(interactionRCS)
```


### Vignette
For a detailed introduction to `interactionRCS` and code examples please refer to this [vignette](https://raw.githack.com/gmelloni/interactionRCS/main/inst/extdata/vignette.html)

### Authors
Giorgio Melloni, Andrea Bellavia

TIMI study group, Department of Cardiovascular Medicine, Brigham and Womens Hospital / Harvard Medical School
