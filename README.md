<!-- badges: start -->
[![R-CMD-check](https://github.com/gmelloni/interactionHR/workflows/R-CMD-check/badge.svg)](https://github.com/gmelloni/interactionHR/actions)
<!-- badges: end -->

# interactionHR
#### A tool to calculate and plot Hazard Ratios after a Cox model in which an interaction between the main predictor and a continuous covariate has been specified.
#### Version 0.1 (October 27, 2021)
---

### Description
`interactionHR` estimates the Hazard Ratio for one covariate of interest (binary or continuous) over levels of a second continuous 
  covariate with which an interaction has been specified in a Cox model with `cph`. In particular, `interactionHR` allows for
  basic interaction assessment (i.e. log-linear interaction model where a product term between the two predictors is included) 
  as well as settings where the second covariate is flexibly modeled with restricted cubic splines. Confidence intervals for 
  the predicted Hazard Ratios can be calculated with either bootstrap or the delta method. Lastly, `interactionHR`
  produces a plot of the hazard ratio over levels of the other covariate.

### Installation
To install the last version of `interactionHR` from GitHub, type
```
devtools::install_github("https://github.com/gmelloni/interactionHR.git")
library(interactionHR)
```
from within a web-aware R.

### Vignette
For an introduction to `interactionHR` please refer to this [vignette](https://raw.githack.com/gmelloni/interactionHR/main/inst/data/vignette.html)

### Authors
Giorgio Melloni, Sabina Murphy, Andrea Bellavia

TIMI study group, Brigham and Womens Hospital / Harvard Medical School
