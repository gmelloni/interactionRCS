---
title: 'interactionRCS: A tool to calculate and plot Hazard Ratios after a Cox model in which an interaction between the main predictor and a continuous covariate has been specified.'
tags:
  - R
  - biostatistics
  - survival analysis
  - interactions
  - splines modeling
authors:
  - name: Giorgio Melloni
    orcid: 0000-0000-0000-0000
    affiliation: "1, 2"
  - name: Andrea Bellavia
    affiliation: "1, 2"
affiliations:
 - name: TIMI Study Group, Brigham and Womens Hospital
   index: 1
 - name: Department of Cardiovascular Medicine, Harvard Medical School
   index: 2
date: 25 January 2022
bibliography: paper.bib
---


# Summary

The `interactionRCS` package is designed to facilitate interpretation and presentation of results from a COX regression model where an interaction between the main predictor of interest $X$ (binary or continuous) and another continuous covariate $Z$ has been specified. Specifically, the package will provide point estimates of the HR for the main predictor $X$ over levels of $Z$, allowing for settings where $Z$ is flexibly modeled with restricted cubic splines, and provide a graphical display of this interaction. Two methods for deriving and plotting confidence intervals are also implemented, including the delta method and bootstrap.

# Introduction


# Mathematics

`interactionRCS` requires results from a Cox model where an interaction between a main predictor (binary or continuous) $X$ and a continuous predictor $Z$ has been specified. This interaction can be included as a simple product term between the 2 predictors, or by flexibly modeling $Z$ with restricted cubic splines. For both interaction settings, the main exposure of interest $X$ has to either be binary or continuous.  

## Log-linear interaction model

A basic Cox model with 2 predictors and their interaction takes the form:

$h(t|x,z)=h_0\cdot\exp(\beta_1x+\beta_2z+\beta_3x\cdot z)$  

After estimation of the model, we want to predict the HR for $X$ over levels of $Z$. This is given by


$HR_{10}=\frac{h(t|x=1,z)}{h(t|x=0,z)}=\frac{h_0\cdot\exp(\beta_1x+\beta_2z+\beta_3x\cdot z)}{h_0\cdot\exp(\beta_1x+\beta_2z+\beta_3x\cdot z)}=\frac{\exp(\beta_1+\beta_2z+\beta_3z)}{\exp(\beta_2z)}$  ,

which will be plotted against $Z$

To estimate $95\%$ confidence intervals $SE(HR_{10})$ is required. This is obtained by focusing on $\log(HR_{10})$ and calculating lower and upper bounds for this quantity, which are then exponentiated.

$SE(\log(HR_{10}))=SE(\log(\frac{\exp(\beta_1+\beta_2z+\beta_3z)}{\exp(\beta_2z)}))=SE(\log(\exp(\beta_1+\beta_2z+\beta_3z)-\log(\exp(\beta_2z)))=SE(\beta_1+\beta_3z)$

This is calculated by using the delta method. Upper and lower bounds are then derived with $\exp(\log(HR_{10})\pm1.96\cdot SE(log(HR_{10})))$

## Restricted cubic splines interaction model

`interactionRCS` allows the continuous covariate $Z$ to be flexibly modeled with restricted cubic splines, with 3 knots ($k_1, k_2, k_3$). The interaction model, in this setting, takes the form:

$h(t|x,z,c)=h_0\cdot\exp(\beta_1x+sp(z)+sp_2(z\cdot x))$  

where

$sp(z)=\alpha_1z+\alpha_2\cdot[\frac{(z-k_1)^3_+-(z-k_2)^3_+\cdot\frac{ k_2-k_1}{ k_3-k_2} +(z-k_3)^3\cdot\frac{k_3-k_1}{k_3-k_2}}{(k_3-k_1)^2}]$

and 

$sp_2(z\cdot x)=\gamma_1xz+\gamma_2x\cdot[\frac{(z-k_1)^3_+-(z-k_2)^3_+\cdot\frac{ k_3-k_1}{ k_3-k_2} +(z-k_3)^3\cdot\frac{k_2-k_1}{k_3-k_2}}{(k_3-k_1)^2}]$


After estimation of the model, we want to predict the HR for $X$ over levels of $Z$. This is given by


$HR_{10}=\frac{h(t|X=x_1,Z)}{h(t|X=x_2,Z)}=\frac{h_0\cdot\exp(\beta_1x_1+sp(z)+sp_2(z)\cdot x_1)}{h_0\cdot\exp(\beta_1x_2+sp(z)+sp_2(z)\cdot x_2)}=\frac{h(t|X=x_1,Z)}{h(t|X=x_2,Z)}=\frac{\exp(\beta_1+sp(z)+sp_2(z))}{\exp(sp(z))}$  ,

which will be plotted against $Z$

Similarly to the previous situation, the $SE$ can be derived by using the delta method to calculate $SE(\beta_1+sp_2(z))$


# How to use the interactionRCS package 

Functions within the `interactionRCS` package require that a Cox model has already been estimated and model results be provided as an object. Because of its flexibility in dealing with restricted cubic splines, `interactionRCS` requires the Cox model to be fitted with the `cph` function of  `rms` package. 

The main functions of `interactionRCS` are `rcsHR` and `loglinHR`. The first one will provide point estimates and confidence intervals for the  HRs of $X$ when $Z$ is modeled with restricted cubic splines (specifically using 3 knots). The second function will instead provide HRs of $X$ in the setting where $Z$ is included in the model as a continuous covariate thus assuming a log-liner effect additive interaction on the log-linear scale. For both functions, the following options must be specified: 

* `model`: the `cph` model previously run (linlogHR accepts `coxph` objects too)
* `var1`: the name of the main predictor of interest ($X$)
* `var2`: the name of the continuous predictor interacting with `var1` ($Z$)
* `var2values`: the values of $var2$ for which the HR of $var1$ should be calculated

Additional options include:

* `data`: the same dataset used for fitting the Cox model (only used for bootstrap CIs). If data=NULL, we will search for model$x
* `ci` (default TRUE) : whether a confidence interval for each HR should also be provided
* `ci.method`: either `"delta"` or `"bootstrap"`. Default `"delta"`
* `ci.boot.method` (default= "percentile" - only if `method="bootstrap"`) : see `boot.ci` type parameter 
* `R` (default=100 - only if `method="bootstrap"`) : number of bootstrap iterations
* `parallel` (default "multicore" - only if `method="bootstrap"`): see `boot.ci` reference 

A third function `plotHR` is also implemented to provide a graphical display of the results and only require the `rcsHR` or `loglinHR` results as object. 


# Illustrative examples

The first example is based on a study on drug relapse among 575 patients enrolled in a clinical trial of residential treatment for drug abuse. The main exposure of interest is the binary indicator of assigned treatment (0/1) and a treatment*age interaction is specified.

The following code provides an estimate of the treatment HR at different ages, when age is modeled with restricted cubic splines. The model is further adjusted for tumor site, race, and previous use of IV drug. It is recommended to check the distribution of $Z$ (here, age) to define a realistic range of `var2values`. Finally, the code is replicated by using the delta method or bootstrap for obtaining confidence intervals. Note that when `ci.method = "bootstrap"` is specified, additional options can be specified. Figure in \autoref{fig:example}

![Caption for example figure.\label{fig:example}](figure1.png)


# Final remarks

# References

