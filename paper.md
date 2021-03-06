---
title: 'interactionRCS: An R package to display flexible interactions from common statistical modeling'
tags:
  - R
  - biostatistics
  - interactions
  - splines modeling
authors:
  - name: Giorgio Melloni
    orcid: 0000-0001-6371-1334
    affiliation: "1, 2"
  - name: Andrea Bellavia
    orcid: 0000-0003-4988-4532
    affiliation: "1, 2"
affiliations:
 - name: TIMI Study Group, Brigham and Womens Hospital
   index: 1
 - name: Division of Cardiovascular Medicine, Harvard Medical School
   index: 2
date: 28 February 2022
bibliography: paper.bib
---

# Summary

Interaction analysis is a critical part of clinical and epidemiological research that allows evaluating how the effects of a given covariate (e.g. a novel treatment, a potential risk factor) on a health outcome changes over levels of other patients characteristics. When assessing interactions with a continuous predictor such as age, the complexity of possible non-linear interactions should be taken into account using flexible tools such as restricted cubic splines (RCS). This flexibility comes with additional complexity, and estimating and presenting interactions based on regression models estimates in this context is not straightforward. With the `interactionRCS` package that we introduce here, we have developed a simple and user-friendly tool that allows estimating and graphically presenting complex interaction effects after estimation of various regression models (linear, logistic, and Cox) in the R statistical software.

# Statement of need

Statistical and biostatistical research are often interested in situations where the main effect of a predictor of interest varies over levels of another predictor. In clinical studies this common situation arises, for example, when the effect of a new treatment in a randomized clinical trial (RCT) changes over levels of participants characteristics such as sex or age. In this context, the joint effect of the two predictors of interest (treatment and age, for example) on a given response is captured by their independent effects and by an interaction effect that quantifies the additional change in the response when both predictors are operating. A thorough description of interaction effects and their interpretation and relevance can be found in @vanderweele2014tutorial. 

In practice, interaction assessment is commonly conducted by including a product term between the two predictors of interest in a statistical model. For instance, in a RCT investigating the effect of a novel treatment on cardiovascular mortality researchers might fit a Cox PH regression that includes a first term for treatment, a second term for the second predictor, and a third term for the interaction between treatment and the other predictor. If the second predictor was a binary variable (e.g. participants sex), the treatment effect would be summarized with two different parameters, one for each value of sex. On the other hand, with a continuous predictor such as age, the effect of the main predictor would require some graphical presentation. Figure 1 provides an example of such graphical presentation, displaying the effect of a binary predictor, in the form of a hazard ratio (HR), over levels of participants age.

![Effect of a binary predictor over levels of participants age, assuming a linear interaciton on the logarithm of the HR.\label{fig:example}](figure2.png){width=80%}

When dealing with continuous predictors such as age, researchers need to take into account the linearity assumptions made by the different statistical models. In particular, the most common statistical approaches used in clinical research share an assumption of linearity either on the response (linear regression), on the logit of the event probability (logistic regression), or on the logarithm of the hazard ratio (Cox PH regression). In the context of interaction analysis, this not only implies that the effect of each continuous predictor will be linear on the underlying scale, but also that the interaction term will be linear. Using the example presented in Figure 1, the treatment effect changes over age levels in a log-linear fashion. This assumption is often not met in real data, and different approaches to relax linearity exist. Among them, modeling continuous predictors with restricted cubic splines represent one of the most flexible approaches, currently recommended by several guidelines and researchers.(@greenland1995dose, @von2007strengthening, @durrleman1989flexible) An introduction to restricted cubic splines can be found in Chapter 2.4 of @harrell2017regression. On our github page we provide a detailed presentation of the mathematical formulation for linear, logistic, and Cox regression models, when an interaction is included and modeled with restricted cubic splines. By using this approach any linearity assumption is relaxed and the interaction effect can be graphically displayed with a smooth and flexible function that more accurately captures how the main effect changes over levels of the other predictors. Figure 2, for example, presents the treatment/age interaction previously described, now flexibly modeled with this approach.

![Effect of a binary predictor over levels of participants age, flexibly modeling the interaction with restricted cubic splines.\label{fig:example}](figure1.png){width=80%}

# Functionalities

To the best of our knowledge, simple procedures to estimate effects and obtain such graphical presentation of interaction effects, while possibly allowing for restricted cubic splines modeling, are unavailable in all major statistical software. The `interactionRCS` package that we have developed includes a set of functions that allow estimating effects in the presence of interactions and deriving these graphical presentations after fitting a linear, logistic, or a Cox PH regression model in `R`. Specifically, the package will provide point estimates of the effect of a predictor $X$ (binary or continuous) over levels of another continuous covariate $Z$, allowing for settings where $Z$ is flexibly modeled with restricted cubic splines, and provide a graphical display of this interaction. Two methods for deriving and plotting confidence intervals are also implemented, including the delta method and bootstrap.

Functions within the `interactionRCS` package require that a regression model has already been estimated and model results be provided as an object. We have implemented several functions to estimate effects in various situations of interest (continuous, binary, and time-to-event outcomes, as well as under different approaches to model $Z$), which can be called by running a unique user-friendly function `intEST`, and an additional function `plotINT` to generate plots such as the ones presented in Figures 1 and 2. Additional details on how the functions are implemented and a theoretical background on how the effects are calculated are presented in the online vignette available on the github page.

# References

