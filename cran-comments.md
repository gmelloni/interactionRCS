## Version: 0.1.3

This is a bug fix update and a feature update

* Replaced use of the removed `PimaIndiansDiabetes` data set (mlbench) with
  `SynthDiabetes`, a synthetic replacement with the same variables, in all
  examples. Added a version requirement `mlbench (>= 2.1-11)` to Suggests,
  since that is the first version of mlbench containing `SynthDiabetes`.
* Added support for Poisson (and quasipoisson) regression in `loglinOR()`,
  `rcsOR()` and the `intEST()` dispatcher. Poisson models return rate ratios
  (RR) instead of odds ratios.
* Fixed plotINT when log=TRUE not showing abline on 1

## Version: 0.1.2

This is a bug fix update

* Fixed dependency to pryr (soon deprecated) by switching to rlang 

## Version: 0.1.1

This is a feature update

* Added support for multi knots rcs

## Version: 0.1.0

* This is the initial release of the package

## R CMD check results

0 errors | 0 warnings | 0 note

## revdepcheck results

No reverse dependencies. We saw 0 new problems. We failed to check 0 packages.
