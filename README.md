Assurance R-package
===================


This R package library simplifies some calculations for the success of a
clinical trial given data from some initial trial.  The approach
throughout is to produce a prior for effect size, use the initial
trial data to produce a posterior distribution for effect size, and
then simulate the later trial using this posterior.  This is an
intrinsically Bayesian approach.


Installation
============


```R
# install.packages('devtools')
devtools::install_github("scientific-computing-solutions/assurance", build_vignettes = TRUE)
```


Contact <paul.metcalfe@astrazeneca.com> or
<david.ruau@astrazeneca.com> for further details.

