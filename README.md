# QPAD: Calibrating indices of avian density from non-standardized survey data

The analysis of large heterogeneous data sets of avian point-count surveys 
compiled across studies is hindered
by a lack of analytical approaches that can deal with detectability 
and variation in survey protocols.

We reformulated removal models of avian singing rates and distance 
sampling models of the effective detection
radius (EDR) to control for the effects of survey protocol and 
temporal and environmental covariates on
detection probabilities. 

* These **estimating procedures** as described 
  in Solymos et al. (2013) are implemented in the
  [`detect`](http://cran.r-project.org/package=detect) R package.

* The **estimates** of singing rates and effective
  detection distances for North American boreal forest songbird species
  is provided as part of the `QPAD` package.

Using offsets derived from these estimates can significantly reduce the 
computational burden when fitting complex models to large data sets and 
can be used with a wide range of statistical techniques for 
inference and prediction of avian densities.

## Install

Get the recent version as:

```{r}
library(devtools)
install_github("borealbirds/QPAD")
```

If you have trouble with **devtools**, try [this](http://peter.solymos.org/drat/) repo and **drat**.

## Use

Functions to calculate QPAD offsets are available in the [`wildrtrax`](https://github.com/ABbiodiversity/wildRtrax) R package or the [`qpad offsets`](https://github.com/borealbirds/qpad-offsets) repository.

## Report a problem

Use the [issue tracker](https://github.com/borealbirds/QPAD/issues)
to report a problem.

## Cite

Solymos et al. 2013. 
Calibrating indices of avian density from non-standardized survey data: 
making the most of a messy situation
_Methods in Ecology and Evolution_, **4**, 1047-1058.
