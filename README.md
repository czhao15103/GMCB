This repository includes an R package for implementing joint mean-covariance estimation using 
the generalized mean-covariance bridge (GMCB) prior and code for replicating the simulations
done in _A Bayesian Generalized Bridge Regression Approach to Covariance Estimation in the Presence of Covariates_.

#### GMCB

This subdirectory contains the files necessary to install the GMCB package, which
can be installed with `devtools::install_github`:

```r
library(devtools)
devtools::install_github(repo = "czhao15103/GMCB/GMCB")
```r

#### Simulations

This subdirectory contains files to replicate the simulations in the paper. The file `HSGHS.m`, which
was used to obtain results from the HS-GHS method, was obtained from the [HS_GHS repostitory](https://github.com/liyf1988/HS_GHS).