# npregfast: Nonparametric Estimation of Regression Models with Factor-by-Curve Interactions

[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/grand-total/npregfast)](https://cran.r-project.org/package=npregfast)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/npregfast)](https://cran.r-project.org/package=npregfast)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13889485.svg)](https://doi.org/10.5281/zenodo.13889485)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/last-month/npregfast?color=ff69b4)](https://cran.r-project.org/package=npregfast)

```npregfast``` is an R package for obtain nonparametric estimates of regression models 
with or without factor-by-curve interactions using local polynomial kernel smoothers or splines. 
Additionally, a parametric model (allometric model) can be estimated.
Particular features of the package are facilities for fast smoothness
estimation, and the calculation of their first and second derivative. Users can 
define the smoothers parameters. Confidence intervals calculation is provided 
by bootstrap methods. Binning techniques were applied to speed up computation 
in the estimation and testing processes.




You can view a live interactive demo to see part of its 
capabilities at http://sestelo.shinyapps.io/npregfast.

## Installation
```npregfast``` is available through both CRAN and GitHub.

Get the released version from CRAN:
```
install.packages("npregfast")
```

Or the development version from GitHub:
```
# install.packages("devtools")
devtools::install_github("sestelo/npregfast")
```
