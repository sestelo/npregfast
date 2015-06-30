#' \code{npregfast}: Nonparametric Estimation of 
#' Regression Models with Factor-By-Curve Interactions. 
#'
#'
#' This package provides a method for obtain nonparametric estimates of regression
#' models using local linear kernel smoothers. Particular features of the package
#' are facilities for fast smoothness estimation, and the calculation of their 
#' first and second derivative. Users can define the smoothers parameters. 
#' Confidences intervals calculation is provided by bootstrap methods. 
#' Binning techniques were applied to speed up computation in the estimation 
#' and testing processes.
#'
#' @name npregfast
#' @docType package
#' @details \tabular{ll}{ Package: \tab npregfast\cr Type: \tab Package\cr
#' Version: \tab 1.0\cr Date: \tab 2015-06-29\cr License: \tab MIT + file LICENSE\cr}
#'
#' CAMBIAR \code{npregfast} provides functions for nonparametric regression models 
#' \code{\link{frfast}}, \code{\link{plot.frfast}}. The term \code{frfast} is 
#' taken to include any nonparametric regression estimated by local lineal 
#' kernel smoothers. A number of other functions such \code{\link{summary.frfast}} 
#' are also provided, for extracting information from a fitted \code{frfast} object.
#'
#' For a listing of all routines in the NPRegfast package type:
#' \code{library(help="npregfast")}. 
#' 
#' 
#' @author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#' 
#' @references 
#' Efron, B. (1979). Bootstrap methods: another look at the jackknife. 
#' Annals of Statistics, 7, 1--26.
#' 
#' Efron, E. and Tibshirani, R. J. (1993). An introduction to the Bootstrap. 
#' Chapman and Hall, London.
#' 
#' Huxley, J. S. (1924). Constant differential growth-ratios and their 
#' significance. Nature, 114:895--896.
#' 
#' Sestelo, M. (2013). Development and computational implementation of 
#' estimation and inference methods in flexible regression models. 
#' Applications in Biology, Engineering and Environment. PhD Thesis, Department
#' of Statistics and O.R. University of Vigo.
#' 
#' Sestelo, M. and Roca-Pardinas, J. (2011). A new approach to estimation of 
#' length-weight relationship of \eqn{Pollicipes}  \eqn{pollicipes} 
#' (Gmelin, 1789) on the Atlantic coast of Galicia (Northwest Spain): some 
#' aspects of its biology and management. Journal of Shellfish Research, 
#' 30(3), 939--948.
#' 
#' Wand, M. P. and Jones, M. C. (1995). Kernel Smoothing. Chapman & Hall, London.
#'
#'
#'


NULL