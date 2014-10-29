\name{NPRegfast-package}
\Rdversion{1.1}
\alias{NPRegfast-package}
\alias{NPRegfast}
\docType{package}
\title{Nonparametric estimation by using local linear kernel smoothers
}
\description{
This package provides a method for obtain nonparametric estimates using local linear kernel smoothers.

Particular features of the package are facilities for fast smoothness estimation, and the calculation of their first and second derivative. Users can define the smoothers parameters. Confidences intervals calculation is provided by bootstrap methods. Binning techniques were applied to speed up computation in the estimation and testing processes.
}
\details{
\tabular{ll}{
Package: \tab NPRegfast\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-01-10\cr}

\code{NPRegfast} provides functions for nonparametric regression models \code{\link{frfast}}, \code{\link{plot.frfast}}. The term \code{frfast} is taken to include any nonparametric regression estimated by local lineal kernel smoothers. A number of other functions such \code{\link{summary.frfast}} are also provided, for extracting information from a fitted \code{frfast} Object. 

For a listing of all routines in the NPRegfast package type:
\code{library(help="NPRegfast")}. For an overview of the NPRegfast package see NPRegfast-package.
}
\author{
Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.

Maintainer: Marta Sestelo <sestelo@uvigo.es>
}
\references{
Efron, B. (1979). Bootstrap methods: another look at the jackknife. \emph{Annals of Statistics}, \bold{7}, 1--26.

Efron, E. and Tibshirani, R. J. (1993). \emph{An introduction to the Bootstrap}. Chapman and Hall, London.

Sestelo, M. and Roca-Pardinas, J. (2011). A new approach to estimation of length-weight relationship of \eqn{Pollicipes}  \eqn{pollicipes} (Gmelin, 1789) on the Atlantic coast of Galicia (Northwest Spain): some aspects of its biology and management. \emph{Journal of Shellfish Research}, \bold{30}(3), 939--948.

Wand, M. P. and Jones, M. C. (1995). \emph{Kernel Smoothing}. Chapman & Hall, London.

}

\keyword{ package}
\keyword{ models}
\keyword{ nonparametric regression}
\keyword{ local linear kernel smoothers}
\keyword{ Interactions factor by curve }



\examples{
## See examples for frfast
}
