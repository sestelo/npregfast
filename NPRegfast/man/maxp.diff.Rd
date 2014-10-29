\name{maxp.diff}
\alias{maxp.diff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Differences between the estimation of maximum points  for two factor's levels
}
\description{
Differences between the estimation of \code{\link{maxp}} for two  factor's levels. \code{maxp}, a returned element  of class \code{\link{frfast}}, is the value of covariate \code{x} which maximizes the  estimate, first or second derivative.
}
\usage{
maxp.diff(model, factor2 = NULL, factor1 = NULL, der = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{parametric or nonparametric regression model obtained by \code{\link{frfast}} function.
}
  \item{factor1}{first factor's level at which to perform the differences between maximum points.
}
  \item{factor2}{second factor's level at which to perform the differences between maximum points.
}
  \item{der}{number which determines any inference process. By default \code{der} is \code{NULL}. If this term is \code{0}, the calculate of the differences for maximum point is for the estimate. If it is \code{1} or \code{2}, it is designed for the first or second derivative, respectively.
}
}
\details{
Differences are calculated by subtracting a factor relative to another (\eqn{factor2 - factor1}).  By default \code{factor2} and \code{factor1} are \code{NULL}, so the differences calculated are for all possible combinations between two factors.
}
\value{
An object is returned with the following element:
  \item{maxp.diff}{a table with a couple of factor's level where it is used to calculate the differences between maximum points, and their 95\% interval confidence (for the estimation, first and second derivative).}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
library(NPRegfast)
data(barnacle)

##############################################
# Nonparametric regression with interactions
##############################################
fit2<-frfast(DW~RC:F,data=barnacle)
maxp.diff(fit2)
maxp.diff(fit2,der=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line