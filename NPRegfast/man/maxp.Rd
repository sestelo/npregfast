\name{maxp}
\alias{maxp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum points for the estimate, first and second derivative, with their 95\% confidence intervals
}
\description{
Value of covariate \code{x} which maximizes the  estimate, first and second derivative, for each level of the factor. 
}
\usage{
maxp(model, der = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{parametric or nonparametric regression out obtained by \code{\link{frfast}} function.
}
  \item{der}{number which determines any inference process. By default \code{der} is \code{NULL}. If this term is \code{0}, the calculate of the maximum point is for the estimate. If it is \code{1} or \code{2}, it is designed for the first or second derivative, respectively.}
}
%\details{
%%  ~~ If necessary, more details than the description above %~~
%}
\value{
An object is returned with the following elements:
  \item{Estimation}{outputs for the estimation where it is included maximum points, and their 95\% confidence intervals (for each level).}
  \item{First_der}{outputs for first derivative with maximum points and their 95\% confidence intervals (for each level).}
  \item{Second_der}{outputs for second derivative. It means, maximum points and 95\% confidence intervals (for each level).}
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

fit<-frfast(x=barnacle$RC,y=barnacle$DW)
maxp(fit)
maxp(fit,der=0)
maxp(fit,der=1)
maxp(fit,der=2)
maxp(fit,der=c(0,1))

fit2<-frfast(x=barnacle$RC,y=barnacle$DW,f=barnacle$F)
maxp(fit2)
maxp(fit2,der=0)
maxp(fit2,der=1)
maxp(fit2,der=2)
maxp(fit2,der=c(0,1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
