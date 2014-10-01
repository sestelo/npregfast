\name{plot.diff}
\alias{plot.diff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of the differences between the estimation of curves for two factor's levels
}
\description{
Useful for drawing the differences between the estimation of curves (initial estimate, first or second derivative) for  two factor's levels. Missing values of factor's levels is not allowed.
}
\usage{
plot.diff(model, factor1, factor2, der = NULL, xlab = "x", 
ylab = "y", col = "black", ICcol = "grey", type = "l", ICtype = "l", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{allometric or nonparametric regression model obtained by \code{\link{frfast}} function.
}
  \item{factor1}{first factor's level at which to perform the differences between curves.
}
  \item{factor2}{
second factor's level at which to perform the differences between curves.
}
  \item{der}{
number or vector which determines any inference process. By default \code{der} is \code{NULL}. If this term is \code{0}, the calculate of the maximum point is for the estimate. If it is \code{1} or \code{2}, it is designed for the first or second derivative, respectively.
}
   \item{xlab}{a title for the x axis. 
}
  \item{ylab}{
a title for the y axis.
}
  \item{col}{
a specification for the default plotting color.
}
  \item{ICcol}{
a specification for the default confidence intervals plotting color.
}
  \item{main}{
an overall title for the plot.
}
  \item{type}{
what type of plot should be drawn. Possible types are, \code{p} for points, \code{l} for lines, \code{o} for overplotted, etc. See details in ?\code{par}.

}
  \item{ICtype}{
what type of plot should be drawn for confidence intervals. Possible types are, \code{p} for points, \code{l} for lines, \code{o} for overplotted. See details in ?\code{par}.
}
  \item{\dots}{
other options.
}
}
\details{
simply produce a plot.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
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

fit2<-diff(x=barnacle$RC,y=barnacle$DW,f=barnacle$F)
plot.diff(fit2)
plot.diff(fit2,2,1,der=0)
plot.diff(fit2,2,1,der=1,col="red",ICcol="green")
plot.diff(fit2,2,1,der=c(0,1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
