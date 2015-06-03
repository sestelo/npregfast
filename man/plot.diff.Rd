\name{plot.diff}
\alias{plot.diff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of the differences between the estimation of curves for two factor's levels
}
\description{
Useful for drawing the differences between the estimation of curves (initial estimate, first or second derivative) for  two factor's levels.  Missing values of factor's levels is not allowed. 
}
\method{plot}{diff}
plot.diff(model, factor2, factor1, der = NULL, est.include = FALSE, xlab = model$name[2], 
ylab = model$name[1], ylim = NULL, main = NULL, col = "black", CIcol = "grey50", ablinecol = "red", abline = TRUE, type = "l", CItype = "l", lwd = 1, CIlwd = 1.5, lty = 1, CIlty = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{allometric or nonparametric regression model obtained by \code{\link{frfast}} function.
}
  \item{factor2}{
second factor's level at which to perform the differences between curves.
}
  \item{factor1}{first factor's level at which to perform the differences between curves.
}
  \item{der}{
number or vector which determines any inference process. By default \code{der} is \code{NULL}. If this term is \code{0}, the calculate of the maximum point is for the estimate. If it is \code{1} or \code{2}, it is designed for the first or second derivative, respectively.
}
   \item{est.include}{ draw the estimates of the model. By default it is \code{FALSE}.
}
   \item{xlab}{a title for the x axis. 
}
  \item{ylab}{
a title for the y axis.
}
   \item{ylim}{the \code{y} limits of the plot.
}
   \item{main}{ an overall title for the plot.
}
  \item{col}{
a specification for the default plotting color.
}
  \item{CIcol}{
a specification for the default confidence intervals plotting color.
}
   \item{ablinecol}{the color to be used for \code{abline}.
}
   \item{abline}{ draw an horizontal line into the plot of the second derivative of the model. By default it is \code{TRUE}.
}
  \item{type}{
what type of plot should be drawn. Possible types are, \code{p} for points, \code{l} for lines, \code{o} for overplotted, etc. See details in \code{\link{par}}.
}
  \item{CItype}{
what type of plot should be drawn for confidence intervals. Possible types are, \code{p} for points, \code{l} for lines, \code{o} for overplotted. 
}
 \item{lwd}{the line width, a positive number, defaulting to 1. See details in \code{\link{par}}.
}
   \item{CIlwd}{the line width for confidence intervals, a positive number, defaulting to 1.
}
   \item{lty}{the line type. Line types can either be specified as an integer (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash). See details in \code{\link{par}}.
}
   \item{CIlty}{the line type for confidence intervals. Line types can either be specified as an integer (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash). 
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

##############################################
# Nonparametric regression with interactions
##############################################
fit2<-frfast(DW~RC:F,data=barnacle)
plot.diff(fit2,factor2=2,factor1=1)
plot.diff(fit2,factor2=2,factor1=1,der=1,col="red",CIcol="green")
plot.diff(fit2,2,1,der=c(0,1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line