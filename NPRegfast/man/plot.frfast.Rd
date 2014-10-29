\name{plot.frfast}
\alias{plot.frfast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Visualization of frfast objects}
%Useful for drawing estimation, first and second derivative by factors.

\description{
Useful for drawing the estimation, first and second derivative  (for each factor)
}
\usage{
\method{plot}{frfast}(model, fac = NULL, der = NULL, points = TRUE, xlab = model$name[2], ylab = model$name[1], ylim = NULL, main = NULL, col = "black", CIcol = "black", ablinecol = "red", abline = TRUE,type = "l", CItype = "l", lwd = 2, CIlwd = 1,lty = 1, CIlty = 2, ...)
}
\arguments{
  \item{model}{\code{frfast} object.
}
  \item{fac}{number or vector which determines the level to take into account in the plot. By default is \code{NULL}.
}
  \item{der}{number or vector which determines any inference process. By default \code{der} is \code{NULL}. If this term is \code{0}, the plot show the initial estimate. If it is \code{1} or \code{2}, it is designed for the first or second derivative, respectively.
}
  \item{points}{draw the original data into the plot. By default it is \code{TRUE}. 
}
  \item{xlab}{
  	a title for the x axis. 
}
  \item{ylab}{
a title for the y axis.
}
  \item{ylim}{the \code{y} limits of the plot.
}
  \item{main}{an overall title for the plot.
}
  \item{col}{
a specification for the default plotting color.
}
  \item{CIcol}{
a specification for the default confidence intervals plotting color.
}
  \item{ablinecol}{the color to be used for \code{abline}.
}
  \item{abline}{draw an horizontal line into the plot of the second derivative of the model.
}
  \item{type}{
what type of plot should be drawn. Possible types are, \code{p} for points, \code{l} for lines, \code{o} for overplotted, etc. See details in \code{\link{par}}.
}
  \item{CItype}{
what type of plot should be drawn for confidence intervals. Possible types are, \code{p} for points, \code{l} for lines, \code{o} for overplotted.
}
  \item{lwd}{the line width, a positive number, defaulting to 1.  See details in \code{\link{par}}.
}
  \item{CIlwd}{the line width for confidence intervals, a positive number, defaulting to 1. 
}
\item{lty}{the line type. Line types can either be specified as an integer (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash).  See details in \code{\link{par}}.
}
\item{CIlty}{the line type for confidence intervals. Line types can either be specified as an integer (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash). 
}
  \item{\dots}{
other options.
}
}
%\details{
%% ~~ If necessary, more details than the description above ~~
%}
\value{
simply produce a plot.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
}
%%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

library(NPRegfast)
data(barnacle)

################################################
# Nonparametric regression without interactions
################################################
fit<-frfast(DW~RC,data=barnacle)
plot(fit,der=0)
plot(fit,der=1,col="red",CIcol="blue",points=FALSE)

##############################################
# Nonparametric regression with interactions
##############################################
fit2<-frfast(DW~RC:F,data=barnacle)
plot(fit2)
plot(fit2,der=0,fac=2)
plot(fit2,der=1,col="red",CIcol="green")
plot(fit2,der=c(0,1),fac=c(1,2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line