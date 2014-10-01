\name{plot.frfast}
\alias{plot.frfast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Visualization of frfast objects}
%Useful for drawing estimation, first and second derivative by factors.

\description{
Useful for drawing the estimation, first and second derivative  (for each factor)
}
\usage{
plot.frfast(model, fac = NULL, der = NULL, points=TRUE, xlab = "x",
 ylab = "y", col = "black", ICcol = "grey", main = "title", type = "l", ICtype = "l", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{\code{frfast} object.
}
  \item{fac}{number or vector which determines the level to take into account in the plot. By default is \code{NULL}.
}
  \item{der}{number or vector which determines any inference process. By default \code{der} is \code{NULL}. If this term is \code{0}, the plot show the initial estimate. If it is \code{1} or \code{2}, it is designed for the first or second derivative, respectively.
}
  \item{points}{draw the original data into the plot. By default it is \code{TRUE}. 
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
fit<-frfast(x=barnacle$RC,y=barnacle$DW)
plot.frfast(fit,der=0)
plot.frfast(fit,der=1,col="red",ICcol="blue",points=FALSE)

fit2<-frfast(x=barnacle$RC,y=barnacle$DW,f=barnacle$F)
plot.frfast(fit2)
plot.frfast(fit2,der=0,fac=2)
plot.frfast(fit2,der=1,col="red",ICcol="green")
plot.frfast(fit2,der=c(0,1),fac=c(1,2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
