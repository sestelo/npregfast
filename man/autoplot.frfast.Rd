% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autoplot.frfast.R
\name{autoplot.frfast}
\alias{autoplot.frfast}
\title{Visualization of \code{frfast} objects with ggplot2 graphics}
\usage{
\method{autoplot}{frfast}(
  object = model,
  fac = NULL,
  der = 0,
  diffwith = NULL,
  points = TRUE,
  xlab = model$name[2],
  ylab = model$name[1],
  ylim = NULL,
  main = NULL,
  col = "black",
  CIcol = "black",
  CIlinecol = "transparent",
  pcol = "grey80",
  abline = TRUE,
  ablinecol = "red",
  lty = 1,
  CIlty = 2,
  lwd = 1,
  CIlwd = 1,
  cex = 1.4,
  alpha = 0.2,
  ...
)
}
\arguments{
\item{object}{\code{frfast} object.}

\item{fac}{Factor's level to be taken into account
in the plot. By default is \code{NULL}.}

\item{der}{Number which determines any inference process. 
By default \code{der} is \code{0}. If this term is \code{0}, the plot 
shows the initial estimate. If it is \code{1} or \code{2},
it is designed for the first or second derivative, respectively.}

\item{diffwith}{Factor's level used for drawing the differences respect to the 
level specified in the \code{fac} argument.  By default, \code{NULL}. 
The differences are computed for the r-th derivative specified 
in the \code{der} argument.}

\item{points}{Draw the original data into the plot. By default it is
\code{TRUE}.}

\item{xlab}{A title for the \code{x} axis.}

\item{ylab}{A title for the \code{y} axis.}

\item{ylim}{The \code{y} limits of the plot.}

\item{main}{An overall title for the plot.}

\item{col}{A specification for the default plotting color.}

\item{CIcol}{A specification for the default confidence intervals
plotting color (for the fill).}

\item{CIlinecol}{A specification for the default confidence intervals
plotting color (for the edge).}

\item{pcol}{A specification for the points color.}

\item{abline}{Draw an horizontal line into the plot of the second derivative 
of the model.}

\item{ablinecol}{The color to be used for \code{abline}.}

\item{lty}{The line type. Line types can either be specified as an integer
(0 = blank, 1 = solid (default), 2 = dashed, 3 = dotted, 4 = dotdash, 
5 = longdash, 6 = twodash).  See details in \code{\link{par}}.}

\item{CIlty}{The line type for confidence intervals. Line types can either 
be specified as an integer (0 = blank, 1 = solid (default), 2 = dashed,
3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash).}

\item{lwd}{The line width, a positive number, defaulting to 1.
See details in \code{\link{par}}.}

\item{CIlwd}{The line width for confidence intervals, a positive number, 
defaulting to 1.}

\item{cex}{A numerical value giving the amount by which plotting symbols
should be magnified relative to the default. See details in \code{\link{par}}.}

\item{alpha}{Alpha transparency for overlapping elements expressed 
as a fraction between 0 (complete transparency) and 1 (complete opacity).}

\item{\ldots}{Other options.}
}
\value{
A ggplot object, so you can use common features from 
ggplot2 package to manipulate the plot.
}
\description{
Useful for drawing the estimated regression function, 
first and second derivative (for each factor's level) using ggplot2 graphics.
 Additionally, with the 
\code{diffwith} argument it is possible to draw the differences between
 two factor's levels.
}
\examples{

library(npregfast)
library(ggplot2)


data(barnacle)

# Nonparametric regression without interactions
fit <- frfast(DW ~ RC, data = barnacle, nboot = 50) 
autoplot(fit)
autoplot(fit, points = FALSE) + ggtitle("Title")
autoplot(fit, der = 1) + xlim(4, 20)
#autoplot(fit, der = 1, col = "red", CIcol = "blue")

# Nonparametric regression with interactions
fit2 <- frfast(DW ~ RC : F, data = barnacle, nboot = 50) 
autoplot(fit2, fac = "barca")
# autoplot(fit2, der = 1, fac = "lens")

# Visualization of the differences between two factor's levels
autoplot(fit2, fac = "barca", diffwith = "lens")
# autoplot(fit2, der = 1, fac = "barca", diffwith = "lens")


#Plotting in the same graphics device
\dontrun{

if (requireNamespace("gridExtra", quietly = TRUE)) {

# For plotting two derivatives in the same graphic windows
ders <- lapply(0:1, function(x) autoplot(fit, der = x))
gridExtra::grid.arrange(grobs = ders, ncol = 2, nrow = 1)

# For plotting two levels in the same graphic windows
facs <- lapply(c("barca", "lens"), function(x) autoplot(fit2, der = 0, fac = x))
gridExtra::grid.arrange(grobs = facs, ncol = 2, nrow = 1)

}

}

}
\author{
Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
}
