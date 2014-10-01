\name{summary.frfast}
\alias{summary.frfast}
\alias{print.frfast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing fits of frfast class
}
\description{
Takes a fitted \code{frfast} object produced by \code{frfast()} and produces various useful
summaries from it.} 
\usage{
summary.frfast(model)
}
\arguments{\item{model}{a fitted \code{frfast} object as produced by \code{frfast()}.}}



\details{\code{print.frfast} tries to be smart about \code{summary.frfast}.}



\value{\code{summary.frfast} computes and returns a list of summary information for a fitted \code{frfast} object.

\item{model}{type of estimate.}
\item{h}{the kernel bandwidth smoothing parameter.}
\item{dp}{degree of a polynomial.}
\item{nboot}{number of bootstrap repeats.}
\item{grid}{number of binning nodes over which the function is to be estimated.}
\item{n}{total number of data.}
\item{fmod}{factor's levels.}

}



%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.}

\examples{
library(NPRegfast)
data(barnacle)
fit<-frfast(x=barnacle$RC,y=barnacle$DW)
summary.frfast(fit)

fit2<-frfast(x=barnacle$RC,y=barnacle$DW,f=barnacle$F)
summary.frfast(fit2)}


