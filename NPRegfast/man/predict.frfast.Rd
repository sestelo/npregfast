\name{predict.frfast}
\alias{predict.frfast}
\alias{predict.frfast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction from fitted \code{frfast} model
}
\description{
Takes a fitted \code{frfast} object  and produces predictions (and optionally estimates standard errors of those predictions) from a fitted model with interactions or without interactions.} 
\method{predict}{frfast}{
predict.frfast(model,newdata,fac=NULL,der=NULL,...)
}
\arguments{
\item{model}{a fitted \code{frfast} object as produced by \code{frfast()}.
}
\item{newdata}{a data frame containing the values of the model covariates at which predictions are required.  If newdata is provided, then it should contain all the variables needed for prediction: a warning is generated if not.
}
\item{fac}{number or vector which determines the level to take into account in the plot. By default is \code{NULL}.
}
\item{der}{number or vector which determines any inference process. By default \code{der} is \code{NULL}. If this term is \code{0}, the plot show the initial estimate. If it is \code{1} or \code{2}, it is designed for the first or second derivative, respectively.
}
}

\details{
	\code{print.frfast} tries to be smart about \code{summary.frfast}.
}
\value{
	\code{predict.frfast} computes and returns a list containing predictions of the estimates, first and second derivative, with their 95 \% confidence intervals. 
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
}
\examples{
library(NPRegfast)
data(barnacle)

################################################
# Nonparametric regression without interactions
################################################
fit<-frfast(DW~RC,data=barnacle)

nd=data.frame(RC=c(10,14,18))
predict.frfast(fit,newdata=nd)


##############################################
# Nonparametric regression with interactions
##############################################
fit2<-frfast(DW~RC:F,data=barnacle)

nd2=data.frame(RC=c(10,15,20))
pred<-predict.frfast(fit2,newdata=nd2)
pred$Level_1$Estimation
predict.frfast(fit2,newdata=nd2,der=0,fac=2)
}