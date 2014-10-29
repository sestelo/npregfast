\name{frfast}
\alias{frfast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting nonparametric models
}
\description{
\code{frfast} is used to fit nonparametric models by using local linear kernel smoothers.
}
\usage{
frfast(formula, data = data, model = 'np', h = -1.0, nh = 30, weights = NULL, kernel = 'epanech', p = 3, kbin = 100, nboot = 500, rankl = NULL, ranku = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class \code{formula}: a sympbolic description of the model to be fitted. The details of model specification are given under 'Details'.%\item{f}{vector of \code{factor} data to take into account in the model.} 
}  
  \item{data}{ a data frame or matrix containing the model response variable and covariates required by the \code{formula}. 
}   
  \item{model}{type model used: \code{model = 'np'}  nonparametric regression model with local linear kernel smoothers, \code{model = 'allo'} the  allometric model.
}  
  \item{h}{the kernel bandwidth smoothing parameter. Large values of bandwidth make smoother estimates, smaller values of bandwidth make less smooth estimates. The default is a bandwidth compute by cross validation.
} %specifically bla bla
\item{nh}{integer number of equally-spaced bandwidth on which the \code{h} is discretised, to speed up computation. %the grid bandwidth´s that it will be used to calculate the grid of bandwidths. %In this grid, it will be selected the optimum bandwidth by cross-validation.If the optimum bandwidth value is close to 0, we will obtain rough estimates; when it is close to 1, we will obtain smooth estimates.
}
\item{weights}{prior weights on the data.
} 
\item{kernel}{character which determines the smoothing kernel. By default \code{kernel = 'epanech' }, this is, the Epanechnikov density function. Also, several types of kernel funcitons can be used:  triangular and Gaussian density function, with 'triang' and 'gaussian' term, respectively.
} 
  \item{p}{degree of polynomial used.  Its value must be greater than or equal to the value of drv. The default value is of degree is drv + 1.
} 
  \item{kbin}{number of binning nodes over which the function is to be estimated. %number of equally-spaced grid points over which the function is to be estimated.
}  
  \item{nboot}{number of bootstrap repeats.
}
  \item{rankl}{number or vector specifying the minimum value for an interval at which to search the \code{x} value which maximizes the estimate, first or second derivative  (for each level). The default is the minimum data value.
}
  \item{ranku}{number or vector specifying the maximum value for an interval at which to search the \code{x} value which maximizes the estimate, first or second derivative  (for each level). The default is the maximum data value.
}}
%  \item{nc}{
%%     ~~Describe \code{nc} here~~
%}
%%  \item{ncmax}{
%%     ~~Describe \code{ncmax} here~~
%}
%  \item{ikernel}{numeric which determines the smoothing kernel. By default \code{ikernel = 1}, this is, the Gaussian density function.   
%}
%%  \item{iopt}{
%%     ~~Describe \code{iopt} here~~
%%}
%%  \item{c2}{
%%     ~~Describe \code{c2} here~~
%%}
\details{
The models fit by \code{frfast} function are specified in a compact symbolic form. The \~operator is basic in the formation of such models. An expression of the form y ~ model is interpreted as a specification that the response y is modelled by a linear predictor specified symbolically by model. The terms themselves consist of variable and factor names separated by : operators. Such a term is interpreted as the interaction of all the variables and factors appearing in the term.
}
\value{
An object is returned with the following elements:
\item{x}{vector of values of the grid points at which model is to be estimate.
}
\item{p}{matrix of values of the grid points at which to compute the estimate, their first and second derivative.
}
\item{pl}{lower values of  95\% confidence interval for the estimate, their first and second derivative.
}
\item{pu}{upper values of  95\% confidence interval for the estimate, their first and second derivative.
}
\item{diff}{differences between the estimation values of a couple of levels (i. e. level 2 - level 1). The same procedure for their first and second derivative.
}
\item{diffl}{lower values of 95\% confidence interval for the differences between the estimation values of a couple of levels. It is performed for their first and second derivative.
}
\item{diffu}{upper values of 95\% confidence interval for the differences between the estimation values of a couple of levels. It is performed for their first and second derivative.
}
\item{nboot}{number of bootstrap repeats.
}
\item{n}{total number of data
}
\item{dp}{degree of polynomial used.
}
\item{h}{the kernel bandwidth smoothing parameter.
}
\item{fmod}{factor's level for each data.
}
\item{xdata}{original x values
}
\item{ydata}{original y values
}
\item{w}{weights on the data.
}
\item{nf}{number of levels.
}
\item{kbin}{number of binning nodes over which the function is to be estimated.%numberofequallyespacedpointsatwhich to estimate the curves.
}
\item{pvalue}{it is \code{NULL} when the nonparamentric model is fitted. However, if the p-value of the allometric test. 
} 
\item{max}{value of covariate \code{x} which maximizes the  estimate, first or second derivative.
}
\item{maxu}{upper value of 95\% confidence interval for the value \code{max}.
}
\item{maxl}{lower value of 95\% confidence interval for the value \code{max}.
}
\item{diffmax}{differences between the estimation of \code{max} for a couple of levels (i. e. level 2 - level 1). The same procedure for their first and second derivative.
}
\item{diffmaxu}{upper value of 95\% confidence interval for the value \code{diffmax}.
}
\item{diffmaxl}{lower value of 95\% confidence interval for the value \code{diffmax}.
}
\item{statistic}{the value of the test statistic.
}
\item{repboot}{matrix of values of the grid points at which to compute the estimate, their first and second derivative for each bootstrap repeat.
}
\item{ranku}{minimum value for an interval at which to search the \code{x} value which maximizes the estimate, first or second derivative  (for each level). The default is the minimum data value.
}
\item{rankl}{maximum value for an interval at which to search the \code{x} value which maximizes the estimate, first or second derivative  (for each level). The default is the maximum data value.
}
\item{nmodel}{type model used: \code{model = 1} the nonparametric model, \code{model = 2} the allometric model. 
} 
\item{label}{labels of the variables in the model.
}
\item{numlabel}{number of labels.
}
\item{kernel}{character which determines the smoothing kernel.
}
\item{name}{name of the variables in the model.
}
\item{formula}{a sympbolic description of the model to be fitted.
}}
  %\item{nh}{}
  %%\item{call}{}
  % \item{maxboot}{values of the covariate \code{x} which maximizes the  estimate, first or second derivative for each bootstrap repeats.}
  %\item{fact}{factor's level.}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
}
%%\note{
%%  ~~further notes~~
%%}
%% ~Make other sections like Warning with \section{Warning }{....} ~
%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
library(NPRegfast)
data(barnacle)


################################################
# Nonparametric regression without interactions
################################################
fit<-frfast(DW~RC,data=barnacle)
fit

summary(fit)

# Change the number of binning nodes and bootstrap replicates
fit<-frfast(DW~RC,data=barnacle,kbin=200,nboot=1000)

##############################################
# Nonparametric regression with interactions
##############################################
fit2<-frfast(DW~RC:F,data=barnacle)
fit2

summary(fit2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
