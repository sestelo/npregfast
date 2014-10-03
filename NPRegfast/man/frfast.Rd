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
frfast(x, y, f = NULL, model = 1, h = 75, w = NULL, p = 2,
       kbin = 100, nc = NULL, ncmax = 5, ikernel = 1, iopt = 1, 
       nboot = 500, c2 = NULL, rankl = NULL, ranku = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of \code{x} data. Missing values are not allowed.
}  
  \item{y}{vector of \code{y} data. Missing values are not accepted.
}  
  \item{f}{vector of \code{factor} data to take into account in the model.
}  
  \item{model}{the nonparametric regression fitting by local linear kernel smoothers (\code{model = 1}). \code{model = 2} is used to fit an allometric model.
}  
  \item{h}{the kernel bandwidth smoothing parameter. Large values of bandwidth make smoother estimates, smaller values of bandwidth make less smooth estimates. The default is a bandwidth compute by cross validation.
} %specifically bla bla
    \item{w}{weights on the data.
}  
  \item{p}{degree of a polynomial.
} 
  \item{kbin}{number of binning nodes over which the function is to be estimated.
}  
%  \item{nc}{
%%     ~~Describe \code{nc} here~~
%}
%%  \item{ncmax}{
%%     ~~Describe \code{ncmax} here~~
%}
  \item{ikernel}{numeric which determines the smoothing kernel. By default \code{ikernel = 1}, this is, the Gaussian density function.   
}
%%  \item{iopt}{
%%     ~~Describe \code{iopt} here~~
%%}
  \item{nboot}{number of bootstrap repeats.
}
%%  \item{c2}{
%%     ~~Describe \code{c2} here~~
%%}
  \item{rankle}{number or vector specifying the minimum value for an interval at which to search the \code{x} value which maximizes the estimate, first or second derivative  (for each level). The default is the minimum data value.
}
  \item{ranku}{number or vector specifying the maximum value for an interval at which to search the \code{x} value which maximizes the estimate, first or second derivative  (for each level). The default is the maximum data value.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
An object is returned with the following elements:
  \item{x}{vector of values of the grid points at which model is to be estimate.}
  \item{p}{matrix of values of the grid points at which to compute the estimate, their first and second derivative.}
  \item{pl}{lower values of  95\% confidence interval for the estimate, their first and second derivative.}
  \item{pu}{upper values of  95\% confidence interval for the estimate, their first and second derivative.}
  \item{diff}{differences between the estimation values of a couple of levels (i. e. level 2 - level 1). The same procedure for their first and second derivative.}
  \item{diffl}{lower values of 95\% confidence interval for the differences between the estimation values of a couple of levels. It is performed for their first and second derivative.}
  \item{diffu}{upper values of 95\% confidence interval for the differences between the estimation values of a couple of levels. It is performed for their first and second derivative.}
  \item{boot}{number of bootstrap repeats.}
  \item{n}{total number of data}
  \item{dp}{degree of a polynomial.}
  \item{h}{the kernel bandwidth smoothing parameter.}
  \item{grid}{the number of equally espaced points at which to estimate the curves.}
  \item{mod}{factor's level for each data.}
  \item{xdata}{original x values}
  \item{data}{original y values}
  \item{w}{weights on the data.}
  \item{fact}{factor's level.}
  \item{nf}{number of levels.}
  %%\item{c2}{}
  %%\item{ncmax}{}
  %%\item{nc}{}
  \item{kbin}{number of binning nodes over which the function is to be estimated.}%%number of equally espaced points at which to estimate the curves. 
  \item{ikernel}{character which determines the smoothing kernel. By default \code{ikernel = 1}, this is, the Gaussian density function.}
  \item{max}{value of covariate \code{x} which maximizes the  estimate, first or second derivative.}
  \item{maxl}{lower value of 95\% confidence interval for the value \code{max}.}
  \item{maxi}{upper value of 95\% confidence interval for the value \code{max}.}
  \item{maxboot}{values of the covariate \code{x} which maximizes the  estimate, first or second derivative for each bootstrap repeats.}
  \item{diffmax}{differences between the estimation of \code{max} for a couple of levels (i. e. level 2 - level 1). The same procedure for their first and second derivative.}
  \item{diffmaxl}{lower value of 95\% confidence interval for the value \code{diffmax}.}
  \item{diffmaxu}{upper value of 95\% confidence interval for the value \code{diffmax}.}
  \item{repboot}{matrix of values of the grid points at which to compute the estimate, their first and second derivative for each bootstrap repeat.}
  \item{ranku}{minimum value for an interval at which to search the \code{x} value which maximizes the estimate, first or second derivative  (for each level). The default is the minimum data value.}
  \item{rankl}{maximum value for an interval at which to search the \code{x} value which maximizes the estimate, first or second derivative  (for each level). The default is the maximum data value.}
  %%\item{call}{}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardi\~{n}as.
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

fit<-frfast(x=barnacle$RC,y=barnacle$DW)
fit

summary(fit)

# Change the number of binning nodes and bootstrap replicates
fit<-frfast(x=barnacle$RC,y=barnacle$DW,kbin=200,nboot=1000)



##############################################
# Nonparametric regression with interactions
##############################################

fit2<-frfast(x=barnacle$RC,y=barnacle$DW,f=barnacle$F)
fit2

summary(fit2)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
