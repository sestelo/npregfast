\name{localtest}
\alias{localtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Testing the equality of critical points
}
\description{
\code{localtest} can be used to test the equality of the \eqn{M} critical points estimated from the respective level-specific curves. Note that, even if the curves and/or their derivatives are different, it is possible for these points to be equal. 
For instance, taking the maxima of the first derivatives into account, interest lies in testing the following null hypothesis

\eqn{H_0: x_{01} = \ldots = x_{0M}}
The above hypothesis is true if \eqn{D = x_{0j} - x_{0k} = 0} where 

\eqn{(j,k) = arg max_{1\le l < m \le M} |x_{0j} - x_{0k}|}

otherwise \eqn{H_0} is false. It is important to highlight that, in practice, the true \eqn{x_{0j}} are not known, and consequently neither is \eqn{D}, so an estimate \eqn{\hat D = \hat x_{0j} − \hat x_{0k}} is used, where, in general, \eqn{\hat x_{0l}} are the estimates of \eqn{x_{0l}} based on the estimated curves \eqn{\hat m_l} with \eqn{l = 1, \ldots , M}.
Needless to say, since \eqn{\hat D} is only an estimate of the true D, the sampling uncertainty of these estimates needs to be acknowledged. Hence, a confidence interval \eqn{(a,b)} is created for \eqn{D} for a specific level of confidence (e.g., 95\%). 
}
\usage{
localtest(formula, data = data, der = NULL, weights = NULL, nboot = 200, h = -1.0, nh = 30, kernel = 1, p = 3, kbin = 100, ranku = NULL, rankl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class \code{formula}: a sympbolic description of the model to be fitted.
}  
  \item{data}{a data frame or matrix containing the model response variable and covariates required by the \code{formula}.
}  
  \item{der}{number which determines any inference process. By default \code{der} is \code{NULL}. If this term is \code{0}, the calculate of the differences for maximum point is for the estimate. If it is \code{1} or \code{2}, it is designed for the first or second derivative, respectively.
}
  \item{weights}{prior weights on the data.
} 
\item{nboot}{number of bootstrap repeats.
} 
  \item{h}{the kernel bandwidth smoothing parameter. Large values of bandwidth make smoother estimates, smaller values of bandwidth make less smooth estimates. The default is a bandwidth compute by cross validation.
} 
\item{nh}{integer number of equally-spaced bandwidth on which the \code{h} is discretised, to speed up computation.
}
\item{kernel}{character which determines the smoothing kernel. By default \code{kernel = 'epanech' }, this is, the Epanechnikov density function. Also, several types of kernel funcitons can be used:  triangular and Gaussian density function, with 'triang' and 'gaussian' term, respectively.
}  
  \item{p}{degree of a polynomial.
} 
  \item{kbin}{number of binning nodes over which the function is to be estimated.
}
\item{rankl}{number or vector specifying the minimum value for an interval at which to search the \code{x} value which maximizes the estimate, first or second derivative  (for each level). The default is the minimum data value.
}
  \item{ranku}{number or vector specifying the maximum value for an interval at which to search the \code{x} value which maximizes the estimate, first or second derivative  (for each level). The default is the maximum data value.
}
}  
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
The estimate of \eqn{D} value is returned and its confidence interval for a specific-level of confidence, i.e. 95\%. Additionally, it is shown the decision, accepted or rejected,  of the local test. Based on the null hypothesis is rejected if a zero value is not within the interval. 
}
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

##############################################
# Nonparametric regression with interactions
##############################################

localtest(DW~RC:F,data=barnacle, der=0)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
