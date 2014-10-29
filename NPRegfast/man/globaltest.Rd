\name{globaltest}
\alias{globaltest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Testing the equality of M curves specific to each level
}
\description{
\code{globaltest}  can be used to test the equality of the \eqn{M} curves specific to each level. This bootstrap based test assumes the  following null hypothesis

\eqn{H_0: m_1 = \ldots = m_M}

Note that, if \eqn{H_0} is not rejected, then the equality of critical points will also accepted. 

To test the null hypothesis, it is used an statistic, \eqn{T}, based on direct nonparametric estimates of the curves. If the null hypothesis is true, the \eqn{T} value should be close to zero but is generally greater. The test rule based on \eqn{T} consists of rejecting the null hypothesis if \eqn{T > T^{1- \alpha}}, where \eqn{T^p} is the empirical \eqn{p}-percentile of \eqn{T} under the null hypothesis. }


\usage{
globaltest(formula,data = data, der = NULL, weights = NULL, nboot = 200, h = -1.0, nh = 30, kernel = 'epanech', p = 3, kbin = 100)
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
  \item{kernel}{character which determines the smoothing kernel. By default \code{kernel='epanech'}, this is, the Epanechnikov density function. Also, several types of kernel funcitons can be used:  triangular and Gaussian density function, with \code{'triang'} and \code{'gaussian'} term, respectively.
}   
  \item{p}{degree of a polynomial.
} 
  \item{kbin}{number of binning nodes over which the function is to be estimated.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
The \eqn{T} valueand the \eqn{p}-value  are returned. Additionally, it is shown the decision, accepted or rejected, of the global test. The null hypothesis is rejected if the \eqn{p}-value\eqn{< 0.05}.  
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
glocaltest(DW~RC:F,data=barnacle, der=0)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
