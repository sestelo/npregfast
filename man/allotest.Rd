\name{allotest} %Revisar mejor
\alias{allotest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap based test for testing an allometric model
}
\description{
In order to facilitate the choice of a model appropriate to the data while at the same time endeavouring to minimise the loss of information,  a bootstrap-based procedure, that test whether the data can be modelled by an allometric model, was developed.  Therefore, \code{allotest} tests the null hypothesis of an allometric model taking into account the logarithm of the original variable (\eqn{X^* = log(X)}{} and \eqn{Y^* =log (Y)}{}). Based on a general model of the type
\deqn{Y^*=m(X^*)+\varepsilon}{} the aim here is to test the null hypothesis of an allometric model 

\deqn{H_0 = m(x^*) =  a^*+ b^* x^*}{} 

\eqn{vs.}{} general hypothesis \eqn{H_1}{}, with \eqn{m}{} being an unknown nonparametric function; or analogously,

\deqn{H_1: m(x^*)= a^*+ b^* x^* + g(x^*)}{}

with \eqn{g(x^*)}{} being an unknown function not equal to zero. To implement this test we have used the wild bootstrap. 
}
\usage{
allotest(formula, data = data, nboot=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class \code{formula}: a sympbolic description of the model to be fitted. 
}
  \item{data}{a data frame or matrix containing the model response variable and covariates required by the \code{formula}.
}
 \item{nboot}{number of bootstrap repeats.
}
}
%\details{
%
%}
\value{
An object is returned with the following elements:
  \item{value}{the p-value of the test.}
  \item{statistic}{the value of the test statistic.}
}
%%\references{
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
allotest(DW~RC:F,data=barnacle)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
