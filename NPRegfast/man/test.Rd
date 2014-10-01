\name{test}
\alias{test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap based test for testing a parametric allometric model
}
\description{
\code{test} is used to test the null hypothesis of an allometric model 
\eqn{H_0 = m(X) =  aX^b} \eqn{vs.} general hypothesis \eqn{H_1} being \eqn{m} an unknown nonparametric function. To implement this test we have used the wild bootstrap.


%%Performs a test of hypotheses where  null hypothesis is that  allometric model fits  data. Alternative hypothesis is a nonparametric model 
}
\usage{
test(x, y, f = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of \code{x} data. Missing values are not allowed.
}
  \item{y}{vector of \code{y} data. Missing values are not accepted.
}
  \item{f}{vector of \code{factor} data  to take into account in the model.
}
}
%%\details{
%%}
\value{
An object is returned with the following elements:
  \item{value}{the p-value of the test (for each factor).}
  \item{statistic}{the value of the statistic obtained by using likelihood ratio test.}
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
test(x=barnacle$RC,y=contrast$DW)
test(x=barnacle$RC,y=barnacle$DW,f=barnacle$F)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
