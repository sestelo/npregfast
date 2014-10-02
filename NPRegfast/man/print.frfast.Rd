\name{print.frfast}
\alias{print.frfast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{hola git, como vas? clase
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.frfast(model, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (model, ...) 
{
    cat("\nCall:\n")
    print(model$call)
    cat("", "\n")
    if (model$modelo == 1) {
        m = "Nonparametric"
    }
    else {
        m = "Allometric"
    }
    cat("*********************************************", "\n")
    cat(m, "Model", "\n")
    cat("*********************************************", "\n")
    cat("\nNumber of Observations: ")
    cat(format(model$n))
    cat("\n")
    cat("\nNumber of Factors: ")
    cat(format(length(unique(model$fmod))))
    cat("\n")
    cat("\nNumber of Bootstrap repeats: ")
    cat(format(model$nboot))
    cat("\n")
    if (model$modelo == 1) {
        cat("\nBandwidth: ")
        cat(format(model$h))
        cat("\n")
        cat("\nKernel: ")
        if (model$kernel == 1) 
            cat("Epanechnikov")
        if (model$kernel == 2) 
            cat("Triangular")
        if (model$kernel == 3) 
            cat("Gaussian")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
