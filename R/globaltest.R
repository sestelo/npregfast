#' Testing the equality of M curves specific to each level
#'@description \code{globaltest}  can be used to test the equality of the 
#'\eqn{M} curves specific to each level.
#'@param formula An object of class \code{formula}: a sympbolic description
#' of the model to be fitted.
#'@param data A data frame or matrix containing the model response variable
#' and covariates required by the \code{formula}.
#' @param der Number which determines any inference process. 
#' By default \code{der} is \code{NULL}. If this term is \code{0}, 
#' the testing procedures is applied for the estimate. If it is \code{1} or
#' \code{2}, it is designed for the first or second derivative, respectively.
#' @param weights Prior weights on the data.
#' @param nboot Number of bootstrap repeats.
#' @param h0 The kernel bandwidth smoothing parameter for the global effect. 
#' Large values of bandwidth make smoother estimates, smaller values of 
#' bandwidth make less smooth estimates. The default is a bandwidth compute by 
#' cross validation.
#' @param h The kernel bandwidth smoothing parameter for the partial effects.
#' @param nh Integer number of equally-spaced bandwidth on which the
#' \code{h} is discretised, to speed up computation.
#' @param kernel Character which determines the smoothing kernel. 
#' By default \code{kernel = "epanech"}, this is, the Epanechnikov
#' density function. Also, several types of kernel funcitons 
#' can be used:  triangular and Gaussian density function, 
#' with \code{"triang"} and \code{"gaussian"} term, respectively.
#' @param p Degree of polynomial used.  Its value must be the value of
#' derivative + 1. The default value is 3 due to the function
#' returns the estimation, first and second derivative.
#' @param kbin Number of binning nodes over which the function 
#' is to be estimated.
#' @param seed Seed to be used in the bootstrap procedure.
#' 
#' 
#' 
#' @details \code{globaltest}  can be used to test the equality of the \eqn{M} 
#' curves specific to each level. This bootstrap based test assumes the  
#' following null hypothesis
#' 
#' \eqn{H_0: m_1 = \ldots = m_M}
#' 
#' Note that, if \eqn{H_0} is not rejected, then the equality of critical points
#' will also accepted. To test the null hypothesis, it is used an statistic, 
#' \eqn{T}, based on direct nonparametric estimates of the curves. 
#' If the null hypothesis is true, the \eqn{T} value should be close to zero 
#' but is generally greater. The test rule based on \eqn{T} consists of 
#' rejecting the null hypothesis if \eqn{T > T^{1- \alpha}}, where \eqn{T^p} 
#' is the empirical \eqn{p}-percentile of \eqn{T} under the null hypothesis. 
#' 
#'@return The \eqn{T} value and the \eqn{p}-value  are returned. Additionally, 
#'it is shown the decision, accepted or rejected, of the global test. 
#'The null hypothesis is rejected if the \eqn{p}-value\eqn{< 0.05}.   
#'
#'@author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#'@examples
#' library(NPRegfast)
#' data(barnacle)
#' globaltest(DW ~ RC : F, data = barnacle, der = 1, seed = 130853)
#' 
#' @export



globaltest <- function(formula, data = data, der, weights = NULL, nboot = 200,
                       h0 = -1, h = -1, nh = 30, kernel = "epanech", p = 3, 
                       kbin = 100, seed = NULL) {
  
  if (missing(der)) {
    stop("Argument \"der\" is missing, with no default")
  }
  
  if (missing(formula)) {
    stop("Argument \"formula\" is missing, with no default")
  }
  if (missing(data)) {
    stop("Argument \"data\" is missing, with no default")
  }
  
  if(!isTRUE(der %in% c(0, 1, 2))) {
    stop("",paste(der)," is not a r-th derivative implemented, only 
         permitted 0, 1 or 2.")
  }
  
  
  if (kernel == "gaussian") 
    kernel <- 3
  if (kernel == "epanech") 
    kernel <- 1
  if (kernel == "triang") 
    kernel <- 2
  
  
  if (is.null(seed)) 
    seed <- -1
  
  
  ffr <- interpret.frfastformula(formula, method = "frfast")
  varnames <- ffr$II[2, ]
  aux <- unlist(strsplit(varnames, split = ":"))
  varnames <- aux[1]
  namef <- aux[2]
  if (length(aux) == 1) {
    f <- NULL
  } else {
    f <- data[, namef]
  }
  newdata <- data
  data <- na.omit(data[, c(ffr$response, varnames)])
  newdata <- na.omit(newdata[, varnames])
  n <- nrow(data)
  
  if (is.null(f)) 
    f <- rep(1, n)
  etiquetas <- unique(f)
  nf <- length(etiquetas)
  
  if(nf == 1) {
    stop("Function not supported.
         There is not factor in the model.")
  }
  
  if (is.null(h0)) {
    h0 <- -1
  }
  if (is.null(h)) {
    h <- rep(-1, nf)
  } else {
    if (length(h) == 1) 
      h <- rep(h, nf)
  }
  
  
  
  if (is.null(weights)) {
    weights <- rep(1, n)
  } else {
    if (sum(weights) <= 0 || any(weights) < 0 || length(weights) != n) 
      stop("The specified weights are not correct")
  }
  
  globaltest <- .Fortran("globaltest", 
                         f = as.integer(f), 
                         x = as.double(data[, varnames]), 
                         y = as.double(data[, ffr$response]), 
                         w = as.double(weights), 
                         n = as.integer(n), 
                         h0 = as.double(h0), 
                         h = as.double(h), 
                         nh = as.integer(nh), 
                         p = as.integer(p), 
                         kbin = as.integer(kbin), 
                         fact = as.integer(c(1:nf)), 
                         nf = as.integer(nf), 
                         kernel = as.integer(kernel), 
                         nboot = as.integer(nboot), 
                         r = as.integer(der), 
                         T = as.double(rep(-1, 1)), 
                         pvalor = as.double(rep(-1, 1)), 
                         seed = as.integer(seed)
                         )
 
  if (globaltest$pvalor < 0.05) {
    decision <- "Rejected"
  } else {
    decision <- "Acepted"
  }
  res <- data.frame(cbind(Statistic = globaltest$T, pvalue = globaltest$pvalor), 
                    Decision = I(decision))
  # res=cbind(Statistic=round(globaltest$T,digits=4),pvalue=round(globaltest$pvalor,digits=4),Decision=I(decision))
  # res=as.numeric(res) res=as.data.frame(res) class(res) <- 'globaltest'
  return(res)
  
} 