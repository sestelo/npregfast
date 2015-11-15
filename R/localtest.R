#' Testing the equality of critical points
#'@description This function can be used to test the equality of the
#' \eqn{M} critical points estimated from the respective level-specific curves.
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
#' @param h0 The kernel bandwidth smoothing parameter for the global effect (see
#' references for more details at the estimation). Large values of the bandwidth lead
#' to smoothed estimates; smaller values of the bandwidth lead lo undersmoothed estimates. 
#' By default, cross validation is used to obtain the bandwidth.
#' @param h The kernel bandwidth smoothing parameter for the partial effects.
#' @param nh Integer number of equally-spaced bandwidth on which the
#' \code{h} is discretised, to speed up computation.
#' @param kernel A character string specifying the desired kernel. 
#' Defaults to \code{kernel = "epanech"}, where the Epanechnikov
#' density function kernel will be used. Also, several types of kernel funcitons 
#' can be used:  triangular and Gaussian density function, 
#' with \code{"triang"} and \code{"gaussian"} term, respectively.
#' @param p Degree of polynomial to be used. Its value must be the value of
#' derivative + 1. The default value is 3 due to the function
#' returns the estimation, first and second derivative.
#' @param kbin Number of binning nodes over which the function 
#' is to be estimated.
#' @param rankl Number or vector specifying the minimum value for the
#' interval at which to search the \code{x} value which maximizes the
#' estimate, first or second derivative  (for each level). The default
#' is the minimum data value.
#' @param ranku Number or vector specifying the maximum value for the
#' interval at which to search the \code{x} value which maximizes the
#' estimate, first or second derivative  (for each level). The default
#' is the maximum data value.
#' @param seed Seed to be used in the bootstrap procedure.
#' 
#' 
#' 
#' @details \code{localtest} can be used to test the equality of the 
#' \eqn{M} critical points estimated from the respective level-specific curves. 
#' Note that, even if the curves and/or their derivatives are different, it is 
#' possible for these points to be equal. 
#' 
#' For instance, taking the maxima of the first derivatives into account, 
#' interest lies in testing the following null hypothesis
#' 
#' \deqn{H_0: x_{01} = \ldots = x_{0M}}
#' 
#' versus the general  alternative  
#' 
#' \deqn{H_1: x_{0i} \ne x_{0j}  \quad {\rm{for}} \quad {\rm{some}} \quad 
#' \emph{i}, \emph{j} \in \{ 1, \ldots, M\}.}
#' 
#' The above hypothesis is true if \eqn{d=x_{0j}-x_{0k}=0} where 
#' \deqn{ (j,k)= argmax \quad (l,m) \quad \{1 \leq l<m \leq M\} \quad |x_{0l}-x_{0m}|, }
#' 
#' otherwise  \eqn{H_0} is false. It is important to highlight that, in practice,
#' the true \eqn{x_{0j}} are not known, and consequently neither is \eqn{d}, 
#' so an estimate \eqn{\hat d = \hat x_{0j}-\hat x_{0k}} is used, where, 
#' in general, \eqn{\hat x_{0l}} are the estimates of \eqn{x_{0l}} based on the 
#' estimated curves \eqn{\hat m_l} with \eqn{l = 1, \ldots , M}. 
#' 
#' Needless to say, 
#' since \eqn{\hat d} is only an estimate of the true \eqn{d}, the sampling 
#' uncertainty of these estimates needs to be acknowledged. Hence, a confidence 
#' interval \eqn{(a,b)} is created for \eqn{d} for a specific level of 
#' confidence (95\%).  Based on this, the null hypothesis is rejected if  
#' zero is not contained in the interval.
#' 
#' Note that if this hypothesis is rejected (and the factor has more than 
#' two levels), one option could be to use the \code{maxp.diff} function in 
#' order to obtain the differences between each pair of factor's levels.
#' 
#'@return The estimate of \eqn{d} value is returned and its confidence interval 
#'for a specific-level of confidence, i.e. 95\%. Additionally, it is shown 
#'the decision, accepted or rejected,  of the local test. Based on the null 
#'hypothesis is rejected if a zero value is not within the interval. 
#'
#'@author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#'
#' @references 
#' Sestelo, M. (2013). Development and computational implementation of 
#' estimation and inference methods in flexible regression models. 
#' Applications in Biology, Engineering and Environment. PhD Thesis, Department
#' of Statistics and O.R. University of Vigo.
#' 
#'@examples
#' library(npregfast)
#' data(barnacle)
#' localtest(DW ~ RC : F, data = barnacle, der = 1, seed = 130853, nboot = 100)
#' 
#' @useDynLib npregfast localtest_
#' @importFrom stats na.omit
#' @export




localtest <- function(formula, data = data, der, weights = NULL, 
                      nboot = 500, h0 = -1.0, h = -1.0, nh = 30, kernel = "epanech", 
                      p = 3, kbin = 100, rankl = NULL, ranku = NULL, seed = NULL) {
  
  if(kernel == "gaussian")  kernel <- 3
  if(kernel == "epanech")   kernel <- 1
  if(kernel == "triang")    kernel <- 2
  
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
  
 
  
  ncmax <- 5
  c2 <- NULL
  if(is.null(seed)) seed <- -1
  
  ffr <- interpret.frfastformula(formula, method = "frfast")
  varnames <- ffr$II[2, ]
  aux <- unlist(strsplit(varnames,split=":"))
  varnames <- aux[1]
  namef <- aux[2]
  if(length(aux) == 1){f <- NULL}else{f <- data[ ,namef]}
  newdata <- data
  data <- na.omit(data[ ,c(ffr$response, varnames)])
  newdata <- na.omit(newdata[ ,varnames])
  n <- nrow(data)
  
  
  
  
  if (is.null(f)) f <- rep(1, n)
  etiquetas <- unique(f)
  nf <- length(etiquetas)
  
  
  if(nf == 1) {
    stop("Function not supported.
         There is not factor in the model.")
  }
  
  if(is.null(h0)){
    h0 <- -1.0
  }
  if(is.null(h)){
    h <- rep(-1.0, nf)
  }else{
    if(length(h) == 1) h <- rep(h, nf)
  }
  
  # Interesaria meter para las derivadas?
  
  if (is.null(weights)) {
    weights <- rep(1, n)
  } else {
    if (sum(weights) <= 0 || any(weights) < 0 || length(weights) != n) 
      stop("The specified weights are not correct")
  }
  
  
  if(is.null(c2)) c2 <- matrix(as.double(-1.0), ncmax, nf) 
  if(is.null(rankl)){
    rankl <- na.omit(as.vector(tapply(data[ ,varnames], f, min)))
  }else{
    if(length(rankl) == 1) rankl <- rep(rankl, nf)
  }
  if(is.null(ranku)){
    ranku <- na.omit(as.vector(tapply(data[ ,varnames], f, max)))
  }else{
    if(length(ranku) == 1) ranku <- rep(ranku, nf)
  } 
  
  localtest  <-.Fortran("localtest_",
                        f = as.integer(f),
                        x = as.double(data[,varnames]),
                        y = as.double(data[,ffr$response]),
                        w = as.double(weights),
                        n = as.integer(n),
                        h0 = as.double(h0),
                        h = as.double(h),
                        nh = as.integer(nh),
                        p = as.integer(p),
                        kbin = as.integer(kbin),
                        #fact = as.integer(c(1:nf)),
                        fact = unique(as.integer(f)),
                        #fact   =as.integer(c(1:nf))
                        nf = as.integer(nf),
                        kernel = as.integer(kernel),
                        nboot = as.integer(nboot),
                        pcmax = as.double(ranku), # rango de busqueda minimo
                        pcmin =as.double(rankl), # rango de busqueda maximo
                        r = as.integer(der),
                        D = as.double(rep(-1.0,1)),
                        Ci = as.double(rep(-1.0,1)),
                        Cs = as.double(rep(-1.0,1)),
                        seed = as.integer(seed)
  )
  

  
  if (localtest$Ci <= 0 & 0 <= localtest$Cs) {
    decision <- "Acepted"
  } else {
    decision <- "Rejected"
  }
  res <- cbind(d = round(localtest$D, digits = 4), Lwr = round(localtest$Ci, digits = 4), 
               Upr = round(localtest$Cs, digits = 4), Decision = decision)
  # class(res) <- 'localtest'
  return(as.data.frame(res))
  
} 
