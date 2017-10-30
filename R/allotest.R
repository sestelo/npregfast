#' Bootstrap based test for testing an allometric model
#'@description Bootstrap-based procedure that tests whether the data 
#' can be modelled by an allometric model.
#'@param formula An object of class \code{formula}: a sympbolic description
#' of the model to be fitted.
#' @param data An optional data frame, matrix or list required by 
#' the formula. If not found in data, the variables are taken from 
#' \code{environment(formula)}, typically the environment from which
#'  \code{frfast} is called.
#' @param na.action A function which indicates what should happen when the 
#' data contain 'NA's. The default is 'na.omit'.
#'@param nboot Number of bootstrap repeats.
#'@param seed Seed to be used in the bootstrap procedure.
#' @param cluster A logical value. If  \code{TRUE} (default), the
#'  bootstrap procedure is  parallelized (only for \code{smooth = "splines"}).
#'   Note that there are cases 
#'  (e.g., a low number of bootstrap repetitions) that R will gain in
#'  performance through serial computation. R takes time to distribute tasks
#'  across the processors also it will need time for binding them all together
#'  later on. Therefore, if the time for distributing and gathering pieces
#'  together is greater than the time need for single-thread computing, it does
#'  not worth parallelize.
#'@param ncores An integer value specifying the number of cores to be used
#' in the parallelized procedure. If \code{NULL} (default), the number of cores 
#' to be used is equal to the number of cores of the machine - 1.
#'@param test Statistic test to be used, based on residuals on the null model
#'  (\code{res}) or based on the likelihood ratio test 
#'  using rss0 and rss1 \code{lrt}.
#' @param \ldots Other options.
#'  
#'  
#' 
#'@details In order to facilitate the choice of a model appropriate
#' to the data while at the same time endeavouring to minimise the 
#' loss of information,  a bootstrap-based procedure, that test whether the 
#' data can be modelled by an allometric model, was developed.  Therefore,
#' \code{allotest} tests the null hypothesis of an allometric model taking 
#' into account the logarithm of the original variable
#'  (\eqn{X^* = log(X)} and \eqn{Y^* = log(Y)}). 
#'  
#' Based on a general model of the type 
#' \deqn{Y^*=m(X^*)+\varepsilon}
#' the aim here is to test the null hypothesis of an allometric model 
#' \deqn{H_0 = m(x^*) =  a^*+ b^* x^*}
#' \eqn{vs.} the general hypothesis 
#' \eqn{H_1}, with \eqn{m}
#' being an unknown nonparametric function; or analogously,
#' \deqn{H_1: m(x^*)= a^*+ b^* x^* + g(x^*)}
#' with \eqn{g(x^*)} being an unknown function not equal to zero. 
#' 
#' To implement this test we have used the wild bootstrap.
#' 
#' 
#' 
#'@return An object is returned with the following elements:
#' \item{statistic}{the value of the test statistic.}
#' \item{value}{the p-value of the test.}
#' 
#'@author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#'
#'@references 
#' Sestelo, M. and Roca-Pardinas, J. (2011). A new approach to estimation of 
#' length-weight relationship of \eqn{Pollicipes}  \eqn{pollicipes} (Gmelin, 1789)
#' on the Atlantic coast of Galicia (Northwest Spain): some aspects of its 
#' biology and management. Journal of Shellfish Research, 30 (3), 939--948.
#' 
#' Sestelo, M. (2013). Development and computational implementation of 
#' estimation and inference methods in flexible regression models. 
#' Applications in Biology, Engineering and Environment. PhD Thesis, Department
#' of Statistics and O.R. University of Vigo.
#' 
#' 
#'@examples
#' library(npregfast)
#' data(barnacle)
#' allotest(DW ~ RC, data = barnacle, nboot = 50, seed = 130853)
#' 
#' @useDynLib npregfast allotest_
#' @importFrom stats na.omit runif
#' @export





allotest <- function(formula, data, na.action = "na.omit",
                     nboot = 500, seed = NULL, cluster = TRUE,
                     ncores = NULL, test = "res", ...) {
  
  if (isTRUE(cluster)) {
    if (is.null(ncores)) {
      num_cores <- detectCores() - 1
    }else{
      num_cores <- ncores
    }
    registerDoParallel(cores = num_cores)
  }
  
  
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(x = c("formula", "data", "subset", "weights", "na.action", "offset"), 
             table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  
  
  
  
  
  mf <- eval(expr = mf, envir = parent.frame())
  
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w)) 
    stop("'weights' must be a numeric vector")
  
  terms <- attr(mt, "term.labels")
  
  
  
  aux <- unlist(strsplit(terms,split = ":"))
  varnames <- aux[1]
  namef <- aux[2]
  response <- as.character(attr(mt, "variables")[2])
  if (unlist(strsplit(varnames,split = ""))[1] == "s") {
    stop("Argument \"formula\" is wrong specified, see details of
         model specification in 'Details' of the frfast help." )
  }
  
  
  #newdata <- data
  data <- mf
  
  if (na.action == "na.omit"){ # ver la f, corregido
    data <- na.omit(data)
  }else{
    stop("The actual version of the package only supports 'na.omit' (observations are removed 
         if they contain any missing values)")
  }
  
  if (length(aux) == 1) {f <- NULL}else{f <- data[ ,namef]}
  n <- nrow(data)
  
  
  
  
  
  
  
  
  
  
  # 
  # 
  # 
  # ffr <- interpret.frfastformula(formula, method = "frfast")
  # varnames <- ffr$II[2, ]
  # aux <- unlist(strsplit(varnames, split = ":"))
  # varnames <- aux[1]
  # namef <- aux[2]
  # if (length(aux) == 1) {
  #   f <- NULL
  # } else {
  #   f <- data[, namef]
  # }
  # newdata <- data
  # data <- data[, c(ffr$response, varnames)]
  # newdata <- newdata[, varnames]
  # 
  # if (na.action == "na.omit") { # ver la f
  #   data <- na.omit(data)
  #   newdata <- na.omit(newdata)
  # }else{
  #   stop("The actual version of the package only supports 'na.omit' (observations are removed 
  #        if they contain any missing values)")
  # }
  # 
  # 
  # n <- nrow(data)
  # 
  # #if (is.null(seed)) {
  # #  set.seed(NULL)
  # #  seed <- .Random.seed[3]
  # #}
  
  if (!is.null(seed)){
    set.seed(seed)
  }
  
  #umatrix <- matrix(runif(n*nboot), ncol = nboot, nrow = n)
  
  if (is.null(f)) 
    f <- rep(1, n)
  etiquetas <- unique(f)
  nf <- length(etiquetas)
  
  
  res <- list()
  
  for (i in etiquetas) {
    yy <- data[, 1][f == i]
    xx <- data[, 2][f == i]
    yy[yy == 0] <- 0.0001
    xx[xx == 0] <- 0.0001
    n <- length(xx)
    w <- rep(1, n)
    
    m <- lm(log(yy) ~ log(xx))
    muhatg <- exp(coef(m)[1]) * xx**coef(m)[2]
    errg <- yy - muhatg
    
   if(test == "res") {t <- sta_res(x = xx, y = yy)}
   if(test == "lrt") {t <- sta_rss(x = xx, y = yy)}
    #print(c(t1, t2))
    
    yboot <- replicate(nboot, muhatg + errg * 
                         sample(c(-sqrt(5) + 1, sqrt(5) + 1)/2, size = n, 
                                replace = TRUE, 
                                prob = c(sqrt(5) + 1, sqrt(5) - 1)/(2 * sqrt(5))))
    
    if(test == "res") {
      tboot <- unlist(foreach(i = 1:nboot) %dopar% {
      sta_res(x = xx, y = yboot[, i])
    })
    }
    
    if(test == "lrt") {
    tboot <- unlist(foreach(i = 1:nboot) %dopar% {
      sta_rss(x = xx, y = yboot[, i])
    })
    }
    pvalue <- mean(tboot>t)


    res[[i]] <- list(statistic = c(t), pvalue = c(pvalue))
    
  }
  
  
  est <- c()
  p <- c()
  for (i in 1:length(res)) {
    est[i] <- res[[i]]$statistic
    p[i] <- res[[i]]$pvalue
  }
  
  kk <- cbind(est, p)
  result <- matrix(kk, ncol = 2, nrow = length(res))
  colnames(result) <- c("Statistic", "pvalue")
  
  
  # factores=paste('Factor',1:length(res))
  factores <- paste("Level", etiquetas[1:length(etiquetas)])
  
  
  for (i in 1:length(res)) {
    rownames(result) <- c(factores)
  }
  
  return(result)
  
} 






sta_res <- function(x, y){
  y[y == 0] <- 0.0001
  model <- lm(log(y) ~ log(x))
  muhat <- exp(coef(model)[1]) * x**coef(model)[2]
  residuo <- y - muhat
  pred <- as.numeric(predict(gam(residuo ~ s(x))))
  #pred <- predict(frfast(residuo ~ x, data = data.frame(x, residuo), nboot = 0), newdata = data.frame(x=x))$Estimation[,1]
  #pred <- pred - mean(pred)
  #rango <- max(x) - min(x)
  #ii <- abs(x) <= (max(x)-0.1*rango)
  t <- sum(abs(pred))
}

sta_rss <- function(x, y){
  y[y == 0] <- 0.0001
  model <- lm(log(y) ~ log(x))
  m0 <- exp(coef(model)[1]) * x**coef(model)[2]
  rss0 <- sum((y - m0)**2)
  m1 <- as.numeric(predict(gam(y ~ s(x))))
  rss1 <- sum((y - m1)**2)
  #rango <- max(x) - min(x)
  #ii <- abs(x) <= (max(x)-0.1*rango)
  t <- (rss0 - rss1)/rss1
}






