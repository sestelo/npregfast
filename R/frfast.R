#' Fitting nonparametric models
#' 
#' @description This function is used to fit nonparametric models by
#' using local linear kernel smoothers. Additionally, a parametric 
#' model (allometric model) can be estimated.
#' @param formula An object of class \code{formula}: a sympbolic 
#' description of the model to be fitted. The details of model 
#' specification are given under 'Details'.
#' @param data A data frame or matrix containing the model response
#' variable and covariates required by the \code{formula}.
#' @param model Type model used: \code{model = "np"}  nonparametric
#' regression model with local linear kernel smoothers, 
#' \code{model = "allo"} the  allometric model.
#' @param h0 The kernel bandwidth smoothing parameter for the global effect. 
#' Large values of bandwidth make smoother estimates, smaller values of 
#' bandwidth make less smooth estimates. The default is a bandwidth compute by 
#' cross validation.
#' @param h The kernel bandwidth smoothing parameter for the partial effects. 
#' @param nh Integer number of equally-spaced bandwidth on which the
#' \code{h} is discretised, to speed up computation.
#' @param weights Prior weights on the data.
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
#' @param nboot Number of bootstrap repeats. Default 500. The wild bootstrap
#' is used when \code{model = "np"} and the simple bootstrap when
#' \code{model = "allo"}.
#' @param rankl Number or vector specifying the minimum value for the
#' interval at which to search the \code{x} value which maximizes the
#' estimate, first or second derivative  (for each level). The default
#' is the minimum data value.
#' @param ranku Number or vector specifying the maximum value for the
#' interval at which to search the \code{x} value which maximizes the
#' estimate, first or second derivative  (for each level). The default
#' is the maximum data value.
#' @param seed Seed to be used in the bootstrap procedure.
#' @details The models fit by \code{frfast} function are specified 
#' in a compact symbolic form. The \~ operator is basic in the formation 
#' of such models. An expression of the form \code{y ~ model}  is interpreted as 
#' a specification that the response \code{y} is modelled by a predictor 
#' specified symbolically by \code{model}. The possible terms consist of a 
#' variable name or a variable name and a factor name separated by : operator. 
#' Such a term is interpreted as the interaction of the continuous variable and 
#' the factor.
#' @return An object is returned with the following elements:
#' \item{x}{Vector of values of the grid points at which model is to 
#' be estimate.}
#' \item{p}{Matrix of values of the grid points at which to compute the 
#' estimate, their first and second derivative.}
#' \item{pl}{Lower values of  95\% confidence interval for the estimate, 
#' their first and second derivative.}
#' \item{pu}{Upper values of  95\% confidence interval for the estimate, 
#' their first and second derivative.}
#' \item{diff}{Differences between the estimation values of a couple of 
#' levels (i. e. level 2 - level 1). The same procedure for their first
#' and second derivative.}
#' \item{diffl}{Lower values of 95\% confidence interval for the differences 
#' between the estimation values of a couple of levels. It is performed 
#' for their first and second derivative.}
#' \item{diffu}{Upper values of 95\% confidence interval for the differences 
#' between the estimation values of a couple of levels. It is performed for 
#' their first and second derivative.}
#' \item{nboot}{Number of bootstrap repeats.}
#' \item{n}{Total number of data}
#' \item{dp}{Degree of polynomial used.}
#' \item{h0}{The kernel bandwidth smoothing parameter for the global effect.}
#' \item{h}{The kernel bandwidth smoothing parameter for the partial effects.}
#' \item{fmod}{Factor's level for each data.}
#' \item{xdata}{Original x values.}
#' \item{ydata}{Original y values.}
#' \item{w}{Weights on the data.}
#' \item{kbin}{Number of binning nodes over which the function is to 
#' be estimated.}
#' \item{nf}{Number of levels.}
#' \item{max}{Value of covariate \code{x} which maximizes the  estimate, 
#' first or second derivative.}
#' \item{maxu}{Upper value of 95\% confidence interval for the 
#' value \code{max}.}
#' \item{maxl}{Lower value of 95\% confidence interval for the 
#' value \code{max}.}
#' \item{diffmax}{Differences between the estimation of \code{max} for a 
#' couple of levels (i. e. level 2 - level 1). The same procedure for their 
#' first and second derivative.}
#' \item{diffmaxu}{Upper value of 95\% confidence interval for the value 
#' \code{diffmax}.}
#' \item{diffmaxl}{Lower value of 95\% confidence interval for the value 
#' \code{diffmax}.}
#' \item{repboot}{Matrix of values of the grid points at which to compute 
#' the estimate, their first and second derivative for each bootstrap repeat.}
#' \item{rankl}{Maximum value for the interval at which to search the 
#' \code{x} value which maximizes the estimate, first or second derivative  
#' (for each level). The default is the maximum data value.}
#' \item{ranku}{Minimum value for the interval at which to search the 
#' \code{x} value which maximizes the estimate, first or second derivative  
#' (for each level). The default is the minimum data value.}
#' \item{nmodel}{Type model used: \code{model = 1} the nonparametric model, 
#' \code{model = 2} the allometric model.}
#' \item{label}{Labels of the variables in the model.}
#' \item{numlabel}{Number of labels.}
#' \item{kernel}{Character which determines the smoothing kernel.}
#' \item{a}{Estimated coefficient in the case of fitting an allometric model.}
#' \item{al}{Lower value of 95\% confidence interval for the value of \code{a}.}
#' \item{au}{Upper value of 95\% confidence interval for the value of \code{a}.}
#' \item{b}{Estimated coefficient in the case of fitting an allometric model.}
#' \item{bl}{Lower value of 95\% confidence interval for the value of \code{b}.}
#' \item{bu}{Upper value of 95\% confidence interval for the value of \code{b}.}
#' \item{name}{Name of the variables in the model.}
#' \item{formula}{A sympbolic description of the model to be fitted.}
#' \item{nh}{Integer number of equally-spaced bandwidth on which the
#' \code{h} is discretised.}
#' \item{r2}{Coefficient of determination.}
#' 
#' 
#' @author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#' 
#' @examples
#' library(NPRegfast)
#' data(barnacle)
#' 
#' # Nonparametric regression without interactions
#' fit <- frfast(DW ~ RC, data = barnacle) 
#' fit
#' summary(fit)
#' 
#' # Change the number of binning nodes and bootstrap replicates
#' fit <- frfast(DW ~ RC, data = barnacle, kbin = 200, nboot = 1000)
#' 
#' # Nonparametric regression with interactions
#' fit2 <- frfast(DW ~ RC : F, data = barnacle)
#' fit2
#' summary(fit2)


frfast <- function(formula, data = data, model = "np", h0 = -1.0, h = -1.0, 
                   nh = 30, weights = NULL, kernel = "epanech", p = 3, 
                   kbin = 100, nboot = 500, rankl = NULL, ranku = NULL, 
                   seed = NULL){
  
  if(kernel == "gaussian")  kernel <- 3
  if(kernel == "epanech")   kernel <- 1
  if(kernel == "triang")    kernel <- 2
  
  if(missing(formula)){
    stop("Argument \"formula\" is missing, with no default")
  }
  if(missing(data)){
    stop("Argument \"data\" is missing, with no default")
  }
  if(!(kernel %in% 1:3)){
    stop("Kernel not suported")
  }
  
  
  ncmax <- 5
  c2 <- NULL
  
  
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
  
  if(is.null(f)) f <- rep(1.0, n)
  etiquetas <- unique(f)
  nf <- length(etiquetas)
  
  if(model == "np") tmodel <- 1
  if(model == "allo") tmodel <- 2
  
  
  
  if(is.null(h0)){
    h0 <- -1.0
  }
  if(is.null(h)){
    h <- rep(-1.0, nf)
  }else{
    if(length(h) == 1) h <- rep(h, nf)
  }
  
  
  if(is.null(weights)) {
    weights <- rep(1.0, n)
  }else{
    if(sum(weights) <= 0 || any(weights < 0) || length(weights) != n) 
      stop("The specified weights are not correct")
  }  
  
  
  if(is.null(c2)) c2 <- matrix(as.double(-1.0), ncmax, nf) 
  if(is.null(rankl)){
    rankl <- as.vector(tapply(data[ ,varnames], f, min))
  }else{
    if(length(rankl) == 1) rankl <- rep(rankl, nf)
    }
  if(is.null(ranku)){
    ranku <- as.vector(tapply(data[ ,varnames], f, max))
  }else{
    if(length(ranku) == 1) ranku <- rep(ranku, nf)
  } 
  
  ipredict2 <- 0
  
  frfast  <- .Fortran("frfast",
                      f = as.integer(f),
                      x = as.double(data[ ,varnames]),
                      y = as.double(data[ ,ffr$response]),
                      w = as.double(weights),
                      n = as.integer(n),
                      h0 = as.double(h0),
                      h = as.double(h),
                      c2 = as.integer(c2),
                      ncmax = as.integer(ncmax),
                      p = as.integer(p),
                      kbin = as.integer(kbin),
                      fact = as.integer(c(1:nf)), 
                      nf = as.integer(nf),
                      nboot = as.integer(nboot),
                      xb = as.double(rep(-1.0, kbin)),
                      pb = array(rep(-1.0), c(kbin, 3, nf)),
                      li = array(as.double(-1.0), c(kbin, 3, nf)),
                      ls = array(as.double(-1.0), c(kbin, 3, nf)),
                      dif = array(as.double(-1.0), c(kbin, 3, nf, nf)),
                      difi = array(as.double(-1.0), c(kbin, 3, nf, nf)),
                      difs = array(as.double(-1.0), c(kbin, 3, nf, nf)),
                      tmodel = as.integer(tmodel), 
                      c = array(as.double(-1.0), c(3, nf)),
                      cs = array(as.double(-1.0), c(3, nf)),
                      ci = array(as.double(-1.0), c(3, nf)),
                      difc = array(as.double(-1.0), c(3, nf, nf)),
                      difcs = array(as.double(-1.0), c(3, nf, nf)),
                      difci = array(as.double(-1.0), c(3, nf, nf)),
                      pboot = array(as.double(-1.0), c(kbin, 3, nf, nboot)),
                      pcmin = as.double(rankl),
                      pcmax = as.double(ranku), 
                      cboot = array(as.double(-1.0), c(3, nf, nboot)), 
                      kernel = as.integer(kernel),
                      nh = as.integer(nh),
                      a = as.double(rep(-1, nf)),
                      ainf = as.double(rep(-1, nf)),
                      asup = as.double(rep(-1, nf)),
                      b = as.double(rep(-1, nf)),
                      binf = as.double(rep(-1, nf)),
                      bsup = as.double(rep(-1, nf)),
                      ipredict = as.integer(ipredict2),
                      predict = array(rep(-1.0), c(kbin, 3, nf)),
                      predictl = array(as.double(-1.0), c(kbin, 3, nf)),
                      predictu = array(as.double(-1.0), c(kbin, 3, nf))
  )
  
  if(tmodel != 2){
    frfast$a <- NULL
    frfast$ainf <- NULL
    frfast$asup <- NULL
    frfast$b <- NULL
    frfast$binf <- NULL
    frfast$bsup <- NULL
    r2 <- NULL
  }
  
  #R-squared
  if(tmodel == 2){
    yhat <- frfast$a * (frfast$x^frfast$b)
    rss <- sum( (frfast$y-yhat)**2 ) / (frfast$n-2)
    tss <- sum(  (frfast$y-mean(frfast$y))**2 ) / (frfast$n-1)
    r2 <- 1-(rss/tss)
  }
  
  frfast$pb[frfast$pb == -1] <- NA
  frfast$li[frfast$li == -1] <- NA
  frfast$ls[frfast$ls == -1] <- NA
  frfast$dif[frfast$dif == -1] <- NA
  frfast$difi[frfast$difi == -1] <- NA
  frfast$difs[frfast$difs == -1] <- NA
  frfast$c[frfast$c == -1] <- NA
  frfast$cs[frfast$cs == -1] <- NA
  frfast$ci[frfast$ci == -1] <- NA
  frfast$difc[frfast$difc == -1] <- NA
  frfast$difcs[frfast$difcs == -1] <- NA
  frfast$difci[frfast$difci == -1] <- NA
  frfast$pboot[frfast$pboot == -1] <- NA
  
  res <- list(x = frfast$xb,
              p = frfast$pb, 
              pl = frfast$li,
              pu = frfast$ls,
              diff = frfast$dif,
              diffl = frfast$difi,
              diffu = frfast$difs,
              nboot = frfast$nboot,
              n = frfast$n,
              dp = frfast$p,
              h = frfast$h,
              h0 = frfast$h0,
              fmod = frfast$f,
              xdata = as.vector(frfast$x),
              ydata = frfast$y,
              w = frfast$w,
              #fact=fact,  # Lo tuve que comentar pq me daba error
              kbin = frfast$kbin,
              nf = frfast$nf,
              max = frfast$c, #
              maxu = frfast$cs, 
              maxl = frfast$ci,
              diffmax = frfast$difc,
              diffmaxu = frfast$difcs,
              diffmaxl = frfast$difci,
              repboot = frfast$pboot,  
              rankl = frfast$pcmin,
              ranku = frfast$pcmax,
              nmodel = frfast$tmodel, 
              label = as.character(etiquetas),
              numlabel = unique(frfast$f),
              kernel = frfast$kernel, 
              a = frfast$a,
              al = frfast$ainf,
              au = frfast$asup,
              b = frfast$b,
              bl = frfast$binf,
              bu = frfast$bsup,
              name = c(ffr$response,varnames),
              formula = formula,
              nh = frfast$nh,
              r2 = r2,
              call = match.call()
  )
  
  class(res) <- "frfast"
  return(res)
}


