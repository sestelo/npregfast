#' Differences between the estimation of maximum points  
#' for two factor's levels
#'@description Differences between the estimation of \code{\link{maxp}} for two 
#' factor's levels. \code{maxp}, a returned element  of class 
#' \code{\link{frfast}}, is the value of covariate \code{x} which maximizes 
#' the  estimate, first or second derivative.
#'@param model Parametric or nonparametric regression model 
#' obtained by \code{\link{frfast}} function.
#' @param factor1 First factor's level at which to perform the differences 
#' between maximum points.
#' @param factor2 Second factor's level at which to perform the differences 
#' between maximum points.
#' @param der Number which determines any inference process. By default 
#' \code{der} is \code{NULL}. If this term is \code{0}, the calculate of the 
#' differences for maximum point is for the estimate. If it is \code{1} or 
#' \code{2}, it is designed for the first or second derivative, respectively.
#' 
#'@details Differences are calculated by subtracting a factor relative to 
#'another (\eqn{factor2 - factor1}).  By default \code{factor2} and 
#'\code{factor1} are \code{NULL}, so the differences calculated are for all 
#'possible combinations between two factors.
#' 
#'@return An object is returned with the following elements:
#' \item{maxp.diff}{a table with a couple of factor's level where it is used 
#' to calculate the differences between maximum points, and their 
#' 95\% interval confidence (for the estimation, first and second derivative).}
#' 
#'@author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#'@examples
#' library(NPRegfast)
#' data(barnacle)
#' fit2 <- frfast(DW ~ RC : F, data = barnacle) # with interactions
#' maxp.diff(fit2)
#' maxp.diff(fit2, der = 1)
#' maxp.diff(fit2, factor1="marta",factor2="nora")
#' 
#' @export

maxp.diff <- function(model, factor1 = NULL, factor2 = NULL, der = NULL) {
  
  if(model$nf == 1) {
    stop("There are not factors in the model.")
  }
  
  if(!is.null(factor1) & !isTRUE(factor1 %in% model$label)) {
    stop("\"",paste(factor1),"\" is not a factor's level.")
  }
  
  if(!is.null(factor2) & !isTRUE(factor2 %in% model$label)) {
    stop("\"",paste(factor2),"\" is not a factor's level.")
  }
  
  if(!is.null(der) & !isTRUE(der %in% c(0, 1, 2))) {
    stop("",paste(der)," is not a r-th derivative implemented, only 
         permitted 0, 1 or 2.")
  }
  
  
  nf <- model$nf
  model$diffmax[model$diffmax == 9999] <- NA
  model$diffmaxl[model$diffmaxl == 9999] <- NA
  model$diffmaxu[model$diffmaxu == 9999] <- NA
  
  a <- t(matrix(combn(nf, 2), nrow = 2))
  nrow(a)
  res <- list()
  
  if (is.null(der) & is.null(factor2) & is.null(factor1)) {
    
    for (i in 1:nrow(a)) {
      res[i] <- list(matrix(ncol = 5, nrow = 3))
      colnames(res[[i]]) <- c("Factor2", "Factor1", "Max point", "Lwr", "Upr")
      rownames(res[[i]]) <- c("Estimation", "First_der", "Second_der")
      for (k in 1:3) {
        res[[i]][k, 1] <- model$label[a[i, 2]]
        res[[i]][k, 2] <- model$label[a[i, 1]]
        res[[i]][k, 3] <- round(c(model$diffmax[k, a[i, 1], a[i, 2]]), 3)
        res[[i]][k, 4] <- round(c(model$diffmaxl[k, a[i, 1], a[i, 2]]), 3)
        res[[i]][k, 5] <- round(c(model$diffmaxu[k, a[i, 1], a[i, 2]]), 3)
      }
    }
    
    return(data.frame(res))
    
  } else if (is.null(der)) {
    
    factor2 <- which(model$label == factor2)
    factor1 <- which(model$label == factor1)
    
    
    res <- matrix(ncol = 5, nrow = 3)
    
    for (k in 1:3) {
      res[k, 1] <- model$label[factor2]
      res[k, 2] <- model$label[factor1]
      if (factor2 < factor1) {
        fac2 <- factor2
        fac1 <- factor1
        factor2 <- fac1
        factor1 <- fac2
      } else {
        fac2 <- factor2
        fac1 <- factor1
      }
      res[k, 3] <- if (fac2 < fac1) {
        -1 * round(c(model$diffmax[k, factor1, factor2]), 3)
      } else {
        round(c(model$diffmax[k, factor1, factor2]), 3)
      }
      res[k, 4] <- if (fac2 < fac1) {
        -1 * round(c(model$diffmaxl[k, factor1, factor2]), 3)
      } else {
        round(c(model$diffmaxl[k, factor1, factor2]), 3)
      }
      res[k, 5] <- if (fac2 < fac1) {
        -1 * round(c(model$diffmaxu[k, factor1, factor2]), 3)
      } else {
        round(c(model$diffmaxu[k, factor1, factor2]), 3)
      }
      if (fac2 < fac1) {
        factor2 <- fac2
        factor1 <- fac1
      }
    }
    
    colnames(res) <- c("Factor2", "Factor1", "Max points Diff.", "Lwr", "Upr")
    rownames(res) <- c("Estimation", "First_der", "Second_der")
    return(data.frame(res))
    
    
    
  } else if (is.null(factor2) & is.null(factor1)) {
    der <- der + 1
    for (i in 1:nrow(a)) {
      res[i] <- list(matrix(ncol = 5, nrow = 1))
      
      colnames(res[[i]]) <- c("Factor2", "Factor1", "Max points Diff.", "Lwr", 
                              "Upr")
      if (der == 1) 
        rownames(res[[i]]) <- c("Estimation")
      if (der == 2) 
        rownames(res[[i]]) <- c("First_der")
      if (der == 3) 
        rownames(res[[i]]) <- c("Second_der")
      res[[i]][1, 1] <- model$label[a[i, 2]]
      res[[i]][1, 2] <- model$label[a[i, 1]]
      res[[i]][1, 3] <- round(c(model$diffmax[der, a[i, 1], a[i, 2]]), 3)
      res[[i]][1, 4] <- round(c(model$diffmaxl[der, a[i, 1], a[i, 2]]), 3)
      res[[i]][1, 5] <- round(c(model$diffmaxu[der, a[i, 1], a[i, 2]]), 3)
    }
    return(data.frame(res))
    
  } else {
    
    factor2 <- which(model$label == factor2)
    factor1 <- which(model$label == factor1)
    
    der <- der + 1
    res <- matrix(ncol = 5, nrow = 1)
    res[1, 1] <- model$label[factor2]
    res[1, 2] <- model$label[factor1]
    if (factor2 < factor1) {
      fac2 <- factor2
      fac1 <- factor1
      factor2 <- fac1
      factor1 <- fac2
    } else {
      fac2 <- factor2
      fac1 <- factor1
    }
    res[1, 3] <- if (fac2 < fac1) {
      -1 * round(c(model$diffmax[der, factor1, factor2]), 3)
    } else {
      round(c(model$diffmax[der, factor1, factor2]), 3)
    }
    res[1, 4] <- if (fac2 < fac1) {
      -1 * round(c(model$diffmaxl[der, factor1, factor2]), 3)
    } else {
      round(c(model$diffmaxl[der, factor1, factor2]), 3)
    }
    res[1, 5] <- if (fac2 < fac1) {
      -1 * round(c(model$diffmaxu[der, factor1, factor2]), 3)
    } else {
      round(c(model$diffmaxu[der, factor1, factor2]), 3)
    }
    
    
    colnames(res) <- c("Factor2", "Factor1", "Max points Diff.", "Lwr", "Upr")
    if (der == 1) 
      rownames(res) <- c("Estimation")
    if (der == 2) 
      rownames(res) <- c("First_der")
    if (der == 3) 
      rownames(res) <- c("Second_der")
    return(data.frame(res))
    
  }
}