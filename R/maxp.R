#' Maximum points for the estimate, first and second derivative, 
#' with their 95\% confidence intervals
#'@description Value of covariate \code{x} which maximizes the  estimate, 
#'first and second derivative, for each level of the factor.
#'
#'@param model Parametric or nonparametric regression out 
#' obtained by \code{\link{frfast}} function.
#'@param der Number which determines any inference process. By default
#' \code{der} is \code{NULL}. If this term is \code{0}, the calculation of the 
#' maximum point is for the estimate. If it is \code{1} or \code{2}, it is 
#' designed for the first or second derivative, respectively.
#' 
#'@return An object is returned with the following elements:
#'\item{Estimation}{outputs for estimation with maximum points and their
#' 95\% confidence intervals (for each level).}
#'\item{First_der}{outputs for first derivative with maximum points and their
#' 95\% confidence intervals (for each level).}
#' \item{Second_der}{outputs for second derivative with maximum points and their
#' 95\% confidence intervals (for each level).}
#' 
#'@author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#'
#'@examples
#'library(NPRegfast)
#'data(barnacle)
#'
#'fit <- frfast(DW ~ RC, data = barnacle) # without interactions
#'maxp(fit)
#'maxp(fit, der = 0)
#'maxp(fit, der = 1)
#'maxp(fit, der = 2)
#'
#'fit2 <- frfast(DW ~ RC : F,data = barnacle) # with interactions
#'maxp(fit2)
#'maxp(fit2, der = 0)
#'maxp(fit2, der = 1)
#'maxp(fit2, der = 2)




maxp <- function(model, der = NULL) {
  
  if(!is.null(der) & !isTRUE(der %in% c(0, 1, 2))) {
    stop("\"der\" is not a r-th derivative implemented")
  }
  
  nf <- model$nf
  jnf <- c()
  model$max[model$max == 9999] <- NA
  model$maxl[model$maxl == 9999] <- NA
  model$maxu[model$maxu == 9999] <- NA
  factores <- c()
  jnf <- c()
  if (is.null(der)) {
    res <- matrix(ncol = 3, nrow = nf)
    k <- 1
    for (i in 1:nf) {
      res[i, 1] <- c(model$max[k, i])
      res[i, 2] <- c(model$maxl[k, i])
      res[i, 3] <- c(model$maxu[k, i])
    }
    res2 <- matrix(ncol = 3, nrow = nf)
    k <- 2
    for (i in 1:nf) {
      res2[i, 1] <- c(model$max[k, i])
      res2[i, 2] <- c(model$maxl[k, i])
      res2[i, 3] <- c(model$maxu[k, i])
    }
    res3 <- matrix(ncol = 3, nrow = nf)
    k <- 3
    
    for (i in 1:nf) {
      res3[i, 1] <- c(model$max[k, i])
      res3[i, 2] <- c(model$maxl[k, i])
      res3[i, 3] <- c(model$maxu[k, i])
      jnf[i] <- which(model$label == model$label[i])
      factores[i] <- paste("Level", model$label[jnf[i]])
    }
    
    colnames(res) <- c("Max point", "Lwr", "Upr")
    colnames(res2) <- c("Max point", "Lwr", "Upr")
    colnames(res3) <- c("Max point", "Lwr", "Upr")
    rownames(res) <- c(factores)
    rownames(res2) <- c(factores)
    rownames(res3) <- c(factores)
    return(list(Estimation = res, First_der = res2, Second_der = res3))
  } else {
    der <- der + 1
    res <- matrix(ncol = 3, nrow = nf * length(der))
    k <- der
    a <- 1
    
    for (j in der) {
      for (i in 1:nf) {
        if (a == 2) {
          ii <- nf + i
        } else {
          ii <- i
        }
        res[ii, 1] <- c(model$max[j, i])
        res[ii, 2] <- c(model$maxl[j, i])
        res[ii, 3] <- c(model$maxu[j, i])
        jnf[i] <- which(model$label == model$label[i])
        factores[i] <- paste("Level", model$label[jnf[i]])
      }
      a <- 2
    }
    colnames(res) <- c("Max point", "Lwr", "Upr")
    rownames(res) <- c(rep(factores, length(der)))
    return(res)
  }
}
