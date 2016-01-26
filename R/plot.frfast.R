#' Visualization of \code{frfast} objects
#' @description Useful for drawing the estimated regression function, 
#' first and second derivative (for each factor's level).
#' @param x \code{frfast} object.
#' @param y NULL.
#' @param fac Vector which determines the level to take into account
#' in the plot. By default is \code{NULL}.
#' @param der Number or vector which determines any inference process. 
#' By default \code{der} is \code{NULL}. If this term is \code{0}, the plot 
#' shows the initial estimate. If it is \code{1} or \code{2},
#' it is designed for the first or second derivative, respectively.
#' @param points Draw the original data into the plot. By default it is
#' \code{TRUE}.
#' @param xlab A title for the \code{x} axis. 
#' @param ylab A title for the \code{y} axis. 
#' @param ylim The \code{y} limits of the plot.
#' @param main An overall title for the plot.
#' @param col A specification for the default plotting color.
#' @param CIcol A specification for the default confidence intervals
#' plotting color.
#' @param ablinecol The color to be used for \code{abline}.
#' @param abline Draw an horizontal line into the plot of the second derivative 
#' of the model.
#' @param type What type of plot should be drawn. Possible types are,
#' \code{p} for points, \code{l} for lines, \code{o} for overplotted, etc. 
#' See details in \code{\link{par}}.
#' @param CItype What type of plot should be drawn for confidence intervals. 
#' Possible types are, \code{p} for points, \code{l} for lines, \code{o} 
#' for overplotted.
#' @param lwd The line width, a positive number, defaulting to 1.  
#' See details in \code{\link{par}}.
#' @param CIlwd The line width for confidence intervals, a positive number, 
#' defaulting to 1.
#' @param lty The line type. Line types can either be specified as an integer
#' (0 = blank, 1 = solid (default), 2 = dashed, 3 = dotted, 4 = dotdash, 
#' 5 = longdash, 6 = twodash).  See details in \code{\link{par}}.
#' @param CIlty The line type for confidence intervals. Line types can either 
#' be specified as an integer (0 = blank, 1 = solid (default), 2 = dashed,
#' 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash).
#' @param pcol  A specification for the points color.
#' @param cex A numerical value giving the amount by which plotting symbols
#' should be magnified relative to the default. See details in \code{\link{par}}.
#' @param \ldots Other options.
#' 
#'@return Simply produce a plot.
#' 
#'@author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#'@examples
#' library(npregfast)
#' data(barnacle)
#' 
#' # Nonparametric regression without interactions
#' fit <- frfast(DW ~ RC, data = barnacle, nboot = 100) 
#' plot(fit)
#' plot(fit, der = 0)
#' plot(fit, der = 0, points = FALSE)
#' plot(fit, der = 1, col = "red", CIcol = "blue")
#' 
#' # Nonparametric regression with interactions
#' fit2 <- frfast(DW ~ RC : F, data = barnacle, nboot = 100) 
#' plot(fit2)
#' plot(fit2, der = 0, fac = "lens")
#' plot(fit2, der = 1, col = "grey", CIcol = "red")
#' plot(fit2, der = c(0,1), fac = c("barca","lens"))
#' 
#' @importFrom graphics lines par plot
#' @export




plot.frfast <- function(x = model, y, fac = NULL, der = NULL, points = TRUE, 
                        xlab = model$name[2], ylab = model$name[1], ylim = NULL,
                        main = NULL, col = "black", CIcol = "black", 
                        ablinecol = "red", abline = TRUE, type = "l", 
                        CItype = "l", lwd = 2, CIlwd = 1, lty = 1, 
                        CIlty = 2, pcol = "grey80", cex = 0.6, ...) {
  # CIcol = 'grey50'

   
  model <- x
  
  nf <- model$nf
  fi <- length(fac)
  co <- length(der)
  facini <- fac
 
  if ((nf == 1) & (fi >= 1)) 
    stop("Argument \"fac\" not suported. 
         There is not factor in the model.")

 if(!is.null(fac) &   !identical(rep(TRUE,fi), fac %in% model$label)) {
   stop("\"",paste(fac[!fac %in% model$label]),"\" is not a factor's level.")
 }
  
  if (sum(der > 2) >= 1) 
    stop("",paste(der[which(der > 2)])," is not a r-th derivative implemented, only 
         permitted 0, 1 or 2.")

  if (missing(model)) 
    stop("Argument \"model\" is missing, with no default. 
         Must be a frfast object.")
  
  

  
#  if (fi == 0 & nf > 1 | fi > 1 | co == 0 & nf > 1 | co > 1 | co == 0 & nf == 1) {
    if (fi == 0) 
      fi <- nf
    if (co == 0) 
      co <- 3
    jnf <- c()
    par(mfrow = c(fi, co))
    if (is.null(fac)) {
     jnf <- c(1:nf)
      fac <- model$label
    } else {
      for (i in 1:length(fac)) {
        jnf[i] <- which(model$label == fac[i])
      }
    }
    
    if (is.null(der)) {
      jder <- c(1:3)
    } else {
      jder <- der + 1
    }
    
    for (j in jnf) {
      for (i in jder) {
        if (i == 1) {
          ylab2 <- ylab
          ylim2 <- c(min(model$ydata[model$fmod == model$numlabel[j]], na.rm = T), 
                     max(model$ydata[model$fmod == model$numlabel[j]], na.rm = T))
        } else {
          ylim2 <- c(min(model$p[, der = i, fac = j], na.rm = T), 
                     max(model$p[, der = i, fac = j], na.rm = T))
        }
        if (i == 2) 
          ylab2 <- "First derivative"
        if (i == 3) 
          ylab2 <- "Second derivative"
        
        if (is.null(main)){
          if(i == jder[1] & model$nf > 1){
            title <- paste("Level", model$label[j])
          }else{
            title <- NULL
          }
        }else{
         # if(length(main) == 1) main <- rep()
          if(i == jder[1]){
            title <- main[j]
          }else{
            title <- NULL
          }
        }
        
        if (is.null(ylim)) ylim <- ylim2  #### ver esto!!!!!
        plot(model$x, model$p[, der = i, fac = j], type = type, xlab = xlab, 
             ylab = ylab2, col = col, main = title, ylim = ylim, lwd = lwd, 
             lty = lty, ...)
        if ((points == TRUE) & (i == 1)) {
          points(model$xdata[model$fmod == model$numlabel[j]], model$ydata[model$fmod == model$numlabel[j]], 
                 col = pcol, cex = cex, ...)
          lines(model$x, model$p[, der = i, fac = j], type = type, xlab = xlab, 
                ylab = ylab2, col = col, main = title, ylim = ylim2, lwd = lwd, 
                lty = lty, ...)
        }
        lines(model$x, model$pl[, der = i, fac = j], lty = CIlty, col = CIcol, 
              type = CItype, lwd = CIlwd, ...)
        lines(model$x, model$pu[, der = i, fac = j], lty = CIlty, col = CIcol, 
              type = CItype, lwd = CIlwd, ...)
        ylim <- NULL # hay que poner el ylim nulo para el siguiente plot
        if (i == 3) {
          if (abline == TRUE) 
            abline(h = 0, col = ablinecol)
        }
      }
    }
    
    
    
    
    
#   } else {
#     # der=co+1 # esto estaba asi, pero no funcionaba plot(fit,der=1) #cuando
#     # formula=DW~RC
#     if (length(der) == 0) {
#       jder <- c(1:3)
#     } else {
#       der <- der + 1
#     }
#     jnf <- c()
#     title <- main
#     if (length(fac) != 0) {
#       for (i in 1:length(fac)) {
#         jnf[i] <- which(model$label == fac[i])
#       }
#     } else {
#       jnf <- 1
#     }
#     
#     if (length(fac) == 0 & is.null(main)) {
#       title <- ""
#     }
#     if (length(fac) != 0 & is.null(main)) 
#       title <- paste("Level", model$label[jnf])
#     
#     if (der == 1) 
#       ylab2 <- ylab
#     if (der == 2) 
#       ylab2 <- "First derivative"
#     if (der == 3) 
#       ylab2 <- "Second derivative"
#     if ((points == TRUE) & (der == 1)) {
#       ylim2 <- c(min(model$ydata[model$fmod == jnf], na.rm = T), 
#                  max(model$ydata[model$fmod == jnf], na.rm = T))
#     } else {
#       ylim2 <- c(min(model$p[, der, jnf], na.rm = T), 
#                  max(model$p[, der, jnf], na.rm = T))
#     }
#     if (is.null(ylim)) 
#       ylim <- ylim2
#     plot(model$x, model$p[, der, jnf], type = type, xlab = xlab, ylab = ylab2, 
#          col = col, main = title, ylim = ylim, lwd = lwd, lty = lty, ...)
#     if ((points == TRUE) & (der == 1)) {
#       points(model$xdata[model$fmod == jnf], model$ydata[model$fmod == jnf], 
#              col = "grey80", cex = 0.6, ...)
#       lines(model$x, model$p[, der, jnf], type = type, xlab = xlab, ylab = ylab2, 
#             col = col, main = title, lwd = lwd, lty = lty, ...)
#     }
#     lines(model$x, model$pl[, der, jnf], lty = CIlty, col = CIcol, type = CItype, 
#           lwd = CIlwd, ...)
#     lines(model$x, model$pu[, der, jnf], lty = CIlty, col = CIcol, type = CItype, 
#           lwd = CIlwd, ...)
#     if (der == 3) {
#       if (abline == TRUE) 
#         abline(h = 0, col = ablinecol)
#     }
#   }
} 