#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot


#' Visualization of \code{frfast} objects with ggplot2 graphics
#' @description Useful for drawing the estimated regression function, 
#' first and second derivative (for each factor's level) using ggplot2 graphics.
#'  Additionally, with the 
#' \code{diffwith} argument it is possible to draw the differences between
#'  two factor's levels.
#' @param object \code{frfast} object.
#' @param fac Factor's level to be taken into account
#' in the plot. By default is \code{NULL}.
#' @param der Number which determines any inference process. 
#' By default \code{der} is \code{0}. If this term is \code{0}, the plot 
#' shows the initial estimate. If it is \code{1} or \code{2},
#' it is designed for the first or second derivative, respectively.
#' @param diffwith Factor's level used for drawing the differences respect to the 
#' level specified in the \code{fac} argument.  By default, \code{NULL}. 
#' The differences are computed for the r-th derivative specified 
#' in the \code{der} argument.
#' @param points Draw the original data into the plot. By default it is
#' \code{TRUE}.
#' @param xlab A title for the \code{x} axis. 
#' @param ylab A title for the \code{y} axis. 
#' @param ylim The \code{y} limits of the plot.
#' @param main An overall title for the plot.
#' @param col A specification for the default plotting color.
#' @param CIcol A specification for the default confidence intervals
#' plotting color (for the fill).
#' @param CIlinecol A specification for the default confidence intervals
#' plotting color (for the edge).
#' @param pcol  A specification for the points color.
#' @param abline Draw an horizontal line into the plot of the second derivative 
#' of the model.
#' @param ablinecol The color to be used for \code{abline}.
#' @param lty The line type. Line types can either be specified as an integer
#' (0 = blank, 1 = solid (default), 2 = dashed, 3 = dotted, 4 = dotdash, 
#' 5 = longdash, 6 = twodash).  See details in \code{\link{par}}.
#' @param CIlty The line type for confidence intervals. Line types can either 
#' be specified as an integer (0 = blank, 1 = solid (default), 2 = dashed,
#' 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash).
#' @param lwd The line width, a positive number, defaulting to 1.  
#' See details in \code{\link{par}}.
#' @param CIlwd The line width for confidence intervals, a positive number, 
#' defaulting to 1.
#' @param cex A numerical value giving the amount by which plotting symbols
#' should be magnified relative to the default. See details in \code{\link{par}}.
#' @param alpha Alpha transparency for overlapping elements expressed 
#' as a fraction between 0 (complete transparency) and 1 (complete opacity).
#' @param \ldots Other options.
#' 
#'@return A ggplot object, so you can use common features from 
#' ggplot2 package to manipulate the plot.
#' 
#'@author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#'@examples
#'
#' library(npregfast)
#' library(ggplot2)
#' 
#' 
#' data(barnacle)
#' 
#' # Nonparametric regression without interactions
#' fit <- frfast(DW ~ RC, data = barnacle, nboot = 100) 
#' autoplot(fit)
#' autoplot(fit, points = FALSE) + ggtitle("Title")
#' autoplot(fit, der = 1) + xlim(4, 20)
#' #autoplot(fit, der = 1, col = "red", CIcol = "blue")
#' 
#' # Nonparametric regression with interactions
#' fit2 <- frfast(DW ~ RC : F, data = barnacle, nboot = 100) 
#' autoplot(fit2, fac = "barca")
#' plot(fit2, der = 1, fac = "lens")
#' 
#' # Visualization of the differences between two factor's levels
#' autoplot(fit2, fac = "barca", diffwith = "lens")
#' autoplot(fit2, der = 1, fac = "barca", diffwith = "lens")
#' 
#' 
#' #Plotting in the same graphics device
#' 
#' if (requireNamespace("gridExtra", quietly = TRUE)) {
#' 
#' # For plotting two derivatives in the same graphic windows 
#' ders <- lapply(0:1, function(x) autoplot(fit, der = x))
#' gridExtra::grid.arrange(grobs = ders, ncol = 2, nrow = 1)
#' 
#' # For plotting two levels in the same graphic windows 
#' facs <- lapply(c("barca", "lens"), function(x) autoplot(fit2, der = 0, fac = x))
#' gridExtra::grid.arrange(grobs = facs, ncol = 2, nrow = 1)
#' 
#' }
#' 
#' @importFrom ggplot2 autoplot geom_point aes_string ylab
#' @importFrom ggplot2 geom_hline geom_ribbon geom_line 
#' @importFrom ggplot2 ggtitle ggplot coord_cartesian
#' @export





autoplot.frfast <- function(object = model, fac = NULL, der = 0, diffwith = NULL, 
                            points = TRUE, xlab = model$name[2], ylab = model$name[1], 
                            ylim = NULL, main = NULL, col = "black", 
                            CIcol = "black", CIlinecol = "transparent", 
                            pcol = "grey80",  abline = TRUE, 
                            ablinecol = "red", lty = 1, CIlty = 2, lwd = 1, 
                            CIlwd = 1, cex = 1.4, alpha = 0.2, ...) {
  
  model <- object
  nf <- model$nf
  
  # Control
  
  if (missing(object)) 
    stop("Argument \"object\" is missing, with no default. 
         Must be a frfast object.")
  
  if ((nf == 1) & (!is.null(fac))) 
    stop("Argument \"fac\" not suported. 
         There is not factor in the model.") 
  
  
  if ((nf > 1) & (is.null(fac))) 
    stop("Argument \"fac\" is missing with no default. 
         Levels supported: ", paste(model$label, collapse = ", "),".")  
  
  
  if (!is.null(fac) &   !identical(rep(TRUE, 1), fac %in% model$label)) {
    stop("\"",paste(fac),"\" is not a factor's level.")
  }
  
  
  if (!(der %in% c(0,1,2))) 
    stop("",paste(der)," is not a r-th derivative implemented, only 
         permitted 0, 1 or 2.")
  
  
  
  
  if (is.null(diffwith)) {
    
    
    if (is.null(fac) & (nf == 1)) {
      jnf <- 1
      fac <- model$label
    } else {
      jnf <- which(model$label == fac)
    }
    
    
    jder <- der + 1
    
    
    if (jder == 1) {
      ylab2 <- ylab
      rgo <-  max(model$ydata[model$fmod == model$numlabel[jnf]], na.rm = TRUE) -
        min(model$ydata[model$fmod == model$numlabel[jnf]], na.rm = TRUE)
      ylim2 <- c(min(model$ydata[model$fmod == model$numlabel[jnf]], na.rm = TRUE) - (rgo*0.05), 
                 max(model$ydata[model$fmod == model$numlabel[jnf]], na.rm = TRUE) + (rgo*0.05))
    } else {
      rgo <-  max(model$p[, der = jder, fac = jnf], na.rm = TRUE) -
        min(model$p[, der = jder, fac = jnf], na.rm = TRUE)
      ylim2 <- c(min(model$p[, der = jder, fac = jnf], na.rm = TRUE) - (rgo*0.05), 
                 max(model$p[, der = jder, fac = jnf], na.rm = TRUE) + (rgo*0.05))
    }
    
    
    if (is.null(ylim)) ylim <- ylim2
    
    if (jder == 2)  ylab2 <- "First derivative"
    if (jder == 3)  ylab2 <- "Second derivative"
    
    
    if (is.null(main) & (nf > 1)) main <- paste("Level", model$label[jnf])
    
    
    #plot
    data_bin <- data.frame(x = model$x,
                           pl = model$pl[, der = jder, fac = jnf],
                           pu = model$pu[, der = jder, fac = jnf],
                           p = model$p[, der = jder, fac = jnf])
    
    data_ori <- data.frame(xdata = model$xdata[model$fmod == model$numlabel[jnf]],
                           ydata = model$ydata[model$fmod == model$numlabel[jnf]])
    
    
    if ((points == TRUE) & (jder == 1)) {
      points_layer <- ggplot2::geom_point(data = data_ori, 
                                          ggplot2::aes_string(x = "xdata",
                                                             y = "ydata"),
                                 colour = pcol, size = cex)
    }else{
      points_layer <- NULL
    }
    
    
    if ((jder == 3) & (abline == TRUE)) {
      abline_layer <- ggplot2::geom_hline(yintercept = 0, colour = ablinecol)
    }else{
      abline_layer <- NULL
    }
    
    
    ggplot2::ggplot() +
      points_layer +
      ggplot2::geom_ribbon(data = data_bin, ggplot2::aes_string(x = "x", 
                                              ymin = "pl", 
                                              ymax = "pu"), 
                  alpha = alpha, fill = CIcol, linetype = lty,
                  size = CIlwd, col = CIlinecol) +
      ggplot2::geom_line(data = data_bin, ggplot2::aes_string(x = "x", 
                                            y = "p"), 
                size = lwd, colour = col, linetype = lty, na.rm = TRUE) +
      abline_layer +
      ggplot2::coord_cartesian(ylim = ylim) +
      ggplot2::ylab(ylab2) +
      ggplot2::xlab(xlab) +
      ggplot2::ggtitle(main)
    
  }else{ # diffwith != NULL
    
    # Control
    if ((nf == 1) & (!is.null(diffwith))) 
      stop("Argument \"diffwith\" not suported. 
           There is not factor in the model.") 
    
    if (fac == diffwith) 
      stop("Argument 'fac' and 'diffwith' are not different")
    
    if (!isTRUE(diffwith %in% model$label)) {
      stop("\"",paste(diffwith),"\" is not a factor's level. Levels supported: ", 
           paste(model$label, collapse = ", "),".")
    }
    
    
    nf <- model$nf
    jnf <- c()
    jnf[1] <- which(model$label == fac)  
    jnf[2] <- which(model$label == diffwith) 
    
    
    jder <- der + 1
    
    
    if (jder == 1) ylab2 <- ylab
    if (jder == 2) ylab2 <- "First derivative"
    if (jder == 3) ylab2 <- "Second derivative"
    
    if (sum(model$diff[, der = jder, jnf[2], jnf[1]], na.rm = T) == 0) { # para ver si -1* o no
      
      if (is.null(ylim)) {
        rgo <- max(-1 * (model$diff[, der = jder, jnf[1], jnf[2]]), na.rm = T) -
          min(-1 * (model$diff[, der = jder, jnf[1], jnf[2]]), na.rm = T)
        
        ylim <- c(min(-1 * (model$diff[, der = jder, jnf[1], jnf[2]]), na.rm = T) - 
                    (rgo * 0.05),
                  max(-1 * (model$diff[, der = jder, jnf[1], jnf[2]]), na.rm = T) + 
                    (rgo * 0.05))
      }
      
      
      
      if (is.null(main)) main <- "Differences"
      
      
      data_bin <- data.frame(x = model$x,
                             p = -1 * model$diff[, der = jder, jnf[1], jnf[2]],
                             pl = -1 * model$diffl[, der = jder, jnf[1], jnf[2]],
                             pu = -1 * model$diffu[, der = jder, jnf[1], jnf[2]])
      
      if (abline == TRUE){
        abline_layer <- ggplot2::geom_hline(yintercept = 0, colour = ablinecol)
      }else{
        abline_layer <- NULL
      }
      
      
      ggplot2::ggplot() +
        ggplot2::geom_ribbon(data = data_bin, ggplot2::aes_string(x = "x", 
                                                ymin = "pl", 
                                                ymax = "pu"), 
                    alpha = alpha, fill = CIcol, linetype = lty,
                    size = CIlwd, col = CIlinecol) +
        ggplot2::geom_line(data = data_bin, ggplot2::aes_string(x = "x", 
                                              y = "p"), 
                  size = lwd, colour = col, linetype = lty, na.rm = TRUE) +
        abline_layer +
        ggplot2::coord_cartesian(ylim = ylim) +
        ggplot2::ylab(ylab2) +
        ggplot2::xlab(xlab) +
        ggplot2::ggtitle(main)
      
      
    } else {
      
      if (is.null(ylim)) {
        rgo <- max(1 * (model$diff[, der = jder, jnf[2], jnf[1]]), na.rm = T) -
          min(1 * (model$diff[, der = jder, jnf[2], jnf[1]]), na.rm = T)
        ylim <- c(min(1 * (model$diff[, der = jder, jnf[2], jnf[1]]), na.rm = T) - 
                    (rgo * 0.05),
                  max(1 * (model$diff[, der = jder, jnf[2], jnf[1]]), na.rm = T) +
                    (rgo * 0.05))
      }
      
      
      if (is.null(main)) main <- "Differences"
      
      data_bin <- data.frame(x = model$x,
                             p = model$diff[, der = jder, jnf[2], jnf[1]],
                             pl = model$diffl[, der = jder, jnf[2], jnf[1]],
                             pu = model$diffu[, der = jder, jnf[2], jnf[1]])
      
      if (abline == TRUE) {
        abline_layer <- ggplot2::geom_hline(yintercept = 0, colour = ablinecol)
      }else{
        abline_layer <- NULL
      }
      
      
      ggplot2::ggplot() +
        ggplot2::geom_ribbon(data = data_bin, ggplot2::aes_string(x = "x", 
                                                ymin = "pl", 
                                                ymax = "pu"), 
                    alpha = alpha, fill = CIcol, linetype = lty,
                    size = CIlwd, col = CIlinecol) +
        ggplot2::geom_line(data = data_bin, ggplot2::aes_string(x = "x", 
                                              y = "p"), 
                  size = lwd, colour = col, linetype = lty, na.rm = TRUE) +
        abline_layer +
        ggplot2::coord_cartesian(ylim = ylim) +
        ggplot2::ylab(ylab2) +
        ggplot2::xlab(xlab) +
        ggplot2::ggtitle(main)
      
    }
  }
} 