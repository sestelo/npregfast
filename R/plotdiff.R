#' Visualization of the differences between the estimated curves 
#' for two factor's levels
#' @description Useful for drawing the differences between the estimation of 
#' curves (initial estimate, first or second derivative) for  two factor's levels.
#' Missing values of factor's levels is not allowed. 
#'@param model Parametric or nonparametric regression out 
#' obtained by \code{\link{frfast}} function.
#'@param level2 Second factor's level at which to perform the 
#'differences between curves.
#'@param level1 First factor's level at which to perform the 
#'differences between curves.
#' @param der Number or vector which determines any inference process. 
#' By default \code{der} is \code{NULL}. If this term is \code{0}, the plot 
#' shows the differences between estimated regression functions. If it is 
#' \code{1} or \code{2}, it is designed for the first or second derivative, 
#' respectively.
#'@param est.include Draws the estimates of the model. 
#'By default it is \code{FALSE}.
#' @param xlab A title for the x axis. 
#' @param ylab A title for the y axis. 
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
#' @param \ldots Other options.
#' @return Simply produce a plot.
#' @author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#' @examples
#' library(npregfast)
#' data(barnacle)
#' 
#' # Nonparametric regression with interactions
#' fit2 <- frfast(DW ~ RC : F, data = barnacle, nboot = 100) 
#' plotdiff(fit2, level2 = "lens", level1 = "barca")
#' plotdiff(fit2, level2 = "lens", level1 = "barca", der = 1, col = "blue", CIcol = "grey")
#' plotdiff(fit2, "lens", "barca", der = c(0, 1), ylim = c(-0.05, 0.05))
#' 
#' 
#' 
#' @importFrom graphics lines par plot
#' @export


plotdiff <- function(model, level2, level1, der = NULL, est.include = FALSE, 
                      xlab = model$name[2], ylab = model$name[1], ylim = NULL, 
                      main = NULL, col = "black", CIcol = "grey50", ablinecol = "red", 
                      abline = TRUE, type = "l", CItype = "l", lwd = 1, CIlwd = 1.5, 
                      lty = 1, CIlty = 2, ...) {
  nf <- model$nf
  # co=length(der)
  jnf <- c()
  jnf[1] <- which(model$label == level1)  #'B' plot.diff(ajus,'A','B');plot.diff(ajus,1,2) 
  jnf[2] <- which(model$label == level2)  #'A'
  # if(length(der)==0) {jder=c(1:3)}else{jder=der+1}
  
  ## Argumentos control
  if (missing(model)) 
    stop("Argument \"model\" is missing, with no default. 
         Must be a frfast object.")
  
  if (missing(level1) & missing(level2)) 
    stop("Argument 'level1' and/or 'level2' are missing, with no default")
  
  if (level1 == level2) 
    stop("Argument 'level1' and 'level2' are not different")
  
  if(!isTRUE(level1 %in% model$label)) {
    stop("\"",paste(level1),"\" is not a factor's level.")
  }
  
  if(!isTRUE(level2 %in% model$label)) {
    stop("\"",paste(level2),"\" is not a factor's level.")
  }
  
  if (sum(der > 2) >= 1) 
    stop("",paste(der[which(der > 2)])," is not a r-th derivative implemented, only 
         permitted 0, 1 or 2.")
  
  if ((nf == 1)) 
    stop("Function \"plot.diff\" not suported. 
         There is not factor in the model.")
  
  
  
  
  if (est.include == FALSE) {
    if (is.null(der)) 
      der <- c(0, 1, 2)
    der <- der + 1
    par(mfrow = c(1, length(der)))
    for (i in der) {
      if (i == 1) 
        ylab2 <- ylab
      if (i == 2) 
        ylab2 <- "First derivative"
      if (i == 3) 
        ylab2 <- "Second derivative"
      if (sum(model$diff[, der = i, jnf[2], jnf[1]], na.rm = T) == 0) { # para ver si -1* o no
        
        if (is.null(ylim)) {
          ylim <- c(min(-1 * (model$diff[, der = i, jnf[1], jnf[2]]), na.rm = T),
                    max(-1 * (model$diff[, der = i, jnf[1], jnf[2]]), na.rm = T))
        }
        if (is.null(main)){
          if(i == der[1] ){
            title <- "Differences"
          }else{
            title <- NULL
          }
        }else{
          if(i == der[1]){
            title <- main
          }else{
            title <- NULL
          }
        }
        plot(model$x, -1 * (model$diff[, der = i, jnf[1], jnf[2]]), 
             type = type, xlab = xlab, ylab = ylab2, col = col, main = title, 
             ylim = ylim, lty = lty, lwd = lwd, ...)
        lines(model$x, -1 * (model$diffl[, der = i, jnf[1], jnf[2]]), 
              lty = CIlty, col = CIcol, type = CItype, lwd = CIlwd, ...)
        lines(model$x, -1 * (model$diffu[, der = i, jnf[1], jnf[2]]), 
              lty = CIlty, col = CIcol, type = CItype, lwd = CIlwd, ...)
        if (abline == TRUE) 
          abline(h = 0, col = ablinecol)
        
        
      } else {
        
        if (is.null(ylim)) {
          ylim <- c(min(1 * (model$diff[, der = i, jnf[2], jnf[1]]), na.rm = T),
                    max(1 * (model$diff[, der = i, jnf[2], jnf[1]]), na.rm = T))
        }
        if (is.null(main)){
          if(i == der[1] ){
            title <- "Differences"
          }else{
            title <- NULL
          }
        }else{
          if(i == der[1]){
            title <- main
          }else{
            title <- NULL
          }
        }
        plot(model$x, model$diff[, der = i, jnf[2], jnf[1]], 
             ylim = ylim, type = type, ylab = ylab2, xlab = xlab, 
             main = title, lty = lty, ...)
        lines(model$x, model$diffl[, der = i, jnf[2], nf[1]], lty = CIlty, 
              col = CIcol, type = CItype, lwd = lwd, ...)
        lines(model$x, model$diffu[, der = i, jnf[2], jnf[1]], lty = CIlty, 
              col = CIcol, type = CItype, lwd = lwd, ...)
        if (abline == TRUE) 
          abline(h = 0, col = ablinecol)
      }
    }
    
    
    
  } else {  # est.include = TRUE
    
    if (is.null(der)) 
      der <- c(0, 1, 2)
    jder <- der + 1  #if(length(der)==0) {jder=c(1:3)}else{jder=der+1}
    par(mfrow = c(nf + 1, length(der)))
    for (i in jder) {
      if (i == 1) 
        ylab2 <- ylab
      if (i == 2) 
        ylab2 <- "First derivative"
      if (i == 3) 
        ylab2 <- "Second derivative"
      if (sum(model$diff[, der = i, jnf[2], jnf[1]], na.rm = T) == 0) {
        if (is.null(ylim)) {
          ylim <- c(min(-1 * (model$diff[, der = i, jnf[1], jnf[2]]), na.rm = T),
                    max(-1 * (model$diff[, der = i, jnf[1], jnf[2]]), na.rm = T))
        }

        
        
        if (is.null(main)){
          if(i == jder[1] ){
            title <- "Differences"
          }else{
            title <- NULL
          }
        }else{
          if(i == jder[1]){
            title <- main
          }else{
            title <- NULL
          }
        }
        

        plot(model$x, -1 * (model$diff[, der = i, jnf[1],  jnf[2]]), type = type, 
             xlab = xlab, ylab = ylab2, col = col, main = title, ylim = ylim, 
             lty = lty, lwd = lwd, ...)
        lines(model$x, -1 * (model$diffl[, der = i, jnf[1], jnf[2]]), 
              lty = CIlty, col = CIcol, type = CItype, lwd = CIlwd, ...)
        lines(model$x, -1 * (model$diffu[, der = i, jnf[1], jnf[2]]), 
              lty = CIlty, col = CIcol, type = CItype, lwd = CIlwd, ...)
        if (abline == TRUE) 
          abline(h = 0, col = ablinecol)
        
        for (j in length(jnf):1) {
          if (jnf[j] == jnf[2]) {
            title <- paste("Level", model$label[jnf[2]])
          } 
          if  (jnf[j] == jnf[1]) {
            title <- paste("Level", model$label[jnf[1]])
          }
          
          plot(model$x, model$p[, der = i, jnf[j]], type = type, 
               xlab = xlab, ylab = ylab2, col = col, main = title, 
               ylim = c(min(model$p[, der = i, jnf[j]], na.rm = T), 
                        max(model$p[, der = i, jnf[j]], na.rm = T)), 
               lty = lty, lwd = lwd, ...)
          lines(model$x, model$pl[, der = i, jnf[j]], 
                lty = CIlty, col = CIcol, type = CItype, 
                lwd = CIlwd, ...)
          lines(model$x, model$pu[, der = i, jnf[j]], 
                lty = CIlty, col = CIcol, type = CItype, 
                lwd = CIlwd, ...)
          if (i == 3) {
            if (abline == TRUE) 
              abline(h = 0, col = ablinecol)
          }
        }
        
      } else {
        
        
        if (is.null(main)){
          if(i == jder[1]){
            title <- "Differences"
          }else{
            title <- NULL
          }
        }else{
          if(i == jder[1]){
            title <- main
          }else{
            title <- NULL
          }
        }
        
        
        
        if (is.null(ylim)) {
          ylim <- c(min(1 * (model$diff[, der = i, jnf[2], jnf[1]]), na.rm = T),
                    max(1 * (model$diff[, der = i, jnf[2], jnf[1]]), na.rm = T))
        }
        plot(model$x, model$diff[, der = i, jnf[2], jnf[1]], 
             ylim = ylim, type = type, ylab = ylab2, xlab = xlab, 
             main = title, lty = lty, lwd = lwd, ...)
        lines(model$x, model$diffl[, der = i, jnf[2], jnf[1]], lty = CIlty, 
              col = CIcol, type = CItype, lwd = CIlwd, ...)
        lines(model$x, model$diffu[, der = i, jnf[2],  jnf[1]], lty = CIlty, 
              col = CIcol, type = CItype, lwd = CIlwd, ...)
        if (abline == TRUE) 
          abline(h = 0, col = ablinecol)
        
        for (j in length(jnf):1) {
          if (jnf[j] == jnf[2]) {
            title <- paste("Level", model$label[jnf[2]])
          } 
          if (jnf[j] == jnf[1]) {
            title <- paste("Level", model$label[jnf[1]])
          } 
         
          
          plot(model$x, model$p[, der = i, jnf[j]], type = type, 
               xlab = xlab, ylab = ylab2, col = col, main = title, 
               ylim = c(min(model$p[, der = i, jnf[j]], 
                            na.rm = T), max(model$p[, der = i, jnf[j]], 
                                            na.rm = T)), lty = lty, lwd = lwd, ...)
          lines(model$x, model$pl[, der = i, jnf[j]], 
                lty = CIlty, col = CIcol, type = CItype, 
                lwd = CIlwd, ...)
          lines(model$x, model$pu[, der = i, jnf[j]], 
                lty = CIlty, col = CIcol, type = CItype, 
                lwd = CIlwd, ...)
          if (i == 3) {
            if (abline == TRUE) 
              abline(h = 0, col = ablinecol)
          }
        }
      }
    }
  }
} 