#' @export plot.frfast

plot.frfast<-function (model, fac = NULL, der = NULL, points = TRUE, xlab = "x", ylab = "y", col = "black",
                       ICcol = "grey50", main=NULL,type = "l",ylim=NULL, ICtype = "l", lwd = 2, ...){
	
    nf = model$nf
    fi = length(fac)
    co = length(der)
    facini=fac
#if(fi==0&nf>1|fi>1|co==0&nf>1|co>=0){
    if (fi == 0) fi = nf # for each of the factor's level
    if (co == 0) co = 1 #only show estimation
    jnf = c()
 
    if(fi>0|co>0){par(mfrow = c(fi, co))}
    if (length(fac) == 0){jnf = c(1:nf);fac = model$etiquetas
    }else{for (i in 1:length(fac)) {jnf[i] = which(model$etiquetas == fac[i])}}
    
    if (length(der) == 0) {jder = c(1)} #only show estimation
    else {jder = der + 1}
    
    for (j in jnf) {
        for (i in jder) {
            if (i == 1) {ylab2 = ylab;ylim2 = c(min(model$ydata[model$fmod == j], na.rm = T), 
                  max(model$ydata[model$fmod == j], na.rm = T))
            }else{ylim2=c(min(model$p[,der = i,fac = j],na.rm=T),max(model$p[,der = i,fac = j],na.rm = T))}
            if (i == 2) ylab2 = "first derivative"
            if (i == 3) ylab2 = "second derivative"
                
            title=main
             
            if(length(fac)!=1&is.null(main) ) title = paste("Level", model$etiquetas[j])
           # if(length(facini)==0&is.null(main)){title=" "}
           # if(is.null(ylim)) ylim=ylim2
            plot(model$x, model$p[, der = i, fac = j], type = type, 
                xlab = xlab, ylab = ylab2, col = col, main = title, ylim = ylim)
          		  if ((points == TRUE) & (i == 1)) {
                  points(model$xdata[model$fmod == j], model$ydata[model$fmod == j], col = "grey80", cex = 0.6)
                  lines(model$x, model$p[, der = i, fac = j], type = type, 
                  xlab = xlab, ylab = ylab2, col = col, main = title,ylim = ylim2)}
            lines(model$x, model$pl[, der = i, fac = j], lty = 2, col = ICcol, type = ICtype)
            lines(model$x, model$pu[, der = i, fac = j], lty = 2, col = ICcol, type = ICtype)
            if (i == 3) abline(h = 0, col = 2)}}
#}

}



# }else{	 
#     	der=der+1
#     	jnf = c()
#     	title=main 
#     	if(length(fac)!=0){ for (i in 1:length(fac)) {jnf[i] = which(model$etiquetas == fac[i])}}else{jnf=1}
#     	
#        if(length(fac)==0&is.null(main)){title=""}
#        if(length(fac)!=0&is.null(main) ) title = paste("Level", model$etiquetas[jnf])  
#     	
#     	if(der==1) ylab2=ylab
#     	if(der==2) ylab2="first derivative"
#     	if(der==3) ylab2="second derivative"
#     	if((points == TRUE) & (der == 1)) {
#     		ylim2=c(min(model$ydata[model$fmod==jnf],na.rm=T),max(model$ydata[model$fmod==jnf],na.rm=T))
#     	}else{ylim2=c(min(model$p[,der,jnf],na.rm=T),max(model$p[,der,jnf],na.rm=T))}
#            if(is.null(ylim)) ylim=ylim2
#     	plot(model$x, model$p[, der, jnf], type = type, xlab = xlab, ylab = ylab2, col = col, main = title,
#     	     ylim = ylim)
#             if ((points == TRUE) & (der == 1)) {
#                 points(model$xdata[model$fmod==jnf],model$ydata[model$fmod==jnf], col = "grey80", cex = 0.6)
#                 lines(model$x,model$p[,der,jnf],type=type,xlab = xlab, ylab = ylab2, col = col, main = title)}
#        lines(model$x, model$pl[, der , jnf], lty = 2, col = ICcol, type = ICtype)
#        lines(model$x, model$pu[, der , jnf], lty = 2, col = ICcol, type = ICtype)
#              if (der == 3) abline(h = 0, col = 2)
# }}