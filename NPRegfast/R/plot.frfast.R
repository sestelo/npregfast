#' @export plot.frfast

plot.frfast<-function (model, fac = NULL, der = NULL, points = TRUE, xlab = model$name[2], ylab = model$name[1],ylim = NULL, main=NULL, col = "black",CIcol = "black",ablinecol="red", abline=TRUE,type = "l", CItype = "l", lwd = 2, CIlwd=1,lty=1, CIlty=2,...){ #CIcol = "grey50"		
	nf = model$nf
    fi = length(fac)
    co = length(der)
    facini=fac
## Argumentos control	
if((nf==1)&(fi>=1)) stop("Argument \"fac\" not suported. There is not factor in the model.") 
if(nf<fi) stop("The specified factor is not correct.")
if(sum(der >2)>=1) stop("Argument \"der\" not suported.")
if(missing(model)) stop("Argument \"model\" is missing, with no default. Must be a frfast object.")

 
if(fi==0&nf>1|fi>1|co==0&nf>1|co>1|co==0&nf==1){
    if (fi == 0) fi = nf 
    if (co == 0) co = 3
    jnf = c()
    par(mfrow = c(fi, co))
    if (length(fac) == 0){jnf = c(1:nf);fac = model$label
    }else{for (i in 1:length(fac)) {jnf[i] = which(model$label == fac[i])}} 
    
    if (length(der) == 0) {jder = c(1:3)
    }else{jder = der + 1}
    
    for (j in jnf) {
        for (i in jder) {
            if (i == 1) {ylab2 = ylab;ylim2 = c(min(model$ydata[model$fmod == j], na.rm = T), 
                  max(model$ydata[model$fmod == j], na.rm = T))
            }else{ylim2=c(min(model$p[,der = i,fac = j],na.rm=T),max(model$p[,der = i,fac = j],na.rm = T))}
            if (i == 2) ylab2 = "First derivative"
            if (i == 3) ylab2 = "Second derivative"
            if(length(main)==0){title=main}else{title=main[j]}     
            if(length(fac)!=0&is.null(main) ){ 
            	if(fi==nf){title=""}else{title = paste("Level", model$label[jnf])}}
            if(length(facini)==0&is.null(main)){title=" "}
            if(is.null(ylim)) ylim=ylim2  #### ver esto!!!!!
            plot(model$x, model$p[, der = i, fac = j], type = type, 
                xlab = xlab, ylab = ylab2, col = col, main = title, ylim = ylim,lwd = lwd,lty = lty,...)
          		  if ((points == TRUE) & (i == 1)) {
                  points(model$xdata[model$fmod==j],model$ydata[model$fmod==j],col="grey80",cex = 0.6,...)
                  lines(model$x, model$p[, der = i, fac = j], type = type, 
                  xlab = xlab, ylab = ylab2, col = col, main = title,ylim = ylim2,lwd=lwd, lty=lty, ...)}
            lines(model$x, model$pl[, der = i, fac = j], lty = CIlty, col = CIcol, type = CItype,lwd=CIlwd,...)
            lines(model$x, model$pu[, der = i, fac = j], lty = CIlty, col = CIcol, type = CItype,lwd=CIlwd,...)
            ylim=NULL
            if (i == 3) {if(abline==TRUE) abline(h = 0, col = ablinecol)}}}
}else{	 
		#der=co+1  # esto estaba asi, pero no funcionaba plot(fit,der=1) #cuando formula=DW~RC
    	if (length(der) == 0) {jder = c(1:3)
  		}else{der = der + 1}
    	jnf = c()
    	title=main 
    	if(length(fac)!=0){ for (i in 1:length(fac)) {jnf[i] = which(model$label == fac[i])}}else{jnf=1}
    	
        if(length(fac)==0&is.null(main)){title=""}
        if(length(fac)!=0&is.null(main) ) title = paste("Level", model$label[jnf])  
    	
    	if(der==1) ylab2=ylab
    	if(der==2) ylab2="First derivative"
    	if(der==3) ylab2="Second derivative"
    	if((points == TRUE) & (der == 1)) {
    		ylim2=c(min(model$ydata[model$fmod==jnf],na.rm=T),max(model$ydata[model$fmod==jnf],na.rm=T))
    	}else{ylim2=c(min(model$p[,der,jnf],na.rm=T),max(model$p[,der,jnf],na.rm=T))}
           if(is.null(ylim)) ylim=ylim2
    	plot(model$x, model$p[, der, jnf], type = type, xlab = xlab, ylab = ylab2, col = col, main = title,
    	     ylim = ylim,lwd=lwd,lty=lty,...)
            if ((points == TRUE) & (der == 1)) {
                points(model$xdata[model$fmod==jnf],model$ydata[model$fmod==jnf], col = "grey80", cex = 0.6,...)
                lines(model$x,model$p[,der,jnf],type=type,xlab = xlab, ylab = ylab2, col = col, main = title,lwd=lwd,lty=lty,...)}
        lines(model$x, model$pl[, der , jnf], lty = CIlty, col = CIcol, type = CItype,lwd=CIlwd,...)
        lines(model$x, model$pu[, der , jnf], lty = CIlty, col = CIcol, type = CItype,lwd=CIlwd,...)
             if (der == 3) {if(abline==TRUE) abline(h = 0, col = ablinecol)}
}}