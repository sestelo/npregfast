plot.diff <-
function(model,factor2,factor1,der=NULL,est.include=FALSE,xlab=model$name[2],ylab=model$name[1],ylim=NULL, main=NULL,col="black",CIcol="grey50",ablinecol="red",abline=TRUE,type="l",CItype="l",lwd=1, CIlwd=1.5,lty=1, CIlty=2,...){
	nf=model$nf
	#co=length(der)
	jnf=c()
	jnf[1]=which(model$label==factor1)#"B" plot.diff(ajus,"A","B");plot.diff(ajus,1,2) 
	jnf[2]=which(model$label==factor2) #"A"
	#if(length(der)==0) {jder=c(1:3)}else{jder=der+1}	
	
## Argumentos control
if(missing(model)) stop("Argument \"model\" is missing, with no default. Must be a frfast object.")
if(sum(der >2)>=1) stop("Argument \"der\" not suported.")
if(missing(factor1)&missing(factor2)) stop("Argument 'factor1' and/or 'factor2' are not set or have length zero")
if(factor1==factor2) stop("Argument 'factor1' and 'factor2' are not different")
	
	
if(est.include==FALSE){
	if(is.null(der)) der=c(0,1,2) 
	der=der+1
	par(mfrow=c(1,length(der)))
	for (i in der){		
		if(i==1) ylab2=ylab 
		if(i==2) ylab2="First derivative"
		if(i==3) ylab2="Second derivative"
		if(sum(model$diff[,der=i,jnf[2],jnf[1]],na.rm=T)==0){
			if(is.null(ylim)){ylim=c(min(-1*(model$diff[,der=i,jnf[1],jnf[2]]),na.rm=T),max(-1*(model$diff[,der=i,jnf[1],jnf[2]]),na.rm=T))}
			if(is.null(main)){ title="Differences"}else{title=main}
			plot( model$x, -1*( model$diff[,der=i, jnf[1],jnf[2]]), type=type, xlab=xlab, ylab=ylab2, col=col,main=title,ylim= ylim,lty = lty,lwd = lwd,...)
			lines(model$x,-1*(model$diffl[,der=i,jnf[1],jnf[2]]),lty=CIlty,col=CIcol,type=CItype,lwd=CIlwd,...)
			lines(model$x,-1*(model$diffu[,der=i,jnf[1],jnf[2]]),lty=CIlty,col=CIcol,type=CItype,lwd=CIlwd,...)
			if(abline==TRUE)abline(h=0,col=ablinecol)		
		}else{
			if(is.null(ylim)){ylim=c(min(-1*(model$diff[,der=i,jnf[2],jnf[1]]),na.rm=T),max(-1*(model$diff[,der=i,jnf[2],jnf[1]]),na.rm=T))}
			if(is.null(main)){ title="Differences"}else{title=main}
			plot(model$x,model$diff[,der=i,jnf[2],jnf[1]],ylim=ylim, 
				type=type,ylab=ylab2, xlab=xlab, main=title,lty = lty,...)
			lines(model$x,model$diffl[,der=i,jnf[2],jnf[1]],lty=CIlty,col=CIcol,type=CItype,lwd=lwd,...)
			lines(model$x,model$diffu[,der=i,jnf[2],jnf[1]],lty=CIlty,col=CIcol,type=CItype,lwd=lwd,...)
			if(abline==TRUE)abline(h=0,col=ablinecol)	}}
}else{
	if(is.null(der)) der=c(0,1,2) 
	jder=der+1  #if(length(der)==0) {jder=c(1:3)}else{jder=der+1}
	par(mfrow=c(nf+1,length(der)))
	for (i in jder){		
		if(i==1) ylab2=ylab 
		if(i==2) ylab2="First derivative"
		if(i==3) ylab2="Second derivative"	
		if(sum(model$diff[,der=i,jnf[2],jnf[1]],na.rm=T)==0){
			if(is.null(ylim)){ylim=c(min(-1*(model$diff[,der=i,jnf[1],jnf[2]]),na.rm=T),max(-1*(model$diff[,der=i,jnf[1],jnf[2]]),na.rm=T))}
			
			if(is.null(main)){ title="Differences"}else{title=main}
			plot(model$x,-1*(model$diff[,der=i,jnf[1],jnf[2]]),type=type,xlab=xlab,ylab=ylab2,col=col,main=title,ylim=ylim,lty = lty,lwd=lwd,...)
			lines(model$x,-1*(model$diffl[,der=i,jnf[1],jnf[2]]),lty=CIlty,col=CIcol,type=CItype,lwd=CIlwd,...)
			lines(model$x,-1*(model$diffu[,der=i,jnf[1],jnf[2]]),lty=CIlty,col=CIcol,type=CItype,lwd=CIlwd,...)
			if(abline==TRUE)abline(h=0,col=ablinecol)	
	
			for(j in length(jnf):1){
				if( (is.null(main)) &	(jnf[j]==jnf[2]) ){title=paste("Level",model$label[jnf[2]])		
				}else if((is.null(main))&	(jnf[j]==jnf[1])){title=paste("Level",model$label[jnf[1]])
				}else{title=main}
				plot(model$x,model$p[,der=i,jnf[j]],type=type,xlab=xlab,ylab=ylab2,col=col,main=title,ylim=c(min(model$p[,der=i,jnf[j]],na.rm=T),max(model$p[,der=i,jnf[j]],na.rm=T)),lty = lty,lwd=lwd,...)
				lines(model$x,model$pl[,der=i,jnf[j]],lty=CIlty,col=CIcol,type=CItype,lwd=CIlwd,...)
				lines(model$x,model$pu[,der=i,jnf[j]],lty=CIlty,col=CIcol,type=CItype,lwd=CIlwd,...)
				if(i==3){if(abline==TRUE)abline(h=0,col=ablinecol)}
			}
	 
		}else{
			if(is.null(main)){ title="Differences"}else{title=main}
			if(is.null(ylim)){ylim=c(min(-1*(model$diff[,der=i,jnf[2],jnf[1]]),na.rm=T),max(-1*(model$diff[,der=i,jnf[2],jnf[1]]),na.rm=T))}
			plot(model$x,model$diff[,der=i,jnf[2],jnf[1]],ylim=ylim,
				type=type,ylab=ylab2,xlab=xlab,main=title,lty = lty, lwd=lwd,lwd=CIlwd,...)
			lines(model$x,model$diffl[,der=i,jnf[2],jnf[1]],lty=CIlty,col=CIcol,type=CItype,lwd=CIlwd,...)
			lines(model$x,model$diffu[,der=i,jnf[2],jnf[1]],lty=CIlty,col=CIcol,type=CItype,lwd=CIlwd,...)
			if(abline==TRUE)abline(h=0,col=ablinecol)	

			for(j in length(jnf):1){
				if( (is.null(main)) &	(jnf[j]==jnf[2]) ){title=paste("Level",model$label[jnf[2]])		
				}else if((is.null(main))&	(jnf[j]==jnf[1])){title=paste("Level",model$label[jnf[1]])
				}else{title=main}
					plot(model$x,model$p[,der=i,jnf[j]],type=type,xlab=xlab,ylab=ylab2,col=col,main=title, ylim=c(min(model$p[,der=i,jnf[j]],na.rm=T),max(model$p[,der=i,jnf[j]],na.rm=T)),lty = lty,lwd=lwd, ...)
					lines(model$x,model$pl[,der=i,jnf[j]],lty=CIlty,col=CIcol,type=CItype,lwd=CIlwd,...)
					lines(model$x,model$pu[,der=i,jnf[j]],lty=CIlty,col=CIcol,type=CItype,lwd=CIlwd,...)
					if(i==3){if(abline==TRUE)abline(h=0,col=ablinecol)}
			}
		}
	}
}}