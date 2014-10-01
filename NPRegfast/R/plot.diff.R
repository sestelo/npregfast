plot.diff <-
function(model,factor2,factor1,der=NULL,est.include=TRUE,xlab="x",ylab="y",col="black",ICcol="grey50",type="l",ICtype="l",main="title",...){ 
	if(missing(factor1)&missing(factor2))
		stop("Argument 'factor1' and/or 'factor2' are not set or have length zero")
	if(factor1==factor2)
		stop("Argument 'factor1' and 'factor2' are not different")
	nf=model$nf
	co=length(der)
	if(co==0) co=3
	jnf=c()
	jnf[1]=which(model$etiquetas==factor1)#"B" plot.diff(ajus,"A","B");plot.diff(ajus,1,2) 
	jnf[2]=which(model$etiquetas==factor2) #"A"
	if(length(der)==0) {jder=c(1:3)}else{jder=der+1}
	if(est.include==FALSE){
	par(mfrow=c(1,co))
	for (i in jder){		
	if(i==1) ylab2=ylab 
	if(i==2) ylab2="first derivative"
	if(i==3) ylab2="second derivative"
	if(sum(model$diff[,der=i,jnf[2],jnf[1]],na.rm=T)==0){
	plot(model$x,-1*(model$diff[,der=i,jnf[1],jnf	[2]]),type=type,xlab=xlab,ylab=ylab2,col=col,main="Differences",ylim=c(min(-1*(model$diffu[,der=i,jnf[1],jnf	[2]]),na.rm=T),max(-1*(model$diffl[,der=i,jnf[1],jnf[2]]),na.rm=T)))
	lines(model$x,-1*(model$diffl[,der=i,jnf[1],jnf[2]]),lty=2,col=col,type=ICtype)
	lines(model$x,-1*(model$diffu[,der=i,jnf[1],jnf[2]]),lty=2,col=col,type=ICtype)
	abline(h=0,col=2)
	}else{plot(model$x,model$diff[,der=i,jnf[2],jnf	[1]],ylim=c(min(model$diffl[,der=i,jnf[2],jnf[1]],na.rm=T),max(model$diffu[,der=i,jnf[2],jnf[1]],na.rm=T)),type=type,ylab=ylab2,xlab=xlab,main="Differences")
	lines(model$x,model$diffl[,der=i,jnf[2],jnf[1]],lty=2,col=col,type=ICtype)
	lines(model$x,model$diffu[,der=i,jnf[2],jnf[1]],lty=2,col=col,type=ICtype)
	abline(h=0,col=2)}}}else{	
	
	par(mfrow=c(nf,co))
	#if(length(der)==0) {jder=c(1:3)}else{jder=der+1}

for (i in jder){		
	if(i==1) ylab2=ylab 
	if(i==2) ylab2="first derivative"
	if(i==3) ylab2="second derivative"
	
	if(sum(model$diff[,der=i,jnf[2],jnf[1]],na.rm=T)==0){
	plot(model$x,-1*(model$diff[,der=i,jnf[1],jnf	[2]]),type=type,xlab=xlab,ylab=ylab2,col=col,main="Differences",ylim=c(min(-1*(model$diffu[,der=i,jnf[1],jnf	[2]]),na.rm=T),max(-1*(model$diffl[,der=i,jnf[1],jnf[2]]),na.rm=T)))
	lines(model$x,-1*(model$diffl[,der=i,jnf[1],jnf[2]]),lty=2,col=col,type=ICtype)
	lines(model$x,-1*(model$diffu[,der=i,jnf[1],jnf[2]]),lty=2,col=col,type=ICtype)
	abline(h=0,col=2)
	
		for(j in length(jnf):1){
		if((main=="title")&	(jnf[j]==jnf[2])){title=paste("Level",model$etiquetas[jnf[2]])
		
		}else{title=paste("Level",model$etiquetas[jnf[1]])}
		plot(model$x,model$p[,der=i,jnf[j]],type=type,xlab=xlab,ylab=ylab2,col=col,main=title,ylim=c(0,max(model$pu[,der=i,jnf[j]],na.rm=T)))
		lines(model$x,model$pl[,der=i,jnf[j]],lty=2,col=ICcol,type=ICtype)
		lines(model$x,model$pu[,der=i,jnf[j]],lty=2,col=ICcol,type=ICtype)}
	 
	}else{
	plot(model$x,model$diff[,der=i,jnf[2],jnf	[1]],ylim=c(min(model$diffl[,der=i,jnf[2],jnf[1]],na.rm=T),max(model$diffu[,der=i,jnf[2],jnf[1]],na.rm=T)),type=type,ylab=ylab2,xlab=xlab,main="Differences")
	lines(model$x,model$diffl[,der=i,jnf[2],jnf[1]],lty=2,col=col,type=ICtype)
	lines(model$x,model$diffu[,der=i,jnf[2],jnf[1]],lty=2,col=col,type=ICtype)
	abline(h=0,col=2)

		for(j in length(jnf):1){
		if((main=="title")&	(jnf[j]==jnf[2])){title=paste("Level",model$etiquetas[jnf[2]])		}else{title=paste("Level",model$etiquetas[jnf[1]])}
		plot(model$x,model$p[,der=i,jnf		[j]],type=type,xlab=xlab,ylab=ylab2,col=col,main=title,ylim=c(0,max(model$pu		[,der=i,jnf[j]],na.rm=T)))
		lines(model$x,model$pl[,der=i,jnf[j]],lty=2,col=ICcol,type=ICtype)
		lines(model$x,model$pu[,der=i,jnf[j]],lty=2,col=ICcol,type=ICtype)
		}}}}}