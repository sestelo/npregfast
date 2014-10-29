#' @export predict.frfast

predict.frfast<-function (model, newdata, fac = NULL, der = NULL,...){
	### tiene que meter a huevo el newdata!!! asi que meter un warning
	newdata<-newdata[,1]
	len<-length(newdata)
	f<-model$fmod
	f<-c(f,rep(1,length(newdata)))	
	x<-model$xdata
	x<-c(x,newdata)
	y<-model$ydata
	y<-c(y,rep(-1.0,length(newdata)))	
	n<-length(x)
	h<-model$h
	weights<-c(model$w,rep(0,length(newdata)))
	p<-model$dp
	kbin<-model$kbin
	ncmax<-5
	ikernel<-1
	iopt<-1
	nboot<-model$nboot
	rankl<-model$rankl
	ranku<-model$ranku
	kernel<-model$kernel
	nh<-model$nh
	nf<-model$nf
	nc <- rep(-1.0,nf)
	c2<-matrix(as.double(-1.0),ncmax,nf) 
	tmodel<-model$nmodel
	ipredict2=1
 frfast  <- .Fortran("frfast",
       f      = as.integer(f),
       x      = as.double(x),
       y      = as.double(y),
       w      = as.double(weights),
       n      = as.integer(n),
       h      = as.double(h),
       c2     = as.integer(c2),
       nc     =as.integer(nc),
       ncmax  =as.integer(ncmax),
       p      =as.integer(p),
       kbin   =as.integer(kbin),
       fact   =as.integer(c(1:nf)),
       nf     =as.integer(nf),
       ikernel=as.integer(ikernel),
       iopt   =as.integer(iopt),
       nboot  =as.integer(nboot),
       xb     = as.double(rep(-1.0,kbin)),
       pb     = array(rep(-1.0),c(kbin,3,nf)),
       li     = array(as.double(-1.0),c(kbin,3,nf)),
       ls     = array(as.double(-1.0),c(kbin,3,nf)),
       dif    = array(as.double(-1.0),c(kbin,3,nf,nf)),
       difi   = array(as.double(-1.0),c(kbin,3,nf,nf)),
       difs	  = array(as.double(-1.0),c(kbin,3,nf,nf)),
       tmodel  = as.integer(tmodel), 
       pvalor = as.double(rep(-1.0,1)),
       c 	  = array(as.double(-1.0),c(3,nf)),
       cs	  = array(as.double(-1.0),c(3,nf)),
       ci	  = array(as.double(-1.0),c(3,nf)),
       difc	  = array(as.double(-1.0),c(3,nf,nf)),
       difcs  = array(as.double(-1.0),c(3,nf,nf)),
       difci  = array(as.double(-1.0),c(3,nf,nf)),
		T 	  = as.double(rep(-1.0,1)),
		pboot =array(as.double(-1.0),c(kbin,3,nf,nboot)),
		pcmin =as.double(rankl), # rango de busqueda minimo
		pcmax =as.double(ranku), # rango de busqueda maximo
		cboot =array(as.double(-1.0),c(3,nf,nboot)), ## max de las bootstrap
		kernel=as.integer(kernel),
		nh    =as.integer(nh),
		a     =as.double(rep(-1.0,nf)),
		ainf  =as.double(rep(-1.0,nf)),
		asup  =as.double(rep(-1.0,nf)),
		b     =as.double(rep(-1.0,nf)),
		binf  =as.double(rep(-1.0,nf)),
		bsup  =as.double(rep(-1.0,nf)),
		ipredict=as.integer(ipredict2),
		predict=array(as.double(-1.0),c(n,3,nf)),
		predictl=array(as.double(-1.0),c(n,3,nf)),
		predictu=array(as.double(-1.0),c(n,3,nf))
	)
frfast$predict[frfast$predict==-1]=NA
frfast$predictl[frfast$predictl==-1]=NA
frfast$predictu[frfast$predictu==-1]=NA



ii=c(rep(FALSE,n-length(newdata)),rep(TRUE,length(newdata)))
#v=vector("list",length=nf)
if(nf==1){
	for(k in 1:nf){
		if(is.null(der)){der=c(0,1,2)}
			cont=0
			der=der+1
			res=array(data=NA,dim=c(len,3,length(der)))
			for(j in der){
				cont=cont+1
				res[,1,cont]=frfast$predict[,,1][ii,j]
				res[,2,cont]=frfast$predictl[,,1][ii,j]
				res[,3,cont]=frfast$predictu[,,1][ii,j]
				}		
		}
		colnames(res)=c("Predict","Lwr","Upr")
		if(length(der)==3) {res<-list(Estimation=res[,,1],
				  First_deriv=res[,,2],
				  Second_deriv=res[,,3])}else{
				  
				  if(der==1) {res<-list(Estimation=res[,,1])}
				 	if(der==2) {res<-list(First_deriv=res[,,1])}
				 	if(der==3) {res<-list(Second_deriv=res[,,1])} }
 
		class(res)<-"predict.frfast"
		return(res)
}else{
	if(is.null(fac)) fac=unique(f)
	factores=c()
	resul=vector("list",length=length(fac))
	zz=1
	for(k in fac){
		factores[zz]=paste("Level_",frfast$fact[k],sep="")
		zz=zz+1}
	names(resul)<-factores
		if(is.null(der)){der=c(0,1,2)}
			der=der+1
			zz=1
		for(k in fac){
			cont=0
			res=array(data=NA,dim=c(len,3,length(der)))
			for(j in der){
				cont=cont+1
				res[,1,cont]=frfast$predict[,,k][ii,j]
				res[,2,cont]=frfast$predictl[,,k][ii,j]
				res[,3,cont]=frfast$predictu[,,k][ii,j]	
				}
			colnames(res)=c("Predict","Lwr","Upr")		
			# res<-list(Estimation=res[,,1],
				  # First_deriv=res[,,2],
				  # Second_deriv=res[,,3])
				  if(length(der)==3) {res<-list(Estimation=res[,,1],
				  First_deriv=res[,,2],
				  Second_deriv=res[,,3])}else{
				  
				  if(der==1) {res<-list(Estimation=res[,,1])}
				 	if(der==2) {res<-list(First_deriv=res[,,1])}
				 	if(der==3) {res<-list(Second_deriv=res[,,1])} }
	
			resul[[zz]]=res
			zz=zz+1
		}
		class(resul)<-"predict.frfast"
		return(resul)
		}
		
	}	

		
		
		
		
		
		
		
