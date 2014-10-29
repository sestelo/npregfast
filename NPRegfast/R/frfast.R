#' @useDynLib NPRegfast
#' @export frfast



frfast <-
function(formula,data=data,model="np",h=-1.0,nh=30,weights=NULL,kernel="epanech",p=3,kbin=100,nboot=500,rankl=NULL,ranku=NULL){	
	
	if(kernel=="gaussian")  kernel=3
	if(kernel=="epanech")   kernel=1
	if(kernel=="triang")    kernel=2
	##for (i in 1:nf){
 	##		if(is.null(pcmin))  pcmin<-rep(min(x[f==i]),nf) ## 
 	##		if(is.null(pcmax)) pcmax <-rep(60,nf)}  ## ver esto... 
 		
 	 if(missing(formula)){
 		 stop("Argument \"formula\" is missing, with no default")
 	 }
 	 if(missing(data)){
 		 stop("Argument \"data\" is missing, with no default")
 	 }
 	 if(!(kernel %in% 1:3)){
 		 stop("Kernel not suported")
 	}
		
	ikernel=1
	nc=NULL
	ncmax=5
	iopt=1
	c2=NULL

	
	ffr <- interpret.frfastformula(formula, method = "frfast")
	varnames <- ffr$II[2,]
	aux<-unlist(strsplit(varnames,split=":"))
	varnames<-aux[1]
	namef<-aux[2]
	if(length(aux)==1){f=NULL}else{f<-data[,namef]}
	newdata <- data
	data <- na.omit(data[,c(ffr$response, varnames)])
	newdata <- na.omit(newdata[,varnames])
	n=nrow(data)
	
	if(is.null(f)) f <- rep(1.0,n)
	etiquetas<-unique(f)
 	nf<-length(etiquetas)
 #	asfactorf<-as.factor(f)
 #	factorf<-factor(asfactorf,levels=etiquetas)
 #	f<-as.numeric(factorf)
 		 	
	if(model=="np") tmodel=1
	if(model=="allo") tmodel=2
	if(model==0) tmodel=0
	
 	if(is.null(h)){h <- rep(-1.0,nf)
 		}else{ h<-rep(h,nf)}#1 h para cada localidad. 
 							#Interesaria meter !=h para las != localidades, y para las derivadas?
 	if(is.null(nc)) nc <- rep(-1.0,nf)
 	if(is.null(weights)) {
 		weights <- rep(1.0, n)
 	}else{
 		if(sum(weights)<=0 || any(weights<0) || length(weights)!= n) 
 			stop("The specified weights are not correct")
 	}   			
 	if(is.null(c2)) c2<-matrix(as.double(-1.0),ncmax,nf) 
 	if(is.null(rankl)) rankl<-as.vector(tapply(data[,varnames],f,min))# si rank son2fact y meto1num casca!!! corregir con repeat. 
 	if(is.null(ranku)) ranku<-as.vector(tapply(data[,varnames],f,max))
 
 	ipredict2=0

 	
 		
 frfast  <- .Fortran("frfast",
       f      = as.integer(f),
       x      = as.double(data[,varnames]),
       y      = as.double(data[,ffr$response]),
       w      = as.double(weights),
       n      = as.integer(n),
       h      = as.double(h),
       c2     = as.integer(c2),
       nc     =as.integer(nc),
       ncmax  =as.integer(ncmax),
       p      =as.integer(p),
       kbin   =as.integer(kbin),
       fact   =as.integer(c(1:nf)), #fact   =as.integer(c(1:nf))
       nf     =as.integer(nf),
       ikernel=as.integer(ikernel),
       iopt   =as.integer(iopt),
       nboot  =as.integer(nboot),
       xb     = as.double(rep(-1.0,kbin)),
       pb     = array(rep(-1.0),c(kbin,3,nf)),
       li     = array(as.double(-1.0),c(kbin,3,nf)),
       ls     = array(as.double(-1.0),c(kbin,3,nf)),
       dif     = array(as.double(-1.0),c(kbin,3,nf,nf)),
       difi    = array(as.double(-1.0),c(kbin,3,nf,nf)),
       difs	= array(as.double(-1.0),c(kbin,3,nf,nf)),
       tmodel	= as.integer(tmodel), 
       pvalor	= as.double(rep(-1.0,1)),
       c 		= array(as.double(-1.0),c(3,nf)),
       cs		= array(as.double(-1.0),c(3,nf)),
       ci		= array(as.double(-1.0),c(3,nf)),
       difc	= array(as.double(-1.0),c(3,nf,nf)),
       difcs	= array(as.double(-1.0),c(3,nf,nf)),
       difci	= array(as.double(-1.0),c(3,nf,nf)),
		T 		= as.double(rep(-1.0,1)),
		pboot	=array(as.double(-1.0),c(kbin,3,nf,nboot)),
		pcmin =as.double(rankl), # rango de busqueda minimo
		pcmax =as.double(ranku), # rango de busqueda maximo
		cboot=array(as.double(-1.0),c(3,nf,nboot)), ## max de las bootstrap
		kernel=as.integer(kernel),
		nh=as.integer(nh),
		a     =as.double(rep(-1,nf)),
		ainf     =as.double(rep(-1,nf)),
		asup     =as.double(rep(-1,nf)),
		b     =as.double(rep(-1,nf)),
		binf     =as.double(rep(-1,nf)),
		bsup     =as.double(rep(-1,nf)),
		ipredict=as.integer(ipredict2),
		predict=array(rep(-1.0),c(kbin,3,nf)),
		predictl=array(as.double(-1.0),c(kbin,3,nf)),
		predictu=array(as.double(-1.0),c(kbin,3,nf)),
    PACKAGE="NPRegfast"
		)

if(tmodel!=2){
	frfast$a=NULL
	frfast$ainf=NULL
	frfast$asup=NULL
	frfast$b=NULL
	frfast$binf=NULL
	frfast$bsup=NULL
	r2=NULL
	}

#R-squared
if(tmodel==2){
	yhat=frfast$a*(frfast$x^frfast$b)
	rss=sum( (frfast$y-yhat)**2 ) / (frfast$n-2)
	tss=sum(  (frfast$y-mean(frfast$y))**2 ) / (frfast$n-1)
	r2=1-(rss/tss)
}




frfast$pb[frfast$pb==-1]=NA
frfast$li[frfast$li==-1]=NA
frfast$ls[frfast$ls==-1]=NA
frfast$dif[frfast$dif==-1]=NA
frfast$difi[frfast$difi==-1]=NA
frfast$difs[frfast$difs==-1]=NA
frfast$c[frfast$c==-1]=NA
frfast$cs[frfast$cs==-1]=NA
frfast$ci[frfast$ci==-1]=NA
frfast$difc[frfast$difc==-1]=NA
frfast$difcs[frfast$difcs==-1]=NA
frfast$difci[frfast$difci==-1]=NA
frfast$pvalor[frfast$pvalor==-1]=NA #"Argument model has to be 0"
frfast$pboot[frfast$pboot==-1]=NA	
	
#for(i in 1:nf){
#colnames(model$p)=c("Estimation","FirstDer","SecondDer")
#	}
	
	
	
	res<- list(x=frfast$xb, 
		p=frfast$pb, 
		pl=frfast$li,
		pu=frfast$ls,
		diff=frfast$dif,
		diffl=frfast$difi,
		diffu=frfast$difs,
		nboot=frfast$nboot,
		n=frfast$n,
		dp=frfast$p,
		h=frfast$h,
		#grid=frfast$kbin, #lo comente pq estaba repetido
		fmod=frfast$f,
		xdata=as.vector(frfast$x),
		ydata=frfast$y,
		w=frfast$w,
		#fact=fact,  # Lo tuve que comentar pq me daba error
		#c2=frfast$c2, #no hace falta sacarlo
		#ncmax=frfast$ncmax, #no hace falta sacarlo
		#nc=frfast$nc,   #no hace falta sacarlo
		kbin=frfast$kbin,
		nf=frfast$nf,
		#ikernel=frfast$ikernel,
		max=frfast$c, #maximo rep boot
		maxu=frfast$cs, 
		maxl=frfast$ci,
		#maxboot=frfast$cboot,  #no hace falta sacarlo
		diffmax=frfast$difc,
		diffmaxu=frfast$difcs,
		diffmaxl=frfast$difci,
		statistic=frfast$T,
		repboot=frfast$pboot,	
		rankl=frfast$pcmin,
		ranku=frfast$pcmax,
        nmodel=frfast$tmodel, #no hace falta sacarlo
        label=as.character(etiquetas),
        numlabel=unique(frfast$f),
        kernel=frfast$kernel, 
        a=frfast$a,
        a_inf=frfast$ainf,
        a_sup=frfast$asup,
        b=frfast$b,
        b_inf=frfast$binf,
        b_sup=frfast$bsup,
        name=c(ffr$response,varnames),
        formula=formula,
        nh=frfast$nh,
        r2=r2,
        pvalue=frfast$pvalor,
        call=match.call())
		
		if(tmodel==0) res=res[-(length(res)-2)]
		
		class(res) <- "frfast"
		return(res)


}

