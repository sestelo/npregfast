frfast <-
function(x,y,f=NULL,model="np",h=-1,w=NULL,p=3,kbin=100,nc=NULL,
ncmax=5,ikernel=1,iopt=1,nboot=500,c2=NULL,rankl=NULL,ranku=NULL,kernel=1,nh=30)

{
	n    <- length(x)
	if(is.null(f)) f <- rep(1.0,n) 
	etiquetas<-unique(f)
 	nf<-length(etiquetas)
 	asfactorf<-as.factor(f)
 	f<-as.numeric(asfactorf)
 	
 	
	if(model=="np") model=1
	if(model=="allo") model=2
 	if(is.null(h)) h <- rep(-1.0,nf)
 	else h<-rep(h,nf)  # meto un mismo h para cada localidad. Igual interesa meter distintos h para las distintas localidades, y para las derivadas... etc. VER QUE SE HACE
 	if(is.null(nc)) nc <- rep(-1.0,nf)
 	if(is.null(w)) w <- rep(1.0, n)  
 	if(is.null(c2)) c2<-matrix(as.double(-1.0),ncmax,nf) 
 	if(is.null(rankl)) rankl<-as.vector(tapply(x,f,min))# si rankl o rank son dos factores y meto un numero solo casca!!! corregir con repeat. 
 	if(is.null(ranku)) ranku<-as.vector(tapply(x,f,max))
 	#for (i in 1:nf){
 	#		if(is.null(pcmin))  pcmin<-rep(min(x[f==i]),nf) ### VER ESTO ...
 	#		if(is.null(pcmax)) pcmax <-rep(60,nf)}  ## ver esto... 
 


 frfast  <- .Fortran("frfast",
       f      = as.integer(f),
       x      = as.double(x),
       y      = as.double(y),
       w      = as.double(w),
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
       model	= as.integer(model), 
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
		bsup     =as.double(rep(-1,nf))
		)

if(model!=2){
	frfast$a=NULL
	frfast$ainf=NULL
	frfast$asup=NULL
	frfast$b=NULL
	frfast$binf=NULL
	frfast$bsup=NULL
	
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
		grid=frfast$kbin,
		fmod=frfast$f,
		xdata=as.vector(frfast$x),
		ydata=frfast$y,
		w=frfast$w,
		#fact=fact,  # Lo tuve que comentar pq me daba error
		c2=frfast$c2,
		ncmax=frfast$ncmax,
		nc=frfast$nc,
		kbin=frfast$kbin,
		nf=frfast$nf,
		ikernel=frfast$ikernel,
		pvalue=frfast$pvalor,
		max=frfast$c, #maximo rep boot
		maxu=frfast$cs, 
		maxl=frfast$ci,
		maxboot=frfast$cboot,
		diffmax=frfast$difc,
		diffmaxu=frfast$difcs,
		diffmaxl=frfast$difci,
		statistic=frfast$T,
		repboot=frfast$pboot,	
		rankl=frfast$pcmin,
		ranku=frfast$pcmax,
        modelo=frfast$model,
        etiquetas=as.character(etiquetas),
        etiquetasnum=unique(frfast$f),
        kernel=frfast$kernel,
        a=frfast$a,
        a_inf=frfast$ainf,
        a_sup=frfast$asup,
        b=frfast$b,
        b_inf=frfast$binf,
        b_sup=frfast$bsup,
        call=match.call())
		
		
		class(res) <- "frfast"
		return(res)


}

