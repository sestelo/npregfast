#' @export localtest

localtest<-
function(formula,data=data,der=NULL,weights=NULL,nboot=200,h=-1.0,nh=30,kernel=1,p=3,kbin=100,ranku=NULL, rankl=NULL){
	if(missing(der)){
 		 stop("Argument \"der\" is missing, with no default")
 	 }
	if(missing(formula)){
 		 stop("Argument \"formula\" is missing, with no default")
 	 }
 	 if(missing(data)){
 		 stop("Argument \"data\" is missing, with no default")
 	 }
	
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
 	if(is.null(h)){h <- rep(-1.0,nf)
 		}else{ h<-rep(h,nf)}#1 h para cada localidad. 
 							#Interesaria meter !=h para las != localidades, y para las derivadas?
	if(is.null(weights)) {
 		weights <- rep(1.0, n)
 	}else{
 		if(sum(weights)<=0 || any(weights)<0 || length(weights)!= n) 
 			stop("The specified weights are not correct")}
	if(is.null(rankl)) rankl<-as.vector(tapply(data[,varnames],f,min))# si rank son2fact y meto1num casca!!! corregir con repeat. 
 	if(is.null(ranku)) ranku<-as.vector(tapply(data[,varnames],f,max))
 	
 	
	localtest	<-.Fortran("localtest",
		f      = as.integer(f),
		x      = as.double(data[,varnames]),
		y      = as.double(data[,ffr$response]),
		w      = as.double(weights),
		n      = as.integer(n),
		h      = as.double(h),
		nh	   = as.integer(nh),
		p      = as.integer(p),
		kbin   = as.integer(kbin),
		fact   = as.integer(c(1:nf)), #fact   =as.integer(c(1:nf))
		nf     = as.integer(nf),
		kernel = as.integer(kernel),
		nboot  = as.integer(nboot),
		pcmax =as.double(ranku), # rango de busqueda minimo
		pcmin =as.double(rankl), # rango de busqueda maximo
		r	   = as.integer(der),
		D	   = as.double(rep(-1.0,1)),
		Ci = as.double(rep(-1.0,1)),
		Cs = as.double(rep(-1.0,1))
		)	
		
		if(localtest$Ci<=0&0<=localtest$Cs){ decision="Acepted"
			}else{decision="Rejected"}
		res=cbind(D=round(localtest$D,digits=4),Lwr=round(localtest$Ci,digits=4),Upr=round(localtest$Cs,digits=4),Decision=decision)
	#class(res) <- "localtest"
	return(as.data.frame(res))
		
}
