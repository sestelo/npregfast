#' @export globaltest

globaltest <-
function(formula,data=data,der=NULL,weights=NULL,nboot=200,h=-1.0,nh=30,kernel="epanech",p=3,kbin=100){
	if(missing(der)){
 		 stop("Argument \"der\" is missing, with no default")
 	 }

	if(missing(formula)){
 		 stop("Argument \"formula\" is missing, with no default")
 	 }
 	 if(missing(data)){
 		 stop("Argument \"data\" is missing, with no default")
 	 }
	
	if(kernel=="gaussian")  kernel=3
	if(kernel=="epanech")   kernel=1
	if(kernel=="triang")    kernel=2
	
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
		
	globaltest	<-.Fortran("globaltest",
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
		r	   = as.integer(der),
		T	   = as.double(rep(-1.0,1)),
		pvalor = as.double(rep(-1.0,1))
		)	
		#res	<-	list(pvalue=globaltest$pvalor,
		#			 T=globaltest$T)
if(globaltest$pvalor<0.05){ decision="Rejected"}else{decision="Acepted"}
		res=data.frame(cbind(Statistic=globaltest$T,pvalue=globaltest$pvalor),Decision=I(decision))
		#res=cbind(Statistic=round(globaltest$T,digits=4),pvalue=round(globaltest$pvalor,digits=4),Decision=I(decision))		
		#res=as.numeric(res)
		#res=as.data.frame(res)
		#class(res) <- "globaltest"
	return(res)
		
}
