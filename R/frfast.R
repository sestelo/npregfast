#  PROBLEMAS EN LA ESTIMACION DE LOS COEFICIENTES DEL MODELO ALOMÉTRICO!!! VER



#' Fitting nonparametric models
#' 
#' @description This function is used to fit nonparametric models by
#' using local linear kernel smoothers. Additionally, a parametric 
#' model (allometric model) can be estimated.
#' @param formula An object of class \code{formula}: a sympbolic 
#' description of the model to be fitted. The details of model 
#' specification are given under 'Details'.
#' @param data A data frame or matrix containing the model response
#' variable and covariates required by the \code{formula}.
#' @param model Type model used: \code{model = 'np'}  nonparametric
#' regression model with local linear kernel smoothers, 
#' \code{model = 'allo'} the  allometric model.
#' @param h The kernel bandwidth smoothing parameter. Large values of
#' bandwidth make smoother estimates, smaller values of bandwidth make
#' less smooth estimates. The default is a bandwidth compute by 
#' cross validation.
#' @param nh Integer number of equally-spaced bandwidth on which the
#' \code{h} is discretised, to speed up computation.
#' @param weights Prior weights on the data.
#' @param kernel Character which determines the smoothing kernel. 
#' By default \code{kernel = 'epanech' }, this is, the Epanechnikov
#' density function. Also, several types of kernel funcitons 
#' can be used:  triangular and Gaussian density function, 
#' with 'triang' and 'gaussian' term, respectively.
#' @param p Degree of polynomial used.  Its value must be the value of
#' derivative + 1. The default value is 3 due to the function
#' returns the estimation, first and second derivative.
#' @param kbin Number of binning nodes over which the function 
#' is to be estimated.
#' @param nboot Number of bootstrap repeats.
#' @param rankl Number or vector specifying the minimum value for the
#' interval at which to search the \code{x} value which maximizes the
#' estimate, first or second derivative  (for each level). The default
#' is the minimum data value.
#' @param ranku Number or vector specifying the maximum value for the
#' interval at which to search the \code{x} value which maximizes the
#' estimate, first or second derivative  (for each level). The default
#' is the maximum data value.
#' @param seed Seed to be used in the bootstrap procedure.
#' @details The models fit by \code{frfast} function are specified 
#' in a compact symbolic form. The \~ operator is basic in the formation 
#' of such models. An expression of the form \code{y ~ model}  is interpreted as 
#' a specification that the response \code{y} is modelled by a predictor 
#' specified symbolically by \code{model}. The possible terms consist of a variable 
#' name or a variable name and a factor name separated by : operator. Such a 
#' term is interpreted as the interaction of the continuous variable and the factor.
#' @return An object is returned with the following elements:
#' \item{x}{vector of values of the grid points at which model is to be estimate.}
#' \item{p}{matrix of values of the grid points at which to compute the estimate,
#' their first and second derivative.}
#' \item{pl}{lower values of  95\% confidence interval for the estimate, 
#' their first and second derivative.}
#' \item{pu}{upper values of  95\% confidence interval for the estimate, 
#' their first and second derivative.}
#' \item{diff}{differences between the estimation values of a couple of 
#' levels (i. e. level 2 - level 1). The same procedure for their first
#' and second derivative.}
#' \item{diffl}{lower values of 95\% confidence interval for the differences 
#' between the estimation values of a couple of levels. It is performed for their 
#' first and second derivative.}
#' \item{diffu}{upper values of 95\% confidence interval for the differences 
#' between the estimation values of a couple of levels. It is performed for their 
#' first and second derivative.}
#' \item{nboot}{number of bootstrap repeats.}
#' \item{n}{total number of data}
#' \item{dp}{degree of polynomial used.}
#' \item{h}{the kernel bandwidth smoothing parameter.}
#' \item{fmod}{factor's level for each data.}
#' \item{xdata}{original x values.}
#' \item{ydata}{original y values.}
#' \item{w}{weights on the data.}
#' \item{kbin}{number of binning nodes over which the function is to be estimated.}
#' \item{nf}{number of levels.}
#' \item{max}{value of covariate \code{x} which maximizes the  estimate, 
#' first or second derivative.}
#' \item{maxu}{upper value of 95\% confidence interval for the value \code{max}.}
#' \item{maxl}{lower value of 95\% confidence interval for the value \code{max}.}
#' \item{diffmax}{differences between the estimation of \code{max} for a couple of 
#' levels (i. e. level 2 - level 1). The same procedure for their first and 
#' second derivative.}
#' \item{diffmaxu}{upper value of 95\% confidence interval for the value \code{diffmax}.}
#' \item{diffmaxl}{lower value of 95\% confidence interval for the value \code{diffmax}.}
#' \item{repboot}{matrix of values of the grid points at which to compute the estimate, 
#' their first and second derivative for each bootstrap repeat.}
#' \item{rankl}{maximum value for the interval at which to search the \code{x} value 
#' which maximizes the estimate, first or second derivative  (for each level). 
#' The default is the maximum data value.}
#' \item{ranku}{minimum value for the interval at which to search the \code{x} value 
#' which maximizes the estimate, first or second derivative  (for each level). 
#' The default is the minimum data value.}
#' \item{nmodel}{type model used: \code{model = 1} the nonparametric model, 
#' \code{model = 2} the allometric model.}
#' \item{label}{labels of the variables in the model.}
#' \item{numlabel}{number of labels.}
#' \item{kernel}{character which determines the smoothing kernel.}
#' \item{a}{}
#' \item{al}{}
#' \item{au}{}
#' \item{b}{}
#' \item{bl}{}
#' \item{bu}{}
#' \item{name}{name of the variables in the model.}
#' \item{formula}{a sympbolic description of the model to be fitted.}
#' \item{nh}{}
#' \item{r2}{}


frfast <-
function(formula,data=data,model="np",h=-1.0,nh=30,weights=NULL,kernel="epanech",p=3,kbin=100,nboot=500,rankl=NULL,ranku=NULL, seed = NULL){	
	
	if(kernel=="gaussian")  kernel=3
	if(kernel=="epanech")   kernel=1
	if(kernel=="triang")    kernel=2
 		
 	 if(missing(formula)){
 		 stop("Argument \"formula\" is missing, with no default")
 	 }
 	 if(missing(data)){
 		 stop("Argument \"data\" is missing, with no default")
 	 }
 	 if(!(kernel %in% 1:3)){
 		 stop("Kernel not suported")
 	}
		
	ncmax=5
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
 		 	
	if(model=="np") tmodel=1
	if(model=="allo") tmodel=2
	#if(model==0) tmodel=0
	
 	if(is.null(h)){h <- rep(-1.0,nf)
 		}else{ h<-rep(h,nf)}#1 h para cada localidad. 
 							#Interesaria meter !=h para las != localidades, y para las derivadas?
 	
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
       ncmax  =as.integer(ncmax),
       p      =as.integer(p),
       kbin   =as.integer(kbin),
       fact   =as.integer(c(1:nf)), 
       nf     =as.integer(nf),
       nboot  =as.integer(nboot),
       xb     = as.double(rep(-1.0,kbin)),
       pb     = array(rep(-1.0),c(kbin,3,nf)),
       li     = array(as.double(-1.0),c(kbin,3,nf)),
       ls     = array(as.double(-1.0),c(kbin,3,nf)),
       dif     = array(as.double(-1.0),c(kbin,3,nf,nf)),
       difi    = array(as.double(-1.0),c(kbin,3,nf,nf)),
       difs	= array(as.double(-1.0),c(kbin,3,nf,nf)),
       tmodel	= as.integer(tmodel), 
       c 		= array(as.double(-1.0),c(3,nf)),
       cs		= array(as.double(-1.0),c(3,nf)),
       ci		= array(as.double(-1.0),c(3,nf)),
       difc	= array(as.double(-1.0),c(3,nf,nf)),
       difcs	= array(as.double(-1.0),c(3,nf,nf)),
       difci	= array(as.double(-1.0),c(3,nf,nf)),
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
		predictu=array(as.double(-1.0),c(kbin,3,nf))
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
frfast$pboot[frfast$pboot==-1]=NA	

	
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
		max=frfast$c, #maximo rep boot
		maxu=frfast$cs, 
		maxl=frfast$ci,
		#maxboot=frfast$cboot,  #no hace falta sacarlo
		diffmax=frfast$difc,
		diffmaxu=frfast$difcs,
		diffmaxl=frfast$difci,
		repboot=frfast$pboot,	
		rankl=frfast$pcmin,
		ranku=frfast$pcmax,
    nmodel=frfast$tmodel, #no hace falta sacarlo
    label=as.character(etiquetas),
    numlabel=unique(frfast$f),
    kernel=frfast$kernel, 
    a=frfast$a,
    al=frfast$ainf,
    au=frfast$asup,
    b=frfast$b,
    bl=frfast$binf,
    bu=frfast$bsup,
    name=c(ffr$response,varnames),
    formula=formula,
    nh=frfast$nh,
    r2=r2,
    call=match.call())
		
	#	if(tmodel==0) res=res[-(length(res)-2)]
		
		class(res) <- "frfast"
		return(res)


}

