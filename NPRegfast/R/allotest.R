#li### opcion cambiar bootstrap y argumentos... ver
allotest <-
function(formula,data=data,nboot=100){
	
	
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
	kbin=100	
				
		res=list()
		for(i in etiquetas){
			yy=data[,1][f==i]
			xx=data[,2][f==i]
			w=rep(1,n)
			#datatest=cbind(yy,xx)
			#formula=yy~xx
			fit<-.Fortran("test_allo",
				x=as.double(xx),
				y=as.double(yy),
				w=as.double(w),
				n=as.integer(n),
				kbin=as.integer(kbin),
				nboot=as.integer(nboot),
				T=as.double(-1.0),
				pvalue = as.double(-1.0)
				)	
				
			#fit<-frfast(formula,data=datatest,model=0,nboot=nboot)
			res[[i]]<- list(statistic=fit$T, 
			pvalue=fit$pvalue)	
		}
		
		
		est=c()
		p=c()
		for(i in 1:length(res)){
			est[i]=res[[i]]$statistic
			p[i]=res[[i]]$pvalue	
		}
			
	kk=cbind(est,p)
	result<-matrix(kk,ncol=2,nrow=length(res))
	colnames(result)=c("Statistic","pvalue")
		
		
		#factores=paste("Factor",1:length(res))
		factores=paste("Level",etiquetas[1:length(etiquetas)])
		
			
		for(i in 1:length(res)){
			rownames(result)=c(factores)
		}
			
	return(result)
		
}
