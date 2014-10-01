test <-
function(x,y,f=NULL){
	n    <- length(x)
	if(is.null(f)) f <- rep(1.0,n)
	
	fact<-unique(f)
 	nf<-length(fact)
 
		res=list()
		for(i in 1:nf){
			fit<-frfast(x=x[f==i],y=y[f==i],model=0)
			res[[i]]<- list(stadistic=fit$stadistic, 
			pvalue=fit$pvalue)	
		}
		
		
		est=c()
		p=c()
		for(i in 1:length(res)){
			est[i]=res[[i]]$stadistic
			p[i]=res[[i]]$pvalue	
		}
			
	kk=cbind(est,p)
	result<-matrix(kk,ncol=2,nrow=length(res))
	colnames(result)=c("Statistic","p-value")
		
		
		factores=paste("Level",1:length(res))
			
		for(i in 1:length(res)){
			rownames(result)=c(factores)
		}
			
	return(result)
		
}

