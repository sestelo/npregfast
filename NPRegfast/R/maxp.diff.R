maxp.diff <-
function(model,factor1=NULL,factor2=NULL,der=NULL){
	nf<-length(unique(model$fmod))
	model$diffmax[model$diffmax==9999]=NA
	model$diffmaxl[model$diffmaxl==9999]=NA
	model$diffmaxu[model$diffmaxu==9999]=NA

a=t(matrix(combn(nf,2), nrow = 2)) 
nrow(a)	
res=list()

if(is.null(der)& is.null(factor2)&is.null(factor1)){

		for(i in 1:nrow(a)){	
			res[i]=list(matrix(ncol=5,nrow=3))
			colnames(res[[i]])=c("Factor2","Factor1","Max point","Lwr","Upr")
			rownames(res[[i]])=c("Estimation","First_der","Second_der")	
						for(k in 1:3){
								res[[i]][k,1]=a[i,2]
								res[[i]][k,2]=a[i,1]
								res[[i]][k,3]= round(c(model$diffmax[k,a[i,1],a[i,2]]),3)
								res[[i]][k,4]= round(c(model$diffmaxl[k,a[i,1],a[i,2]]),3)
								res[[i]][k,5]= round(c(model$diffmaxu[k,a[i,1],a[i,2]]),3)					
								}			
			}
		
		return(res)
		
}else if(is.null(der)){
	
			res=matrix(ncol=5,nrow=3)	
			
						for(k in 1:3){
								res[k,1]=factor2
								res[k,2]=factor1
								if(factor2<factor1) {fac2=factor2;fac1=factor1;factor2=fac1;factor1=fac2}else{fac2=factor2;fac1=factor1}
								res[k,3]= if(fac2<fac1){-1*round(c(model$diffmax[k,factor1,factor2]),3)}else{round(c(model$diffmax[k,factor1,factor2]),3)}
								res[k,4]= if(fac2<fac1){-1*round(c(model$diffmaxl[k,factor1,factor2]),3)}else{round(c(model$diffmaxl[k,factor1,factor2]),3)}
								res[k,5]= if(fac2<fac1){-1*round(c(model$diffmaxu[k,factor1,factor2]),3)}else{round(c(model$diffmaxu[k,factor1,factor2]),3)}
								if(fac2<fac1){factor2=fac2;factor1=fac1}
								}			
		
		colnames(res)=c("Factor2","Factor1","Max points Diff.","Lwr","Upr")
		rownames(res)=c("Estimation","First_der","Second_der")
		return(res)
		
		
		
}else if(is.null(factor2)&is.null(factor1)){
	der=der+1
	for(i in 1:nrow(a)){	
		res[i]=list(matrix(ncol=5,nrow=1))
		
		colnames(res[[i]])=c("Factor2","Factor1","Max points Diff.","Lwr","Upr")
		if(der==1)		rownames(res[[i]])= c("Estimation")
		if(der==2)		rownames(res[[i]])=c("First_der")
		if(der==3)		rownames(res[[i]])=c("Second_der")
		res[[i]][1,1]=a[i,2]
		res[[i]][1,2]=a[i,1]
		res[[i]][1,3]=round(c(model$diffmax[der,a[i,1],a[i,2]]),3)
		res[[i]][1,4]=round(c(model$diffmaxl[der,a[i,1],a[i,2]]),3)
		res[[i]][1,5]=round(c(model$diffmaxu[der,a[i,1],a[i,2]]),3)
			}		
		return(res)

	}else{
			der=der+1
			res=matrix(ncol=5,nrow=1)
			res[1,1]=factor2
			res[1,2]=factor1
			if(factor2<factor1) {fac2=factor2;fac1=factor1;factor2=fac1;factor1=fac2}else{fac2=factor2;fac1=factor1}
			res[1,3]= if(fac2<fac1){-1*round(c(model$diffmax[der,factor1,factor2]),3)}else{round(c(model$diffmax[der,factor1,factor2]),3)}
			res[1,4]= if(fac2<fac1){-1*round(c(model$diffmaxl[der,factor1,factor2]),3)}else{round(c(model$diffmaxl[der,factor1,factor2]),3)}
			res[1,5]= if(fac2<fac1){-1*round(c(model$diffmaxu[der,factor1,factor2]),3)}else{round(c(model$diffmaxu[der,factor1,factor2]),3)}
			
								
			colnames(res)=c("Factor2","Factor1","Max points Diff.","Lwr","Upr")
			if(der==1)		rownames(res)= c("Estimation")
			if(der==2)		rownames(res)=c("First_der")
			if(der==3)		rownames(res)=c("Second_der")
			return(res)
		
		}
}

