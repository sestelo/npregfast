summary.frfast <-
function(model)
{
	if (missing(model)) 
        stop("Argument 'model' is missing with no default")
        
	 if(model$nmodel==1){
	 
    		cat("\nCall:\n")
	 		print(model$call)
	 		cat("", "\n")
     		if(model$nmodel==1){m="Nonparametric"}else{m="Alometric"}
     	cat("*********************************************", "\n")
  		cat(m,"Model","\n")
  		cat("*********************************************", "\n")
  		 if (model$kernel==1) cat("Kernel: Epanechnikov \n")
		 if (model$kernel==2) cat("Kernel: Triangular \n")
		 if (model$kernel==3) cat("Kernel: Gaussian \n")
   		 cat("Bandwidth:", model$h,"\n")
   		 cat("Degree of polinomium:", model$dp,"\n") 
   		 cat("Number of bootstrap repeats:", model$nboot,"\n")
    	 cat("Number of binning nodes", model$kbin,"\n")   
   		 cat("", "\n")
   		 cat("", "\n")
   		 cat("The number of data is: ",model$n, "\n")
    	
    	 #cat("The factor's levels are: ",etiquetas<-model$etiquetas, "\n")
    
   		 nf<-length(model$label)
   		 if (nf != 1){
   		 	cat("The factor's levels are: ",etiquetas<-model$label, "\n")  	
    		for (factor in c(1:nf)){
    			  cat("The number of data for the level", etiquetas[factor], "is:",length(model$xdata[model$fmod==factor]), "\n")
    		}
    		cat("", "\n") 
    	    cat("Summaries for the response variable (for each level): ")
    	    for (factor in c(1:nf)){
    			cat("", "\n")
    			cat("Level", etiquetas[factor],":", "\n")
		    	print(summary(model$ydata[model$fmod==factor]))
		    } 
    	    	     
    	  }else{
    	  cat("", "\n") 
    	  cat("Summaries for the response variable: \n")  
     	  print(summary(model$ydata))}
    	    
    	 	  
  }
        
        
if(model$nmodel==2){
	etiquetas<-model$label
	   nf<-length(etiquetas)


  	cat("\nCall:\n")
	 print(model$call)
	 cat("", "\n")
     if(model$nmodel==1){m="Nonparametric"}else{m="Allometric"}
     cat("*********************************************", "\n")
     cat(m,"Model","\n")
    cat("*********************************************", "\n")
    
    
    
    
    cat("", "\n")
    if(nf==1){cat("Coefficients:", "\n")}else{cat("Coefficients (for each level):", "\n")}
    cat("", "\n")
     for (factor in c(1:nf)){
    	if(nf>1)	cat("Level", etiquetas[factor],":", "\n")
    	aux=c(model$a[factor],model$al[factor],model$au[factor],model$b[factor],model$bl[factor],model$bu[factor])
    	aux=round(aux,6)
    	tabla=matrix(aux,ncol=3,nrow=2,byrow=T)
    	colnames(tabla)=c("","2.5 %","97.5 %")
    	rownames(tabla)=c("a","b")
    	print(tabla)
    	cat("", "\n")
    	}
    	cat("Adjusted R-squared: ",model$r2, "\n")
    	cat("", "\n")
    cat("*********************************************", "\n")
    
    
    cat("", "\n")
    cat("", "\n")
    cat("Degree of polinomium:", model$dp,"\n") 
    cat("Number of bootstrap repeats:", model$nboot,"\n")
    cat("Number of binning nodes", model$kbin,"\n")   
    cat("", "\n")
    cat("", "\n")
    cat("The number of data is: ",model$n, "\n")
    cat("The factor's levels are: ",etiquetas<-model$label, "\n")
    

    if (nf != 1){	
    for (factor in c(1:nf)){
    	  cat("The number of data for the level", etiquetas[factor], "is:",length(model$xdata[model$fmod==factor]), "\n")}
    	  }
    	    
    	 cat("", "\n") 
    if(nf>1){cat("Summaries for the variable y (for each level): ")}else{"Summaries for the variable y:" }
    	  
  for (factor in c(1:nf)){
    	cat("", "\n")
    if(nf>1)	cat("Level", etiquetas[factor],":", "\n")
    	
    	print(summary(model$ydata[model$fmod==factor]))} 	  

 }
}
