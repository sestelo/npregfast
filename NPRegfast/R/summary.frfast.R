summary.frfast <-
function(model)
{
	if (missing(model)) 
        stop("Argument 'model' is missing with no default")
        
	 if(model$modelo==1){
	 
    		cat("\nCall:\n")
	 		print(model$call)
	 		cat("", "\n")
     		if(model$modelo==1){m="Nonparametric"}else{m="Alometric"}
     		cat("*********************************************", "\n")
  		   cat(m,"Model","\n")
  		  cat("*********************************************", "\n")
  		  if (model$kernel==1) cat("Kernel: Epanechnikov \n")
			if (model$kernel==2) cat("Kernel: Triangular \n")
			if (model$kernel==2) cat("Kernel: Gaussian \n")
   		 cat("Bandwidth:", model$h,"\n")
   		 cat("Degree of polinomium:", model$dp,"\n") 
   		 cat("Number of bootstrap repeats:", model$nboot,"\n")
    	cat("Number of binning nodes", model$grid,"\n")   
   		 cat("", "\n")
   		 cat("", "\n")
   		 cat("The number of data is: ",model$n, "\n")
    	cat("The factor's levels are: ",etiquetas<-model$etiquetas, "\n")
    
   		 nf<-length(etiquetas)
   		 if (nf != 1){  	
    		for (factor in c(1:nf)){
    			  cat("The number of data for the level", etiquetas[factor], "is:",length(model$xdata[model$fmod==factor]), "\n")}  
    	  }
    	    
    	 cat("", "\n") 
    	 cat("Summaries for the variable y (for each level): ")  
    	  
    for (factor in c(1:nf)){
    	cat("", "\n")
    	cat("Level", etiquetas[factor],":", "\n")
    	
    	print(summary(model$ydata[model$fmod==factor]))} 	  
  }
        
        
if(model$modelo==2){
	etiquetas<-model$etiquetas
	   nf<-length(etiquetas)

  	cat("\nCall:\n")
	 print(model$call)
	 cat("", "\n")
     if(model$modelo==1){m="Nonparametric"}else{m="Allometric"}
     cat("*********************************************", "\n")
     cat(m,"Model","\n")
    cat("*********************************************", "\n")
    
    
    
    
    cat("", "\n")
    if(nf==1){cat("Coefficients:", "\n")}else{cat("Coefficients (for each level):", "\n")}
    cat("", "\n")
     for (factor in c(1:nf)){
    	if(nf>1)	cat("Level", etiquetas[factor],":", "\n")
    	aux=c(model$a[factor],model$a_inf[factor],model$a_sup[factor],model$b[factor],model$b_inf[factor],model$b_sup[factor])
    	aux=round(aux,6)
    	tabla=matrix(aux,ncol=3,nrow=2,byrow=T)
    	colnames(tabla)=c("","2.5 %","97.5 %")
    	rownames(tabla)=c("a","b")
    	print(tabla)
    	cat("", "\n")
    	}
    cat("*********************************************", "\n")
    
    
    cat("", "\n")
    cat("", "\n")
    cat("Degree of polinomium:", model$dp,"\n") 
    cat("Number of bootstrap repeats:", model$nboot,"\n")
    cat("Number of binning nodes", model$grid,"\n")   
    cat("", "\n")
    cat("", "\n")
    cat("The number of data is: ",model$n, "\n")
    cat("The factor's levels are: ",etiquetas<-model$etiquetas, "\n")
    

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
