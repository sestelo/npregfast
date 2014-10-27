#' @export print.frfast

print.frfast <-
function(model, ...){ #  print.frfast2(model) es igual escribir model
	cat("\nCall:\n")
	 print(model$call)
	 cat("", "\n")
     if(model$modelo==1){m="Nonparametric"}else{m="Allometric"}
     cat("*********************************************", "\n")
     cat(m,"Model","\n")
    cat("*********************************************", "\n")
		cat("\nNumber of Observations: ")
		cat(format(model$n))
		cat("\n")
		
		cat("\nNumber of Factors: ")
		cat(format(length(unique(model$fmod))))
		cat("\n")
	
		#cat("\nResidual Standar Error: ")
		#cat("\n")
		
		cat("\nNumber of Bootstrap repeats: ")
		cat(format(model$nboot))
		cat("\n")
		
		if(model$modelo==1){cat("\nBandwidth: ")
		cat(format(model$h))
		cat("\n")
		
		cat("\nKernel: ")
		if (model$kernel==1) cat("Epanechnikov")
		if (model$kernel==2) cat("Triangular")
		if (model$kernel==3) cat("Gaussian")
		}}

