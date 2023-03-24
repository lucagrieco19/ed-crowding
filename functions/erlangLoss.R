erlangLoss <- function(l,m,c){
	
	res <- l

	if((m*c)>0){
		num <- ((l/m)^c)/factorial(c)
		den <- sum(sapply(0:c,function(x){((l/m)^x)/factorial(x)}))
		res <- num / den
	}

	return(res)

}