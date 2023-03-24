queueSize <- function(l,m,c){
	
  if(c==0 | m==0){
    res <- Inf
  }else{
  
    r <- l/(c*m)
  
    den <- (1/factorial(c)) * ((l/m)^c) * (1/(1-r))
    for(i in 0:(c-1)) den <- den + (1/factorial(i)) * ((l/m)^i)
    
    p0 <- 1/den
  
	  res <- (1/factorial(c)) * ((l/m)^c) * (r/(1-r)^2) * p0

  }

	return(res)

}