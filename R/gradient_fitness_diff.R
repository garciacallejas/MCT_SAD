# returns a list with lambda values
gradient_fitness_diff <- function(lambda, A = NULL, steps){
  
  res <- list()
  
  mean.v <- mean(lambda)
  final.v <- rep(mean.v,length(lambda))
  
  if(steps > 1){
    v_director = final.v - lambda
    for (i in 1:steps){
      res[[i]] <- lambda + (i-1)*v_director/(steps-1)
    }
  }else{
    res <- final.v
  }
  
  return(res)
  
}
