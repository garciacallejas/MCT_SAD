homogeneize_vector <- function(v, steps = 2){
  
  res <- matrix(nrow = steps, ncol = length(v))
  
  mean.v <- mean(v)
  final.v <- rep(mean.v,length(v))
  
  if(steps > 1){
    v_director = final.v - v
    for (i in 1:steps){
      res[i,] <- v + (i-1)*v_director/(steps-1)
    }
  }else{
    res <- final.v
  }
  
  return(res)
}
