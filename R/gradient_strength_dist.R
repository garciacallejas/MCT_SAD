# returns a list of matrices with decreasing heterogeneity in their elements
gradient_strength_dist <- function(lambda = NULL,A,steps){
  
  mean.interaction <- mean(A)
  
  res <- list()
  
  for(i.step in 1:steps){
    step.matrix <- A
    
    for(i in 1:nrow(A)){
      for(j in 1:ncol(A)){
        my.gradient <- seq(from = A[i,j],to = mean.interaction,
                           length.out = steps)
        step.matrix[i,j] <- my.gradient[i.step]
      }
    }
    
    res[[i.step]] <- step.matrix
    
  }# for i.step
  
  return(res)
  
}