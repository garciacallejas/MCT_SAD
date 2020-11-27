
# returns a list of matrices, with decreasing pairwise asymmetry
gradient_asymmetry <- function(lambda = NULL,A,steps){
  
  res <- list()

  # go step by step
  for(i.step in 1:steps){
    
    step.matrix <- A
    
    # each position
    for(i in 1:nrow(A)){
      for(j in 1:ncol(A)){
        my.gradient <- seq(from = A[i,j],to = mean(c(A[i,j],A[j,i])),
                           length.out = steps)
        step.matrix[i,j] <- my.gradient[i.step]
      }
    }
    
    res[[i.step]] <- step.matrix
    
  }# for each step
  
  return(res)
  
}