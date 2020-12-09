# returns a list of matrices with decreasing heterogeneity 
# in their DIAGONAL elements: if they are dominant, they gradually converge
# towards the overall mean

gradient_diag_dominance <- function(lambda = NULL,A,steps){
  
  mean.interaction <- mean(A)
  
  res <- list()
  
  for(i.step in 1:steps){
    step.matrix <- A
    
    for(i in 1:nrow(A)){
        my.gradient <- seq(from = A[i,i],to = mean.interaction,
                           length.out = steps)
        step.matrix[i,i] <- my.gradient[i.step]
    }
    
    res[[i.step]] <- step.matrix
    
  }# for i.step
  
  return(res)
  
}