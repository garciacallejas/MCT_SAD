# returns a list of matrices with decreasing heterogeneity 
# in their DIAGONAL elements: if they are dominant, they gradually converge
# towards the overall mean
# if decreasing = FALSE, the diagonal elements gradually increase up to max.value

gradient_diag_dominance <- function(A,steps,
                                    min.value = 0.01,
                                    max.value = 0.99,
                                    decreasing = TRUE){
  
  mean.interaction <- mean(A)
  
  res <- list()
  
  if(decreasing){
    
    for(i.step in 1:steps){
      step.matrix <- A
      
      for(i in 1:nrow(A)){
        my.gradient <- seq(from = A[i,i],to = mean.interaction,
                           length.out = steps)
        step.matrix[i,i] <- my.gradient[i.step]
      }
      
      res[[i.step]] <- step.matrix
      
    }# for i.step
    
  }else{
    
    for(i.step in 1:steps){
      step.matrix <- A
      
      for(i in 1:nrow(A)){
        my.gradient <- seq(from = A[i,i],to = max.value,
                           length.out = steps)
        step.matrix[i,i] <- my.gradient[i.step]
      }
      
      res[[i.step]] <- step.matrix
      
    }# for i.step
    
  }# if-else decreasing
  
  return(res)
  
}