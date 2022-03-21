
# returns a list of matrices, with decreasing pairwise asymmetry
gradient_asymmetry <- function(A,steps,
                               min.value = 0.01,
                               max.value = 0.99,
                               decreasing = TRUE){
  
  res <- list()

  if(decreasing){
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
  
  }else{
    
    for(i.step in 1:steps){
      
      step.matrix <- A
      
      # each pair
      for(i in 1:nrow(A)){
        for(j in (i+1):ncol(A)){
          
          my.low <- min(c(A[i,j],A[j,i]))
          my.high <- max(c(A[i,j],A[j,i]))
          
          low.gradient <- seq(from = my.low,to = min.value,
                              length.out = steps)
          high.gradient <- seq(from = my.high, to = max.value,
                               length.out = steps)
          if(A[i,j]<A[j,i]){
            step.matrix[i,j] <- low.gradient[i.step]
            step.matrix[j,i] <- high.gradient[i.step]
          }else{
            step.matrix[i,j] <- high.gradient[i.step]
            step.matrix[j,i] <- low.gradient[i.step]
          }
          
        }# for j
      }# for i
      
      res[[i.step]] <- step.matrix
      
    }# for each step
    
  }# if-else decreasing
    
  return(res)
  
}