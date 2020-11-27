
# returns a list with matrices
gradient_niche_diff <- function(lambda = NULL, A, steps){
  
  # auxiliary function
  removeNiches_nsp <- function(m, steps){
    pairs <- combn(colnames(m),m = 2)
    nm <- m
    
    res <- list()
    
    # for each step, generate a matrix and store it in res[[i.step]]
    for(i.step in 1:steps){
    for(i.pair in 1:ncol(pairs)){
      pairm <- m[pairs[,i.pair],pairs[,i.pair]]
      
      muIntra=sqrt(pairm[1,1]*pairm[2,2])
      muInter=sqrt(pairm[2,1]*pairm[1,2])
      rescale=muIntra/muInter
      new.alphas=pairm
      
      # weighting term
      if(steps > 1){
        w <- seq(from = 1/rescale,to = 1,length.out = steps)
      }else{
        w <- 1
      }
      
      new.alphas[2,1]=pairm[2,1]*(rescale*w[i.step])
      new.alphas[1,2]=pairm[1,2]*(rescale*w[i.step])
      
      nm[pairs[1,i.pair],pairs[2,i.pair]] <- new.alphas[pairs[1,i.pair],pairs[2,i.pair]]
      nm[pairs[2,i.pair],pairs[1,i.pair]] <- new.alphas[pairs[2,i.pair],pairs[1,i.pair]]
    }# for each pair
    
    nm[which(is.na(nm))] <- 0
    res[[i.step]] <- nm
    }# for each step
    
    return(res)
  }
  
  # generate matrices with progressively less niche diff
  # using the auxiliary function
  res <- removeNiches_nsp(A,steps)
  
  return(res)
}