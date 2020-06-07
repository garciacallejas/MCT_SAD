# taken from Chu and Adler (2015)
# alphas is a 2 x2 matrix, with alpha_ij the effect of j on i

removeNiches<-function(alphas){
  muIntra=sqrt(alphas[1,1]*alphas[2,2])
  muInter=sqrt(alphas[2,1]*alphas[1,2])
  rescale=muIntra/muInter
  new.alphas=alphas
  new.alphas[2,1]=alphas[2,1]*rescale
  new.alphas[1,2]=alphas[1,2]*rescale
  return(new.alphas)
}

removeNiches_nsp <- function(m){
  pairs <- combn(colnames(m),m = 2)
  nm <- m
  
  for(i.pair in 1:ncol(pairs)){
    pairm <- m[pairs[,i.pair],pairs[,i.pair]]
    
    muIntra=sqrt(pairm[1,1]*pairm[2,2])
    muInter=sqrt(pairm[2,1]*pairm[1,2])
    rescale=muIntra/muInter
    new.alphas=pairm
    new.alphas[2,1]=pairm[2,1]*rescale
    new.alphas[1,2]=pairm[1,2]*rescale
    
    nm[pairs[1,i.pair],pairs[2,i.pair]] <- new.alphas[pairs[1,i.pair],pairs[2,i.pair]]
    nm[pairs[2,i.pair],pairs[1,i.pair]] <- new.alphas[pairs[2,i.pair],pairs[1,i.pair]]
  }# for each pair
  
  nm
}

