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
