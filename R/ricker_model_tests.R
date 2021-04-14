
# testing the ricker model with carrying capacity
library(tidyverse)

n.sp <- 10
# k <- rep(50,n.sp)

lambda <- runif(10,1,10)
names(lambda) <-  paste("sp",1:n.sp,sep="")

init.ab <- sample(1:100,n.sp,TRUE)
names(init.ab) <- paste("sp",1:n.sp,sep="")

alpha.matrix <- matrix(runif(n.sp*n.sp,0,0.01),nrow = n.sp,dimnames = list(names(lambda),names(lambda)))
diag(alpha.matrix) <- runif(n.sp,0.01,0.1)

timesteps <- 100

ab <- data.frame(timestep = 1:timesteps)
for(i.sp in 1:n.sp){
  ab[,ncol(ab)+1] <- rep(NA_real_,timesteps)
}
names(ab)[2:ncol(ab)] <- names(lambda)
ab[1,2:(n.sp+1)] <- init.ab
  
for(i.timestep in 2:timesteps){
  for(i.sp in 1:n.sp){
    focal.lambda <- lambda[i.sp]
    alpha_intra <- alpha.matrix[i.sp,i.sp]
    names(alpha_intra) <- names(lambda)[i.sp]
    alpha_inter <- alpha.matrix[i.sp,-i.sp]
    
    abund_intra <- ab[i.timestep-1,i.sp+1]
    abund_inter <- ab[i.timestep-1,2:(n.sp+1)]
    abund_inter <- as.numeric(abund_inter[-i.sp])
    
    # model
    alpha_inter <- alpha_inter * -1
    alpha_intra <- alpha_intra * -1
    
    term <- 0
    for(z in 1:length(alpha_inter)){
      term <- term + alpha_inter[z]*abund_inter[z]
    }
    
    intra_term <- alpha_intra*abund_intra
    term <- term + intra_term

    temp.abund <- (focal.lambda + term) * abund_intra
    
    ab[i.timestep,i.sp+1] <- ifelse(temp.abund < 0, 0, temp.abund)
      
    # ab[i.timestep,i.sp+1] <- (focal.lambda + term) *abund_intra
    # ab[i.timestep,i.sp+1] <- abund_intra + focal.lambda * abund_intra * ((k[i.sp] + abund_intra)/k[i.sp]) 
    
  }# for i.sp
}# for i.timestep

ab.long <- pivot_longer(ab,2:(n.sp+1))

ggplot(ab.long, aes(x = timestep, y = value)) + 
  geom_point(aes(color = name)) +
  NULL



