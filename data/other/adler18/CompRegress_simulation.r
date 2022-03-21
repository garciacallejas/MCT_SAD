
# compare estimates of known rho from Effect and Response comparisons

###
### 1. Set up alpha matrices and calculate rho and log ratios
###

N <- 10000   # how many pairwise comparisons
intras <- matrix(runif(N*2,0.0001,0.01),N,2)
inters <- matrix(runif(N*2,0.00005,0.005),N,2)
rawD <- data.frame(cbind(intras, inters))
names(rawD) <- c("a11","a22","a12","a21")
rawD$rho <- sqrt((rawD$a12*rawD$a21)/(rawD$a11*rawD$a22))
rawD$Effect1 <- rawD$a12/rawD$a22
rawD$Effect2 <- rawD$a21/rawD$a11
rawD$Response1 <- rawD$a12/rawD$a11
rawD$Response2 <- rawD$a21/rawD$a22

###
### 2. sample from the raw data and estimate means
###

sampleN <- seq(10,500,10)  # sample sizes to compare
trials <- 500  # number of trials to conduct at each sample size
means <- array(NA,dim=c(length(sampleN),3,trials))

for(i in 1:length(sampleN)){
  doN <- sampleN[i]
  for(j in 1:trials){
    k <- sample(1:N,doN)
    # sample rho
    means[i,1,j] <- mean(log(rawD$rho[k]))
    # sample effects
    means[i,2,j] <- mean(log(rawD$Effect1[k]))
    # sample responses
    means[i,3,j] <- mean(log(rawD$Response1[k]))
  }
}

###
### 3. compare estimates
###

means.mu <- apply(means,MARGIN=c(1,2),FUN=mean)
means.sd <- apply(means,MARGIN=c(1,2),FUN=sd)

png("simulation-fig.png",height=3,width=8,units="in",res=100)

par(mfrow=c(1,2), tcl=-0.2,mgp=c(1.75,0.5,0),mgp=c(2,0.5,0),cex.axis=0.95,mar=c(3.5,3.5,1.5,1))

myCol=c("black","cornflowerblue","red3")

matplot(sampleN,means.mu,xlab="Samples",ylab=expression(paste("Mean of log ",rho)),type="l",lty=1,col=myCol)
abline(h=mean(log(rawD$rho)),col="black",lty="dashed")

matplot(sampleN,means.sd,xlab="Samples",ylab=expression(paste("Standard deviation of log ",rho)),type="l",lty=1,col=myCol)
legend(x=0.5*max(sampleN),y=0.9*max(means.sd),legend=c(expression(paste(rho)),"Effects","Responses"),lty=1,col=myCol,bty="n",cex=0.9)

dev.off()

