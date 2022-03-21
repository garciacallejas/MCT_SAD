

setwd("C:/Users/adler/Box Sync/activeProjects/CompetitionReview/datapluscode")

#-----------------------------------------------------------------------------
## Godoy et al 2014 Ecol Letts
# data from: https://figshare.com/articles/Dataset_Godoy_et_al_2014_Ecology_Letters/4793488

# Note that the data file below is not included in data package, 
# download from url above and save. 
datafile <- "../post2014/godoy2014/alpha_estimates_row_is_target.csv"  # change path as needed
D <- read.csv(datafile,row.names = 1)
spp_names <- row.names(D)
combos <- outer(spp_names,spp_names,FUN="paste0")
combos <- combos[upper.tri(combos,diag=F)]
spp_i <- substring(combos,1,4)
spp_j <- substring(combos,5,8)
rhoD_godoy <- data.frame(spp_i,spp_j,rho=NA)
for(row_i in 1:NROW(rhoD_godoy)){
  i <- which(spp_names==spp_i[row_i])
  j <- which(spp_names==spp_j[row_i])
  intras <- c(D[i,i],D[j,j])
  inters <- c(D[i,j],D[j,i])
  rhoD_godoy$rho[row_i] <- sqrt((inters[1]*inters[2])/(intras[1]*intras[2]))
}
keepers <- which(!is.na(rhoD_godoy$rho))
log_rho_mean <- mean(log(rhoD_godoy$rho[keepers]))
log_rho_sd <- sd(log(rhoD_godoy$rho[keepers]))
rho_mean <- exp(log_rho_mean)
print(length(keepers))
print(log_rho_mean)
print(rho_mean)

# Now re-do overall rho analysis with these data.
# to do this, first run CompRegress_wrapper up at least through the line:
# 'source("CompRegress_analyzeRho.r")'
# Then run the chunk above to read in Godoy and calculate rho

# combine data
rhoD2 <- data.frame(rhoD[,c("Source","Densities","rho")],stringsAsFactors = F)
tmp <-  data.frame(Source=rep("Godoy2014",length(keepers)),
                    Densities=rep("Manipulated",length(keepers)),
                    rho=rhoD_gody$rho[keepers],stringsAsFactors = F)
rhoD2 <- rbind(rhoD2,tmp)

# rerun overall model
library(lme4)
lm.rho.addGody <- lmer(log(rho)~1 + (1|Source),data=rhoD2)
print(summary(lm.rho.addGody))
rho.CI=confint(lm.rho.addGody,parm="(Intercept)")
print(rho.CI)
# look at back-transformed values
print(exp(fixef(lm.rho.addGody)[1]))
print(exp(rho.CI))

# rerun "Design" model
rho.design.addGodoy <- lmer(log(rho)~Densities + (1|Source),data=rhoD2)
summary(rho.design.addGodoy)  
aov.rho.design.addGodoy <- anova(lm.rho.addGody,rho.design.addGodoy) # not significant
# LL ratio test p value drops from 0.046 to 0.025

#-----------------------------------------------------------------------------
## Lamanna et al. 2017 Science, pooled estimates of CNDD and HNDD
# Data from Table S3: http://science.sciencemag.org/content/sci/suppl/2017/06/28/356.6345.1389.DC1/aam5678_LaManna_SM.pdf
# use means, focus on adults (not heterospecific sapling effects)

# Note that the data file below is not included in data package, 
# download from url above, convert Table S3 to csv format, and save.
datafile <- "../post2014/lamanna_S3/lamanna_S3_sheet.csv" # adjust path as needed
D <- read.csv(datafile)

CNDD <- mean(D$Mean_CNDD)
CNDD_sd <-  sd(D$Mean_CNDD)

HNDD <- mean(D$Mean_adult_HNDD)
HNDD_sd <- sd(D$Mean_adult_HNDD)

print(NROW(D))
print(c(CNDD,CNDD_sd,HNDD,HNDD_sd))
      
boxplot(D[,c("Mean_CNDD","Mean_adult_HNDD")],ylab="Coefficient")
abline(h=0,lty="dotted")


## Usinowiz et al 2017
# rho = sqrt(AijAji/1), since intraspecific is set equal to 1
# eyeballing:
low=sqrt(0.35)
hi=sqrt(0.52)

