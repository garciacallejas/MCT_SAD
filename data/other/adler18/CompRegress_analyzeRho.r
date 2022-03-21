
# This script is called from CompRegress_wrapper.r.
# It assumes the interD and intraD dataframes are loaded

# clean up treatment column to prevent multiple matches
interD$treatment=as.character(interD$treatment)
intraD$treatment=as.character(intraD$treatment)
interD$treatment[is.na(interD$treatment)]<-"Control"
intraD$treatment[is.na(intraD$treatment)]<-"Control"

# for Adler 2010 paste fitness component in as treatment
tmp=which(interD$First.Author=="Adler" & interD$Year=="2010")
interD$treatment[tmp]<- interD$Response2[tmp]
tmp=which(intraD$First.Author=="Adler" & intraD$Year=="2010")
intraD$treatment[tmp]<-intraD$Response2[tmp]

# put together full 2x2 pairwise competition matrices
rhoD <- data.frame(Key=character(),treatment=character(),Response2=character(),Focal.Lifestage2=character(),
          spp1=character(),spp2=character(),spp1.origin=character(),spp2.origin= character(),
          a.11=numeric(),a.22=numeric(),a.12=numeric(),a.21=numeric(),stringsAsFactors = T)
rho.inters<-numeric()
studyList <- sort(unique(intraD$Key))
# studyList<-studyList[-which(studyList=="KTP2IQBR")] # remove Zarnetske
for(iStudy in 1:length(studyList)){
  tmpD1<-subset(intraD,Key==studyList[iStudy])
  trtList<-unique(tmpD1$treatment)
  for(iTrt in 1:length(trtList)){
    tmpD2 <-subset(tmpD1,treatment==trtList[iTrt])
    sppList<-unique(tmpD2$Target.species)
    if(length(sppList)>1){
      for(iSpp1 in 1:(length(sppList)-1)){
        for(iSpp2 in (iSpp1+1):(length(sppList))){
          spp1<-sppList[iSpp1]; spp2<-sppList[iSpp2]
          intraIndex <- which(intraD$Key==studyList[iStudy] &
                    intraD$treatment==trtList[iTrt] & 
                   (intraD$Target.species==spp1 | intraD$Target.species==spp2))
          interIndex <- which(interD$Key==studyList[iStudy] & 
                  interD$treatment==trtList[iTrt] & 
                  (interD$Target.species==spp1 | interD$Target.species==spp2) &
                  (interD$Neighbor.species.inter==spp1 | interD$Neighbor.species.inter==spp2))
          if(length(intraIndex)>2 | length(interIndex)>2) stop("Too many alphas") # this shouldn't happen, investigate if it does
          if(length(intraIndex)<2 | length(interIndex)<2) next  # incomplete alpha matrix, exit loop
          #grab the alphas
          a.11=intraD$alpha.intra[intraIndex[1]]; a.22=intraD$alpha.intra[intraIndex[2]]; 
          a.12=interD$alpha.inter[interIndex[1]]; a.21=interD$alpha.inter[interIndex[2]]
          # if all negative, add to data frame
          if(sum(c(a.11,a.22,a.12,a.21)<rep(0,4))==4){
            rho.inters=c(rho.inters,interIndex)
            tmp=data.frame(Key=studyList[iStudy],treatment=trtList[iTrt],
                           Response2=intraD$Response2[intraIndex[1]],Focal.Lifestage2=intraD$Focal.Lifestage2[intraIndex[1]],
                           spp1.origin=intraD$Focal.Origin[intraIndex[1]],spp2.origin=intraD$Focal.Origin[intraIndex[1]],
                           spp1,spp2,a.11,a.22,a.12,a.21)
            rhoD=rbind(rhoD,tmp)
          } # if all alphas negative
          
        } # next spp2 
      } # next spp1
    } # if > 2 intra spp observed
  } # next treatment within study
} # next study

print(NROW(rhoD))
table(rhoD$Key)
rhoD$rho<-sqrt((rhoD$a.12*rhoD$a.21)/(rhoD$a.11*rhoD$a.22))

# ANALYZE RESULTS -----------------------------------

rhoD <- merge(rhoD,covars,all.x=T)

rhoD$Source <- rhoD$Key

sort(table(rhoD$Key))

# library(lme4) # this is loaded in CompRegress_wrapper.R

# simple model
lm.rho <- lmer(log(rho)~1 + (1|Source),data=rhoD)
# plot(predict(lm.rho),residuals(lm.rho)) # looks good
print(summary(lm.rho))
rho.CI=confint(lm.rho,parm="(Intercept)")
print(rho.CI)
# look at back-transformed values
print(exp(fixef(lm.rho)[1]))
print(exp(rho.CI))

# test greenhouse vs field
rho.gh <- lmer(log(rho)~lab.OR.field + (1|Source),data=rhoD)
summary(rho.gh)  
aov.rho.gh <- anova(lm.rho,rho.gh) # not significant

# ecosystem type
rho.eco <- lmer(log(rho)~Ecosystem2 + (1|Source),data=rhoD)
summary(rho.eco)  
aov.rho.eco <- anova(lm.rho,rho.eco) # not significant

# Experimental vs Observational
rho.design <- lmer(log(rho)~Densities + (1|Source),data=rhoD)
summary(rho.design)  
aov.rho.design <- anova(lm.rho,rho.design) # not significant

# Fitness components
rho.fit <- lmer(log(rho)~Response2 + (1|Source),data=rhoD)
summary(rho.fit)  
aov.rho.fit <- anova(lm.rho,rho.fit) # marginally significant

# Focal lifestage
rho.stage <- lmer(log(rho)~Focal.Lifestage2 + (1|Source),data=rhoD)
summary(rho.stage)  
aov.rho.stage <- anova(lm.rho,rho.stage) 

# given sample sizes, no point to analyze native vs non-native
#table(paste0(rhoD$spp1.origin,rhoD$spp2.origin))



