
# This script is called from CompRegress_wrapper.r.
# It assumes the effectD and responseP, responseS data.frames are loaded.

# library(lme4)  # loaded by CompRegress_wrapper.R

covars <- cleanD[,c("Key","lab.OR.field","Ecosystem2","Study.length","Densities")]
covars <- unique(covars,MARGIN=2)
effectD <- merge(effectD,covars,all.x=T)
responseS <- merge(responseS,covars,all.x=T)


# analyze EFFECT data --------------------------------------------------------------

effectD$Source<-effectD$Key

# first, count different combinations of facilitation (F) and competition (C)
FC <- which(effectD$alpha.intra>=0 & effectD$alpha.inter <0) 
CF <- which(effectD$alpha.intra<0 & effectD$alpha.inter >=0)
FF <- which(effectD$alpha.intra>=0 & effectD$alpha.inter >= 0)
CC <- which(effectD$alpha.intra<0 & effectD$alpha.inter<0)
CFcounts<-c(length(FC),length(CF),length(FF),length(CC))
names(CFcounts) <- c("FC","CF","FF","CC")
print(CFcounts)

# look at facilitation records
effectD[FC,]  # all Greenhouse studies

# compare intras and inters with negative effects
eff.subset <- c(CC)

# calculate proportion of these cases in which intra > inter
sum(effectD$logratio[eff.subset]<0)/length(eff.subset)

###
### simplest model
###
effect.LR1 <- lmer(logratio~1+(1|Source/Intra.species),data=effectD,subset=eff.subset)
# plot(predict(effect.LR1),residuals(effect.LR1)) # looks ok, but some heteroscedasticity
effect.CI=confint(effect.LR1,parm="(Intercept)")
# look at back-transformed values
print(exp(fixef(effect.LR1)[1]))
print(exp(effect.CI))

# repeat simple model without the 66 records from Zarnetske 2013 -> little change
tmp=which(effectD$alpha.intra<0 & effectD$alpha.inter<0 & effectD$Source!="Zarnetske2013")
effect.LR1b <- lmer(logratio~1+(1|Source/Intra.species),data=effectD,subset=tmp)
myCIb=confint(effect.LR1b,parm="(Intercept)")

# repeat simple model without the 67 records from Geijzendorffer2011 -> little change
tmp=which(effectD$alpha.intra<0 & effectD$alpha.inter<0 & effectD$Source!="Geijzendorffer2011")
effect.LR1b <- lmer(logratio~1+(1|Source/Intra.species),data=effectD,subset=tmp)
myCIb=confint(effect.LR1b,parm="(Intercept)")

###
### consider additional covariates
###

# test greenhouse vs field
effect.LR.gh <- lmer(logratio~lab.OR.field+(1|Source/Intra.species),data=effectD,subset=eff.subset)
summary(effect.LR.gh)  
aov.eff.gh <- anova(effect.LR1,effect.LR.gh) # likelihood ratio test not-significant

# test different ecosystem types
effect.LR.eco <- lmer(logratio~Ecosystem2+(1|Source/Intra.species),data=effectD,subset=eff.subset)
summary(effect.LR.eco)  
aov.eff.eco <- anova(effect.LR1,effect.LR.eco) # likelihood ratio test not-significant

# test experimental vs observational 
effect.LR.design <- lmer(logratio~Densities+(1|Source/Intra.species),data=effectD,subset=eff.subset)
summary(effect.LR.design)  # marginal effect?
aov.eff.design <- anova(effect.LR1,effect.LR.design) # likelihood ratio test marginal

# test fitness component
effect.LR.fit <- lmer(logratio~Response2+(1|Source/Intra.species),data=effectD,subset=eff.subset)
summary(effect.LR.fit)  # sig. effect
aov.eff.fit <- anova(effect.LR1,effect.LR.fit) # likelihood ratio test not-significant

# test lifestage
effect.LR.stage <- lmer(logratio~Focal.Lifestage2 +(1|Source/Intra.species),data=effectD,subset=eff.subset)
summary(effect.LR.stage)  # no effect
aov.eff.stage <- anova(effect.LR1,effect.LR.stage) # likelihood ratio test not-significant

# test native vs exotic
# Note: We removed any mention of this test from the manuscript, in the interest
# of simplicity, but we kept it here for the sake of completenes
# effectD$Origin <- as.factor(paste0(effectD$Focal.Origin,effectD$Neighbor.Origin))
# effect.LR.origin <- lmer(logratio~Origin +(1|Source/Intra.species),data=effectD,subset=eff.subset)
# summary(effect.LR.origin)  # no effect
# aov.eff.stage <- anova(effect.LR1,effect.LR.origin) # likelihood ratio test not-significant


### just for fun: 

# take a peek at intras and inters with different signs
cbind(effectD$Source[CF],effectD$alpha.intra[CF],effectD$alpha.inter[CF])

# calculate  ratio of absolute magnitudes (facilitative inter / competitive intra)
absLR <- effectD$alpha.inter[CF] / (-1* effectD$alpha.intra[CF] )
sum(absLR==0)
hist(absLR[absLR>0],xlab="ratio",main="")
abline(v=mean(absLR[absLR>0]),col="red")
print(mean(absLR[absLR>0]))  
print(median(absLR[absLR>0]))

# look at facilitation intras and inters
cbind(effectD$Source[FF],effectD$alpha.intra[FF],effectD$alpha.inter[FF])
print(c(effectD$alpha.inter[FF]/effectD$alpha.intra[FF]))


# analyze RESPONSE data --------------------------------------------------------------

responseS$Source <- responseS$Key

# first, count different combinations of facilitation (F) and competition (C)
FC <- which(responseS$alpha.intra>=0 & responseS$alpha.inter <0) #none!
CF <- which(responseS$alpha.intra<0 & responseS$alpha.inter >=0)
FF <- which(responseS$alpha.intra>=0 & responseS$alpha.inter >= 0)
CC <- which(responseS$alpha.intra<0 & responseS$alpha.inter<0)
CFcounts<-c(length(FC),length(CF),length(FF),length(CC))
names(CFcounts) <- c("FC","CF","FF","CC")
print(CFcounts)

resp.subset <- CC

resp.LR1 <- lmer(logratio~1+(1|Source/Target.species),data=responseS,subset=resp.subset)
summary(resp.LR1)  # log ratio significantly negative
plot(predict(resp.LR1),residuals(resp.LR1)) # looks ok, but some heteroscedasticity
resp.CI=confint(resp.LR1,parm="(Intercept)")
print(resp.CI)
# look at back-transformed values
print(exp(fixef(resp.LR1)[1]))
print(exp(resp.CI))

# test greenhouse vs field
resp.LR.gh <- lmer(logratio~lab.OR.field+(1|Source/Target.species),data=responseS,subset=resp.subset)
summary(resp.LR.gh)  
aov.resp.gh <- anova(resp.LR1,resp.LR.gh) # likelihood ratio marginally significant

# test different ecosystem types
resp.LR.eco <- lmer(logratio~Ecosystem2+(1|Source/Target.species),data=responseS,subset=resp.subset)
summary(resp.LR.eco)  
aov.resp.eco<- anova(resp.LR1,resp.LR.eco) # likelihood ratio test not-significant

# test experimental vs observational 
resp.LR.design <- lmer(logratio~Densities+(1|Source/Target.species),data=responseS,subset=resp.subset)
summary(resp.LR.design)  # marginal effect?
aov.resp.design <- anova(resp.LR1,resp.LR.design) # likelihood ratio test marginal

# test fitness component
resp.LR.fit <- lmer(logratio~Response2+(1|Source/Target.species),data=responseS,subset=resp.subset)
summary(resp.LR.fit)  # sig. effect
aov.resp.fit <- anova(resp.LR1,resp.LR.fit) # likelihood ratio test significant

# test lifestage
resp.LR.stage <- lmer(logratio~Focal.Lifestage2 +(1|Source/Target.species),data=responseS,subset=resp.subset)
summary(resp.LR.stage)  # no effect
aov.resp.stage <- anova(resp.LR1,resp.LR.stage) # likelihood ratio test not-significant

# test native vs exotic
# Note: We removed any mention of this test from the manuscript, in the interest
# of simplicity, but we kept it here for the sake of completeness
# responseS$Origin <- as.factor(paste0(responseS$Focal.Origin,responseS$Neighbor.Origin))
# resp.LR.origin <- lmer(logratio~Origin +(1|Source/Target.species),data=responseS,subset=resp.subset)
# summary(resp.LR.origin)  # no effect
# aov.resp.origin <- anova(resp.LR1,resp.LR.origin) # likelihood ratio test not-significant


## just for fun:

# take a peek at intras and inters with different signs
cbind(responseS$Source[CF],responseS$alpha.intra[CF],responseS$alpha.inter[CF])

# calculate  ratio of absolute magnitudes (facilitative inter / competitive intra)
absLR <- responseS$alpha.inter[CF] / (-1* responseS$alpha.intra[CF] )
sum(absLR==0)
hist(absLR[absLR>0],xlab="ratio",main="")
abline(v=mean(absLR[absLR>0]),col="red")
print(mean(absLR[absLR>0]))  
print(median(absLR[absLR>0]))

# take a peek at intra (-) and inter (+) with different signs
cbind(responseS$Source[FC],responseS$alpha.intra[FC],responseS$alpha.inter[FC])

# look at facilitation intras and inters
cbind(responseS$Source[FF],responseS$alpha.intra[FF],responseS$alpha.inter[FF])
print(c(responseS$alpha.inter[FF]/responseS$alpha.intra[FF]))


