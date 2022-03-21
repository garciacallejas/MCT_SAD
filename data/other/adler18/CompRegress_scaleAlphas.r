
# This script is called from CompRegress_wrapper.r.
# It assumes the cleanD data.frame is loaded.

# see where the problems are
table(cleanD$Negative.coefs.mean.competition)

# remove Lankau 2009--it's a relative, not per capita, measure of intra vs interspecific effects
tmp <- which(cleanD$First.Author=="Lankau" & cleanD$Year==2009)
cleanD <- cleanD[-tmp,]

# remove alpha = 0 values from Farrer paper (non-significant effects coded as zero)
tmp<-which(cleanD$First.Author=="Farrer" & cleanD$competition.coefficient==0)
cleanD <- cleanD[-tmp,]

# rescale Turnbull 2004 values. These are exponential decay parameters--they give
# the impact of a neighbor at a certain distance. We will rescale using:
#     exp(-alpha*5)
# which is the impact of a plant at 5 cm from the target. A value of zero
# means no competitive effect, more positive values mean more competition.
# I picked 5 cm after looking at the parameter estimates and the corresponding curves,
# 5 cm is often near the maximum difference between curves (difference goes to zero at large distance)
tmp <- which(cleanD$First.Author=="Turnbull" & cleanD$Year==2004)
cleanD$Negative.coefs.mean.competition[tmp] <- "No"
cleanD$competition.coefficient[tmp] <- exp(cleanD$competition.coefficient[tmp]*5)

# make sure competition always takes negative values
table(cleanD$Negative.coefs.mean.competition)
cleanD$Negative.coefs.mean.competition[cleanD$Negative.coefs.mean.competition=="yes"]<-"Yes"
cleanD$Negative.coefs.mean.competition[cleanD$Negative.coefs.mean.competition=="no"]<-"No"
table(cleanD$Negative.coefs.mean.competition)
tmp <- which(cleanD$Negative.coefs.mean.competition=="No")
cleanD$competition.coefficient[tmp] <- -1*cleanD$competition.coefficient[tmp]

# split intra and interspecific effects
cleanD$Focal.Species <- as.character(cleanD$Focal.Species)
cleanD$Comp.Species <- as.character(cleanD$Comp.Species)
cleanD$comp.type <- ifelse(cleanD$Focal.Species==cleanD$Comp.Species,"Intra","Inter")
cbind(cleanD$Focal.Species,cleanD$Comp.Species,cleanD$comp.type) # make sure this worked

# join intra and inter based on effects, not responses-----------------------

intraD<-subset(cleanD,comp.type=="Intra")
names(intraD)[which(names(intraD)=="competition.coefficient")]<-"alpha.intra"
names(intraD)[which(names(intraD)=="Focal.Species")]<-"Intra.species"
names(intraD)[which(names(intraD)=="Comp.Species")]<-"Neighbor.species"
intraD<-intraD[,c("Key","First.Author","Year","Response2","Focal.Lifestage2","Focal.Origin","Intra.species","Neighbor.species","treatment","alpha.intra")]

interD<-subset(cleanD,comp.type=="Inter")
names(interD)[which(names(interD)=="competition.coefficient")]<-"alpha.inter"
names(interD)[which(names(interD)=="Focal.Species")]<-"Target.species"
names(interD)[which(names(interD)=="Comp.Species")]<-"Neighbor.species"
names(interD)[which(names(interD)=="Comp.Origin")]<-"Neighbor.Origin"
interD<-interD[,c("Key","First.Author","Year","Response2","Target.species","Neighbor.species","Neighbor.Origin","treatment","alpha.inter")]

effectD<-merge(intraD,interD,all.y=T)
effectD <- subset(effectD,!is.na(effectD$alpha.intra))  # remove NAs

effectD$logratio <- log(effectD$alpha.inter/effectD$alpha.intra)

# join intra and inter based on responses, not effects--------------------------

intraD<-subset(cleanD,comp.type=="Intra")
names(intraD)[which(names(intraD)=="competition.coefficient")]<-"alpha.intra"
names(intraD)[which(names(intraD)=="Focal.Species")]<-"Target.species"
names(intraD)[which(names(intraD)=="Comp.Species")]<-"Neighbor.species.intra"
intraD<-intraD[,c("Key","First.Author","Year","Response2","Focal.Lifestage2","Focal.Origin","Target.species","Neighbor.species.intra","treatment","alpha.intra")]

interD<-subset(cleanD,comp.type=="Inter")
names(interD)[which(names(interD)=="competition.coefficient")]<-"alpha.inter"
names(interD)[which(names(interD)=="Focal.Species")]<-"Target.species"
names(interD)[which(names(interD)=="Comp.Species")]<-"Neighbor.species.inter"
names(interD)[which(names(interD)=="Comp.Origin")]<-"Neighbor.Origin"
interD<-interD[,c("Key","First.Author","Year","Response2","Target.species","Neighbor.species.inter","Neighbor.Origin","treatment","alpha.inter")]

responseD<-merge(intraD,interD,all.y=T)

# put studies which pool species responses into their own data frame

tmp <- grep("All",responseD$Target.species)
responseP <- responseD[tmp,]
responseP$Neighbor.species.inter[responseP$Neighbor.species.inter=="conspecifics"] <- "Conspecifics"
responseC <- subset(responseP,Neighbor.species.inter=="Conspecifics")
responseH <- subset(responseP,Neighbor.species.inter=="Heterospecifics")
rm(responseP)
responseC <- responseC[,c("Key","First.Author","Year","Response2","treatment","alpha.inter")]
names(responseC)[names(responseC)=="alpha.inter"] <- "alpha.conspp"
responseH <- responseH[,c("Key","First.Author","Year","Response2","treatment","alpha.inter")]
names(responseH)[names(responseH)=="alpha.inter"] <- "alpha.heterospp"
responseP <- merge(responseC,responseH)
rm(responseC,responseH)

# remove pooled species from species-specific data set
responseS <- responseD[-tmp,]
responseS <- subset(responseS,!is.na(alpha.intra))  # remove NAs
rm(responseD)

# flag "composite" heterospp records
tmp <- grep("all",responseS$Neighbor.species.inter)
tmp <- c(tmp,grep("other",responseS$Neighbor.species.inter,ignore.case=T))
responseS$composite.neighbor <- 0
responseS$composite.neighbor[tmp]<-1

responseS$logratio <- log(responseS$alpha.inter/responseS$alpha.intra)

  