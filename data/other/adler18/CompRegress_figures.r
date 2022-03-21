
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(viridis)
library(gridExtra)
library(grid)
library(gtable)
library(cowplot)
library(texreg)

###
### stats tables
###

statsOutput <- "stats_tables.doc"

# save model summary to file
htmlreg(list(lm.rho,effect.LR1,resp.LR1),statsOutput,ci.force = c(TRUE), ci.test = 0,ci.force.level = 0.95, bold = 0.05,
        override.ci.low = list(rho.CI[1],effect.CI[1],resp.CI[1]),
        override.ci.up = list(rho.CI[2],effect.CI[2],resp.CI[2]),
        caption="Summary of the mixed effect models",single.row=T,caption.above = T,
        custom.model.names=c("Rho","Effects","Responses"),append=F)

###
### build table showing studies and observations for each data set
###

tmpD1 <- as.data.frame(table(rhoD$Key))
names(tmpD1)[2]="Rho"

tmpD2 <- as.data.frame(table(effectD$Key))
names(tmpD2)[2]="Effects"

tmpD3 <- as.data.frame(table(responseS$Key))
names(tmpD3)[2]="Response"

tmpD4 <- as.data.frame(table(responseP$Key))
names(tmpD4)[2]="Pooled"

studyD <- merge(tmpD1,tmpD2)
studyD <- merge(studyD,tmpD3)
studyD <- merge(studyD,tmpD4)
names(studyD)[1] <- "Key"
studyD<-subset(studyD,Key!="BQPV6ZZF") # remove Lankau (never used, but stored in Key levels)
studyD<-subset(studyD,Key!="T4IP3JQW") # remove Dormann (never used, but stored in Key levels)

tmpD <- unique(cleanD[,c("Key","First.Author","Year","Ecosystem2")],MARGIN=2)

studyD <- merge(studyD,tmpD,all.x=T)
studyD <- studyD[,c(6:8,2:5)]
studyD <- studyD[order(studyD$First.Author,studyD$Year),]
write.csv(studyD,"StudyTable.csv",row.names=F)

rm(tmpD,tmpD1,tmpD2,tmpD3,tmpD4,studyD) # clean up

###
### Fig. S2:  pie charts Effects
###

myColors=c("cornflowerblue","violet","beige","gold","darkseagreen2")

png("EffectsPie.png",height=4,width=8,res=400,units="in")

par(mfrow=c(2,3),mar=c(0.5,0.5,1.5,0.5))

# Ecosystem type
pie(table(effectD$Ecosystem2),col=myColors,main="Ecosystem")

# setting
pie(table(effectD$lab.OR.field)[2:3],col=myColors,main="Setting") # all field

# observational vs experimental
pie(table(effectD$Densities)[3:4],col=myColors,main="Neighbor densities")

# response variables
pie(table(effectD$Response2),col=myColors,main="Fitness component")

# life stage
pie(table(effectD$Focal.Lifestage2),col=myColors,main="Life stage")  # not much variation, no seedlings only

# origin (we decided not to discuss native/non-native status in the manuscript)
# pie(table(effectD$Focal.Origin)[1:2],col=myColors,main="Origin")  # not much variation, no seedlings only

dev.off()

###
### Fig. S3: pie charts Response
###

png("ResponsePie.png",height=4,width=8,res=400,units="in")

par(mfrow=c(2,3),mar=c(0.5,0.5,1.5,0.5))

# Ecosystem type
pie(table(responseS$Ecosystem2),col=myColors,main="Ecosystem")

# setting
pie(table(responseS$lab.OR.field)[2:3],col=myColors,main="Setting") # all field

# observational vs experimental
pie(table(responseS$Densities)[3:4],col=myColors,main="Neighbor densities")

# response variables
pie(table(responseS$Response2),col=myColors,main="Fitness component")

# life stage
pie(table(effectD$Focal.Lifestage2),col=myColors,main="Life stage")  # not much variation, no seedlings only

# origin (we decided not to discuss native/non-native status in the manuscript)
# pie(table(responseS$Neighbor.Origin)[1:2],col=myColors,main="Origin")  # not much variation, no seedlings only

dev.off()

###
### Fig 1: pie charts rho dataset
###

png("RhoPie.png",height=4,width=8,res=400,units="in")

par(mfrow=c(2,3),mar=c(0.5,0.5,1.5,0.5))

# Ecosystem type
pie(table(rhoD$Ecosystem2),col=myColors,main="Ecosystem")

# setting
pie(table(rhoD$lab.OR.field)[2:3],col=myColors,main="Setting") # all field

# observational vs experimental
pie(table(rhoD$Densities)[3:4],col=myColors,main="Neighbor densities")

# response variables
pie(table(rhoD$Response2),col=myColors,main="Fitness component")

# life stage
pie(table(rhoD$Focal.Lifestage2),col=myColors,main="Life stage")  # not much variation, no seedlings only

# origin (we decided not to discuss native/non-native status in the manuscript)
# pie(table(responseS$Neighbor.Origin)[1:2],col=myColors,main="Origin")  # not much variation, no seedlings only

dev.off()

###
### Fig. 2: distribution rho, log scale
###
png("density-plots-log.png",height=2.75,width=7.5,res=400,units="in")
par(mfrow=c(1,3),tcl=-0.25,mgp=c(2,0.5,0),mar=c(3,1,2,1),oma=c(0,2.5,0,0),cex.lab=1.1)

#rho
d<-density(log(rhoD$rho),bw=1)
plot(d,xlab=expression(paste(log(rho))),ylab="",main="",ylim=c(0,0.24),yaxs="i")
polygon(d,col="slategray1")
abline(v=fixef(lm.rho)[1],col="red3",lwd=1.5)
abline(v=rho.CI[1],col="red3",lwd=0.5)
abline(v=rho.CI[2],col="red3",lwd=0.5)
abline(v=0,lty="dashed")
mtext(side=3,"(A)",line=0.5,adj=0)

#effects
d<-density(effectD$logratio[eff.subset],bw=1)
plot(d,xlab=expression(paste(log(alpha[ji]/alpha[ii]))),ylab="",main="",ylim=c(0,0.22),yaxs="i")
polygon(d,col="slategray1")
abline(v=fixef(effect.LR1)[1],col="red3",lwd=1.5)
abline(v=effect.CI[1],col="red3",lwd=0.5)
abline(v=effect.CI[2],col="red3",lwd=0.5)
abline(v=0,lty="dashed")
mtext(side=3,"(B)",line=0.5,adj=0)

#responses
d<-density(responseS$logratio[resp.subset],bw=1)
plot(d,xlab=expression(paste(log(alpha[ij]/alpha[ii]))),ylab="",main="",ylim=c(0,0.22),yaxs="i")
polygon(d,col="slategray1")
abline(v=fixef(resp.LR1)[1],col="red3",lwd=1.5)
abline(v=resp.CI[1],col="red3",lwd=0.5)
abline(v=resp.CI[2],col="red3",lwd=0.5)
abline(v=0,lty="dashed")
mtext(side=3,"(C)",line=0.5,adj=0)

mtext(side=2,"Density",outer=T,line=1,cex=0.9)

dev.off()


###
### distribution rho, one panel, color (not used in the paper)
###
png("density-plots-log-color.png",height=3.5,width=3.5,res=400,units="in")
par(mfrow=c(1,1),tcl=-0.25,mgp=c(2,0.5,0),mar=c(3,3,2,1),cex.lab=1.1)

#rho
d<-density(log(rhoD$rho),bw=1)
plot(d,xlab=expression(paste(log(rho))),ylab="Density",main="",ylim=c(0,0.24),yaxs="i")
polygon(d,col="wheat")
abline(v=fixef(lm.rho)[1],col="red3",lwd=1.5)
abline(v=rho.CI[1],col="red3",lwd=0.5)
abline(v=rho.CI[2],col="red3",lwd=0.5)
abline(v=0,lty="dashed")

dev.off()

###
### distribution rho, backtransformed (not used in the paper)
###

png("histograms-backtransformed.png",height=2.75,width=7.5,res=400,units="in")
par(mfrow=c(1,3),tcl=-0.25,mgp=c(2,0.5,0),mar=c(3,1,2,1),oma=c(0,2.5,0,0),cex.lab=1.1)

#rho
hist(rhoD$rho,xlab=expression(paste(rho)),breaks=20,ylab="",main="",yaxs="i",col="gray")
abline(v=exp(fixef(lm.rho)[1]),col="red3",lwd=1.5)
abline(v=exp(rho.CI[1]),col="red3",lwd=0.5)
abline(v=exp(rho.CI[2]),col="red3",lwd=0.5)
abline(v=1,lty="dashed")
mtext(side=3,"(A)",line=0.5,adj=0)

#effects
hist(exp(effectD$logratio[eff.subset]),breaks=20,xlab=expression(paste(alpha[ji]/alpha[ii])),ylab="",main="",yaxs="i",col="gray")
abline(v=exp(fixef(effect.LR1)[1]),col="red3",lwd=1.5)
abline(v=exp(effect.CI[1]),col="red3",lwd=0.5)
abline(v=exp(effect.CI[2]),col="red3",lwd=0.5)
abline(v=1,lty="dashed")
mtext(side=3,"(B)",line=0.5,adj=0)

#responses
hist(exp(responseS$logratio[resp.subset]),breaks=20,xlab=expression(paste(alpha[ji]/alpha[ii])),ylab="",main="",yaxs="i",col="gray")
abline(v=exp(fixef(resp.LR1)[1]),col="red3",lwd=1.5)
abline(v=exp(resp.CI[1]),col="red3",lwd=0.5)
abline(v=exp(resp.CI[2]),col="red3",lwd=0.5)
abline(v=1,lty="dashed")
mtext(side=3,"(C)",line=0.5,adj=0)

mtext(side=2,"Frequency",outer=T,line=1,cex=0.9)

dev.off()

###
### Fig. 3: plot covariate effects along with results of likelihood ratio tests
###

# function to get coefficients
getCoefs<- function(model){
  out <- fixef(model) + c(0,rep(fixef(model)[1],(length(fixef(model))-1)))
  return(out)
}

# custom barplot function
makeBarplot <- function(coefData,aov1,aov2,aov3,plot_title=NULL){
  
  #myCols<-terrain.colors(NROW(coefData))
  
  coefData<-gather(coefData,"Test","Coefficient",c(2:4))
  
  # next line is to force ggplot to order the bars the way want
  tmp <- length(unique(coefData$Class))
  coefData$Test <- as.factor(paste0(sort(rep(c("a","b","c"),tmp)),coefData$Test))
  
  tabD<-data.frame(Test=c("rho","Eff.","Resp."),
                   Chisq=round(c(aov1$Chisq[2],aov2$Chisq[2],aov3$Chisq[2]),2),
                   Pvals=round(c(aov1$'Pr(>Chisq)'[2],aov2$'Pr(>Chisq)'[2],aov3$'Pr(>Chisq)'[2]),3))
  colnames(tabD) <- c("Test","chi^2","italic(p)")
  
 
  mycex <- 0.6
  tt <- gridExtra::ttheme_minimal(
    core = list(fg_params=list(cex = mycex, parse=TRUE)),
    colhead = list(fg_params=list(cex = mycex, parse=TRUE)),
    rowhead = list(fg_params=list(cex = mycex)))
  table_element <- tableGrob(tabD, rows = NULL, theme=tt)
  table_element$widths <- unit(rep(0.8/ncol(table_element), ncol(table_element)), "npc")
  table_element$heights <- unit(rep(0.19/nrow(table_element), nrow(table_element)), "npc")
  table_element <- gtable_add_grob(table_element,
                                   grobs = segmentsGrob( # line across the bottom
                                     x0 = unit(0,"npc"),
                                     y0 = unit(0,"npc"),
                                     x1 = unit(1,"npc"),
                                     y1 = unit(0,"npc"),
                                     gp = gpar(lwd = 1.0)),
                                   t = 1, b = 1, l = 1, r = 3)
  
  gout <- ggplot(coefData,aes(x=Test,y=Coefficient,fill=Class)) + 
    geom_col(position = position_dodge(0.9), color="grey35") +
    scale_x_discrete(labels = c(expression(rho),"Eff.","Resp."), name = NULL) +
    scale_y_continuous(limits = c(-2.5,0.17))+
    scale_fill_viridis(end = 0.8, discrete = TRUE, name = NULL) +
    theme_classic() +
    theme(legend.justification = "top",
          legend.text=element_text(size=10),
          legend.key.size = unit(0.3, "cm"),
          plot.margin = unit(c(0.2,0,0.1,0), "cm"),
          plot.title = element_text(size = 10),
          axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black")) +
    annotation_custom(grob = table_element, xmin = 3.4, xmax = 6.3, ymin = -4) +
    coord_fixed(ratio=2/1)+
    ggtitle(plot_title)

  #add stats table

  #addtable2plot(max(myplot)*1.1,0.8*min(coefData),tabD,cex=0.8,bty="n",hlines=T,lwd=0.1)
  return(gout)

}

masterCols <- c("cornflowerblue","violet","beige","gold","darkseagreen2")

# compare ecosystems
tmpD <- data.frame(Class=sort(unique(effectD$Ecosystem2)),
                   Rho=getCoefs(rho.eco),
                   Effects=getCoefs(effect.LR.eco),
                   Responses=getCoefs(resp.LR.eco))
eco_plot <- makeBarplot(tmpD,aov.rho.eco,aov.eff.eco,aov.resp.eco,plot_title = "(E)")

# compare lab vs greenhouse
tmpD <- data.frame(Class=sort(unique(effectD$lab.OR.field)),
                   Rho=getCoefs(rho.gh),
                   Effects=getCoefs(effect.LR.gh),
                   Responses=getCoefs(resp.LR.gh))
gh_plot <- makeBarplot(tmpD,aov.rho.gh,aov.eff.gh,aov.resp.gh,plot_title = "(B)")

# compare observational vs experimental
tmpD <- data.frame(Class=sort(unique(effectD$Densities)),
                   Rho=getCoefs(rho.design),
                   Effects=getCoefs(effect.LR.design),
                   Responses=getCoefs(resp.LR.design))
design_plot <- makeBarplot(tmpD,aov.rho.design,aov.eff.design,aov.resp.design,plot_title = "(A)")

# compare fitness components
tmpD <- data.frame(Class=sort(unique(effectD$Response2)),
                   Rho=getCoefs(rho.fit),
                   Effects=getCoefs(effect.LR.fit),
                   Responses=getCoefs(resp.LR.fit))
tmpD[,"Rho"]<-tmpD[c(2,1,4,3),"Rho"] # order rho same as the others
fitness_plot <- makeBarplot(tmpD,aov.rho.fit,aov.eff.fit,aov.resp.fit,plot_title = "(C)")

# compare stage
tmpD <- data.frame(Class=sort(unique(effectD$Focal.Lifestage2)),
                   Rho=getCoefs(rho.stage),
                   Effects=getCoefs(effect.LR.stage),
                   Responses=getCoefs(resp.LR.stage))
tmpD[,"Rho"] <-tmpD[c(2,1,3),"Rho"] # order this one same as the others
stage_plot <- makeBarplot(tmpD,aov.rho.stage,aov.eff.stage,aov.resp.stage,plot_title = "(D)")

# Vertical
# covariates_plot <- plot_grid(eco_plot,gh_plot,design_plot,fitness_plot,stage_plot,ncol=2, nrow=3, align = "v")
# ggsave(filename = "ggplot_covariates.pdf", plot = covariates_plot, width = 6.5, height = 9, units = "in", dpi=120)

# Horizontal
covariates_plot <- plot_grid(design_plot,gh_plot,fitness_plot,stage_plot,eco_plot,ncol=3, nrow=2, align = "v")
ggsave(filename = "covariates.png", plot = covariates_plot, width = 8.5, height = 5.5, units = "in", dpi=400)


### Fig. S3: look at studies with pooled conspecific and heterospecific effects --------------------
responseP

png("pooled-boxplot.png",height=3,width=4,res=400,units="in")
par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,3,1,1))
boxplot(responseP[,c("alpha.conspp","alpha.heterospp")],ylab="Coefficient",
        names=c("Intraspecific","Interspecific"))
abline(h=0,lty="dotted",col="gray")
dev.off()
