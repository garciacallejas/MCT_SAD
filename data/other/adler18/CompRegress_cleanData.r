
# this script performs a number of misc. data clean-up operations
# on the "rawD" data frame and returns a "cleanD" data frame

# recode column names
names(rawD)[which(names(rawD)=="Paper.Key..")] <- "Key"
names(rawD)[which(names(rawD)=="Greenhouse.or.field.experiment")] <- "lab.OR.field"
names(rawD)[which(names(rawD)=="Manipulated.or.natural.densities.")] <- "Densities"
names(rawD)[which(names(rawD)=="Longterm...1.season.")] <- "Study.length"
names(rawD)

# fix drag and drop mistakes in spreadsheet for Key column
table(rawD$Key)
rawD$Key[grep("I8R6XFA",rawD$Key)]<-"I8R6XFA9"
rawD$Key[grep("XKG6JWW",rawD$Key)]<-"XKG6JWW4"
table(as.character(rawD$Key))

# check for odd categorical values
table(rawD$Response.variable)
# comments: diameter = growth, plant size = growth, Population growth = fitness
rawD$Response2 <-as.character(rawD$Response.variable)
rawD$Response2[grep("Population",rawD$Response2)] <- "Fitness"
rawD$Response2[grep("Fitness?",rawD$Response2)] <- "Fitness"
rawD$Response2[grep("Diameter",rawD$Response2)] <- "Growth"
rawD$Response2[grep("growth rate",rawD$Response2)] <- "Growth"
rawD$Response2[grep("Biomass",rawD$Response2)] <- "Growth"
rawD$Response2[grep("Height",rawD$Response2)] <- "Growth"
rawD$Response2[grep("Tiller count",rawD$Response2)] <- "Growth"
rawD$Response2[grep("Seed production",rawD$Response2)] <- "Fitness"
rawD$Response2[grep("Shoot biomass",rawD$Response2)] <- "Growth"
rawD$Response2[grep("Tiller count",rawD$Response2)] <- "Growth"
rawD$Response2[grep("Total biomass",rawD$Response2)] <- "Growth"
rawD$Response2[grep("Aboveground biomass",rawD$Response2)] <- "Growth"
# the "Plant size" response comes from Coomes 2002; but the response is actually seed production (fitness)
rawD$Response2[grep("Plant Size",rawD$Response2)] <- "Fitness"
rawD$Response2[grep("survival",rawD$Response2,ignore.case=TRUE)] <- "Survival"
table(rawD$Response2)

table(rawD$lab.OR.field)
# "Both" comes from Zarnetske; but the mesocosm data were used to describe intrinsic growth,
# not competitive effects, so I am changing this one to "field"
rawD$lab.OR.field[rawD$lab.OR.field=="Both"] <- "Field"
# Suter 2007 is coded as "field" but it is pots grown outdoors, left as is
table(rawD$lab.OR.field)

table(rawD$Ecosystem)
rawD$Ecosystem2 <- as.character(rawD$Ecosystem)
rawD$Ecosystem2[rawD$Ecosystem2=="Sub-boreal forest"] <- "Forest"
rawD$Ecosystem2[rawD$Ecosystem2=="Coastal Dune"] <- "Dune"
rawD$Ecosystem2[rawD$Ecosystem2=="Unknown"] <- "Agriculture"  # Dormann = "managed lawn"
rawD$Ecosystem2[rawD$Ecosystem2=="Annual plants"] <- "Grassland"  # Lankau = California annual grassland
rawD$Ecosystem2[rawD$Ecosystem2=="Grassland Perennials"] <- "Grassland" # Suter 2007 Swiss fen species
rawD$Ecosystem2[rawD$Ecosystem2=="Pasture"] <- "Agriculture"  
rawD$Ecosystem2[rawD$Ecosystem2=="Agricultural fields"] <- "Agriculture"
rawD$Ecosystem2[rawD$Ecosystem2=="Mediterranean Shrubland"] <- "Steppe" 
rawD$Ecosystem2[rawD$Ecosystem2=="Tropical forest"] <- "Forest"
rawD$Ecosystem2[rawD$Ecosystem2=="Sagebrush steppe"] <- "Steppe"
rawD$Ecosystem2[rawD$Ecosystem2=="Prairie"] <- "Grassland"
rawD$Ecosystem2[rawD$Ecosystem2=="Annual grassland"] <- "Grassland"
rawD$Ecosystem2[rawD$Ecosystem2=="Temperate forest"] <- "Forest"
rawD$Ecosystem2[rawD$Ecosystem2=="Serpentine grassland"] <- "Grassland"

table(rawD$Ecosystem2)

table(rawD$Study.length)
rawD$Study.length[rawD$Study.length=="No"] <- "Short term"
rawD$Study.length[rawD$Study.length=="yes"] <- "Longterm"
# "Both" comes from Zarnetske; but the mesocosm data were used to describe intrinsic growth,
# not competitive effects, so I am changing this one to "Longterm"
rawD$Study.length[rawD$Study.length=="Both"] <- "Longterm"
table(rawD$Study.length)

table(rawD$Densities)
rawD$Densities[rawD$Densities=="Both"] <- "Natural"  # see notes above on Zarnetske
rawD$Densities[rawD$Densities=="Experimental"] <- "Manipulated"
levels(rawD$Densities)[4] <- c("Observed") # rename this category
table(rawD$Densities)

table(rawD$Focal.Lifestage)
rawD$Focal.Lifestage2<-as.character(rawD$Focal.Lifestage)
# reclassify into "seedlings","mature" and "all"
rawD$Focal.Lifestage2[rawD$Focal.Lifestage2=="All sizes"] <- "All" 
rawD$Focal.Lifestage2[rawD$Focal.Lifestage2=="all sizes"] <- "All" 
rawD$Focal.Lifestage2[grep("seedling",rawD$Focal.Lifestage2)] <- "Seedling"
rawD$Focal.Lifestage2[which(rawD$Focal.Lifestage2!="Seedling" & rawD$Focal.Lifestage2!="All")] <- "Mature"
table(rawD$Focal.Lifestage2)

table(rawD$Comp.Lifestage)
rawD$Comp.Lifestage2<-as.character(rawD$Comp.Lifestage)
# reclassify into "seedlings","mature" and "all"
rawD$Comp.Lifestage2[rawD$Comp.Lifestage2=="All sizes"] <- "All" 
rawD$Comp.Lifestage2[rawD$Comp.Lifestage2=="all sizes"] <- "All" 
rawD$Comp.Lifestage2[grep("seedling",rawD$Comp.Lifestage2)] <- "Seedling"
rawD$Comp.Lifestage2[which(rawD$Comp.Lifestage2!="Seedling" & rawD$Comp.Lifestage2!="All")] <- "Mature"
table(rawD$Comp.Lifestage2)

table(rawD$Focal.Form)
# reclassify into "annual","perennial.herb" and "woody"
lookup<-data.frame(Focal.Form=names(table(rawD$Focal.Form)),
    Focal.Form2=c("annual","annual","annual","annual","perennial.herb","perennial.herb",
                  "perennial.herb","perennial.herb","perennial.herb","perennial.herb",
                  "perennial.herb","perennial.herb", "perennial.herb", "perennial.herb",
                  "perennial.herb","woody","woody","woody","woody","perennial.herb",
                  "perennial.herb"))
rawD<-merge(rawD,lookup,all.x=T)
table(rawD$Focal.Form2)

table(rawD$Comp.Form)
# reclassify into "annual","perennial.herb" and "woody"
lookup<-data.frame(Comp.Form=names(table(rawD$Comp.Form)),
    Comp.Form2=c("annual","annual","annual","annual","perennial.herb","perennial.herb",
                 "perennial.herb","perennial.herb","perennial.herb","perennial.herb",
                 "perennial.herb","perennial.herb", "perennial.herb", "perennial.herb",
                 "perennial.herb","woody","woody","woody","woody","woody","perennial.herb",
                 "perennial.herb"))
rawD<-merge(rawD,lookup,all.x=T)
table(rawD$Comp.Form2)
rm(lookup)

table(rawD$Focal.Origin)  
rawD$Focal.Origin[rawD$Focal.Origin=="Unknown"] <- "Native"  # PBA: I am pretty sure Suter 2007 is all natives

table(rawD$Comp.Origin)   # looks ok
rawD$Comp.Origin[rawD$Comp.Origin=="Unknown"] <- "Native"  # PBA: I am pretty sure Suter 2007 is all natives

cleanD<-rawD[,c("Key","First.Author","Year","Response2","lab.OR.field","Ecosystem2",
  "Study.length","Densities","Focal.Species","Focal.Form2","Focal.Lifestage2","Focal.Origin",
  "Comp.Species","Comp.Form2","Comp.Lifestage2","Comp.Origin","treatment","nonsig.zero",
  "competition.coefficient","Negative.coefs.mean.competition","standard.error.of.coefficient",    
  "Description.of.Uncertainty.Metric", "P.Value","DF","CI.lower","CI.upper","Comments"   )]



