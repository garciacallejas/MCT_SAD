
# filter the raw competition/abundances data to include in the parameter estimation

library(tidyverse)

abund <- read.csv2("../Caracoles/data/abundances.csv",
                   header = TRUE,stringsAsFactors = FALSE)
comp <- read.csv2("../Caracoles/data/competition_wide.csv",
                  header = TRUE,stringsAsFactors = FALSE)
sp.rates <- read.csv2("../Caracoles/data/plant_species_traits.csv",
                      header = TRUE,stringsAsFactors = FALSE)
sp.valid <- sp.rates$species.code[which(!is.na(sp.rates$germination.rate))]

base.abund <- abund %>% 
  filter(species %in% sp.valid) 
# group_by(species) %>% 
# summarise(num = sum(individuals))
# year.off <- 2019

comp <- comp[,which(names(comp) %in% c("year","plot","subplot",
                                       "focal","fruit","seed",sp.valid))]
comp <- subset(comp, seed > 0)

all.sp <- names(comp)[7:length(names(comp))]
focal.sp <- unique(comp$focal[which(comp$focal %in% all.sp)])


# -------------------------------------------------------------------------
# filters:
# 1 - leave 2018 out

base.abund.1 <- base.abund %>% filter(year != 2018)
comp.1 <- comp %>% filter(year != 2018)

# 2 - filter by abundance? e.g. leave out very rare species, or those that don't appear all years

abund.table <- base.abund.1 %>%
  group_by(species, year) %>%
  summarise(min.abund = min(individuals),
            max.abund = max(individuals),
            total = sum(individuals))

# e.g. LYTR, CHMI have <100 focal observations
# SUSP does not have focal observations
# table(comp.1$focal) 

to.remove <- c("LYTR","CHMI","SUSP")
base.abund.2 <- subset(base.abund.1,!species %in% to.remove)
comp.2 <- subset(comp.1,!focal %in% to.remove)
comp.2$LYTR <- NULL
comp.2$CHMI <- NULL
comp.2$SUSP <- NULL

# 3 - combine melilotus sp. into one
mesp.abund <- subset(base.abund.2, species %in% c("MEEL","MESU"))
mesp.comp <- subset(comp.2,focal %in% c("MEEL","MESU"))
others.comp <- subset(comp.2, !focal %in% c("MEEL","MESU"))

mesp.abund.2 <- mesp.abund %>% group_by(year,plot,subplot) %>% summarise(species = "MESP",sum.ind = sum(individuals))
names(mesp.abund.2)[5] <- "individuals"

others.abund <- subset(base.abund.2[,c("year","plot","subplot","species","individuals")], 
                       !species %in% c("MEEL","MESU"))
base.abund.3 <- rbind(others.abund,mesp.abund.2)

mesp.comp$focal <- "MESP"
mesp.comp$MESP <- mesp.comp$MEEL + mesp.comp$MESU
mesp.comp$MEEL <- NULL
mesp.comp$MESU <- NULL

others.comp$MESP <- others.comp$MEEL + others.comp$MESU
others.comp$MEEL <- NULL
others.comp$MESU <- NULL

comp.3 <- rbind(others.comp,mesp.comp)
comp.3 <- comp.3[,c("year","plot","subplot","focal","fruit","seed","BEMA",
                    "CETE","CHFU","HOMA","LEMA","MESP","PAIN","PLCO","POMA",
                    "POMO","PUPA","SASO","SCLA","SOAS","SPRU")]

write.csv2(comp.3,file = "results/neigh_filtered.csv",row.names = FALSE)
write.csv2(base.abund.3,file = "results/abund_filtered.csv",row.names = FALSE)






