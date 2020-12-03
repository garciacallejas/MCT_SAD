# obtain the three components of RADs 
library(tidyverse)
source("R/hill_diversity.R")

abund.proj <- read.csv2("./results/projected_abundances_subplot.csv",header = TRUE,stringsAsFactors = FALSE)
abund.obs <- read.csv2("../Caracoles/data/abundances.csv",header = TRUE,stringsAsFactors = FALSE)

# average over plots
# plot-level values are averages over subplots or sums of all subplots?

plot.obs <- abund.obs %>%
  filter(!is.na(individuals) &
           individuals > 0) %>%
  group_by(year,plot,species) %>%
  summarise(mean.abund = mean(individuals),sum.abund = sum(individuals))
plot.obs$SAD.model <- "observed"
plot.obs <- plot.obs[,c("year","plot","SAD.model","species","mean.abund","sum.abund")]

plot.proj <- abund.proj %>%
  filter(!is.na(individuals)) %>%
  group_by(year,plot,SAD.model,sp) %>%
  summarise(mean.abund = mean(individuals),sum.abund = sum(individuals))
names(plot.proj)[4] <- "species"
plot.proj <- plot.proj[,c("year","plot","SAD.model","species","mean.abund","sum.abund")]

plot.all <- rbind(plot.obs,plot.proj)

# extract abundance, richness, evenness -----------------------------------
years <- sort(unique(plot.all$year))
plots <- sort(unique(plot.all$plot))
models <- sort(unique(plot.all$SAD.model))
components <- c("abundance","richness","evenness")

rad.components <- expand.grid(year = years,
                              plot = plots,
                              SAD.model = models,richness = NA_real_,abundance = NA_real_,evenness = NA_real_)

for(i.year in 1:length(years)){
  for(i.plot in 1:length(plots)){
    for(i.model in 1:length(models)){
      my.rad <- subset(plot.all,year == years[i.year] & 
                         plot == plots[i.plot] & 
                         SAD.model == models[i.model])
      
      pos <- which(rad.components$year == years[i.year] & rad.components$plot == plots[i.plot] & 
                     rad.components$SAD.model == models[i.model])
      
      rad.components$abundance[pos] <- sum(my.rad$sum.abund)
      rad.components$richness[pos] <- nrow(my.rad)
      rad.components$evenness[pos] <- hill.diversity(my.rad$sum.abund)
    }
  }
}

rad.components.long <- gather(rad.components, key = "component", value = "value", abundance,richness,evenness)

write.csv2(rad.components.long,"./results/rad_components.csv",row.names = FALSE)

# obtain RAD for plotting -------------------------------------------------

# 1 - observed
# rank.obs <- abund.obs %>%
#   filter(!is.na(individuals) &
#            individuals > 0) %>%
#   group_by(year,plot,species) %>%
#   
#   # summarise(mean.abund = mean(individuals)) %>%
#   summarise(mean.abund = sum(individuals)) %>%
#   
#   group_by(year,plot) %>%
#   mutate(plot.abund = sum(mean.abund)) %>%
#   group_by(year,plot,species) %>%
#   summarise(rel.abund = mean.abund/plot.abund) %>%
#   mutate(sp.rank = rank(-rel.abund,ties.method = "first"))
# 
# rank.obs <- arrange(rank.obs,year,plot,desc(rel.abund))
# rank.obs$SAD.model <- "observed"
# names(rank.obs)[which(names(rank.obs) == "species")] <- "sp"
# 
# # 2 - projected
# rank.abundances <- abund.proj %>%
#   filter(!is.na(individuals)) %>%
#   group_by(year,plot,SAD.model,sp) %>%
#   
#   # summarise(mean.abund = mean(individuals)) %>%
#   summarise(mean.abund = sum(individuals)) %>%
#   
#   group_by(year,plot,SAD.model) %>%
#   mutate(plot.abund = sum(mean.abund)) %>%
#   group_by(year,plot,SAD.model,sp) %>%
#   summarise(rel.abund = mean.abund/plot.abund) %>%
#   mutate(sp.rank = rank(-rel.abund,ties.method = "first"))
# 
# rank.abundances <- arrange(rank.abundances,SAD.model,year,plot,desc(rel.abund))
# rank.all <- bind_rows(rank.obs,rank.abundances)
# 
# write.csv2(rank.all,file = "./results/rank_abundances_complete.csv",row.names = FALSE)

