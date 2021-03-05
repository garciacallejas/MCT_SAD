
library(tidyverse)
source("R/hill_diversity.R")

# read data ---------------------------------------------------------------

abund.obs <- read.csv2("results/abund_filtered.csv",header = TRUE,
                       stringsAsFactors = FALSE)

# this just to constrain abundances to the same set of sp as the other files
alpha.df <- read.csv2("results/alpha_RK.csv",stringsAsFactors = FALSE,
                      row.names = 1)
alpha.matrix <- as.matrix(alpha.df)

# sp.rates <- read.csv2("../Caracoles/data/plant_species_traits.csv",
#                       header = TRUE,stringsAsFactors = FALSE)
# sp.valid <- sp.rates$species.code[which(!is.na(sp.rates$germination.rate))]

# extract components ------------------------------------------------------

plot.obs <- abund.obs %>%
  filter(species %in% rownames(alpha.matrix)) %>%
  filter(!is.na(individuals) &
           individuals > 0) %>%
  group_by(year,plot, subplot,species) %>%
  summarise(mean.abund = mean(individuals),sum.abund = sum(individuals))

plot.obs <- plot.obs[,c("year","plot","subplot","species",
                        "mean.abund","sum.abund")]

years <- sort(unique(plot.obs$year))
plots <- sort(unique(plot.obs$plot))
subplots <- sort(unique(plot.obs$subplot))
components <- c("abundance","richness","evenness")

rad.components <- expand.grid(year = years,
                              plot = plots,
                              subplot = subplots,
                              richness = NA_real_,
                              abundance = NA_real_,
                              evenness = NA_real_)

for(i.year in 1:length(years)){
  for(i.plot in 1:length(plots)){
    for(i.subplot in 1:length(subplots)){
      my.rad <- subset(plot.obs,year == years[i.year] & 
                         plot == plots[i.plot] &
                         subplot == subplot[i.subplot])
      
      pos <- which(rad.components$year == years[i.year] & 
                     rad.components$plot == plots[i.plot] &
                     rad.components$subplot == subplots[i.subplot])
      
      rad.components$abundance[pos] <- sum(my.rad$sum.abund)
      rad.components$richness[pos] <- nrow(my.rad)
      rad.components$evenness[pos] <- hill.diversity(my.rad$sum.abund)
    }
  }
}

# store results -----------------------------------------------------------

rad.components.long <- gather(rad.components, 
                              key = "metric", 
                              value = "value", 
                              abundance,richness,evenness)

rad.components.long$type <- "observed"

write.csv2(rad.components.long,file = "results/observed_SAD_components_subplot.csv",
           row.names = FALSE)
