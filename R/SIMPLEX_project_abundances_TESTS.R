# obtain alpha and lambda values for all species, 
# considering all years together except one

# read and prepare data ---------------------------------------------------
library(tidyverse)
library(cxr)

abund <- read.csv2("../Caracoles/data/abundances.csv",header = TRUE,stringsAsFactors = FALSE)
comp <- read.csv2("../Caracoles/data/competition_wide.csv",header = TRUE,stringsAsFactors = FALSE)
sp.rates <- read.csv2("../Caracoles/data/plant_species_traits.csv",header = TRUE,stringsAsFactors = FALSE)
sp.valid <- sp.rates$species.code[which(!is.na(sp.rates$germination.rate))]

# remove LYTR, which has too few observations/cannot be fitted properly
abund <- subset(abund,species != "LYTR")
comp <- subset(comp, focal != "LYTR")
comp$LYTR <- NULL
sp.valid <- sp.valid[which(sp.valid != "LYTR")]

base.abund <- abund %>% 
  filter(species %in% sp.valid) 
# group_by(species) %>% 
# summarise(num = sum(individuals))
# year.off <- 2019

comp <- comp[,which(names(comp) %in% c("year","plot","subplot","focal","fruit","seed",sp.valid))]
comp <- subset(comp, seed > 0)

all.sp <- names(comp)[7:length(names(comp))]
focal.sp <- sort(unique(comp$focal[which(comp$focal %in% all.sp)]))

plots <- sort(unique(base.abund$plot))
subplots <- sort(unique(base.abund$subplot))

# initial values ----------------------------------------------------------

# load model
# this is an Annual Plant model, 
# using the Ricker model from Mayfield and Stouffer 2017
# We constrain intraspecific alphas to be positive,
# but interspecific ones can take any value

source("./R/AP_pm_alpha_pairwise_lambdacov_none_alphacov_none.R")
source("./R/AP_project_alpha_pairwise_lambdacov_none_alphacov_none.R")

model_family <- "AP"

optimization_method <- "bobyqa"
alpha_form <- "pairwise"
lambda_cov_form <- "none"
alpha_cov_form <- "none"

fixed_terms <- NULL
bootstrap_samples <- 0 # WHEN OK, SET TO 100

initial_values <- list(lambda = 10, alpha_intra = 0.01, alpha_inter = 0.01)
# to remove niche diff, we need positive coefficients
lower_bounds <- list(lambda = 0, alpha_intra = 0, alpha_inter = 0)
upper_bounds <- list(lambda = 1e3, alpha_intra = 1, alpha_inter = 1)

timesteps <- 2

# project one year based on all others ------------------------------------
# except for 2015, for which we don't have initial abunances

pred.years <- 2016:2019

projected.abundances <- NULL
for(i.year in 1:length(pred.years)){
  
  comp.training <- subset(comp, year == pred.years[i.year])
  all.sp.training <- names(comp.training)[7:length(names(comp.training))]
  
  # create multifit list
  neigh.list <- list()
  for(i.sp in 1:length(focal.sp)){
    my.data <- subset(comp.training,focal == focal.sp[i.sp])
    if(nrow(my.data)>10){
    my.data <- my.data[,6:length(my.data)]
    names(my.data)[1] <- "fitness"
    neigh.list[[i.sp]] <- my.data
    names(neigh.list)[i.sp] <- focal.sp[i.sp]
    }# if observations
  }
  
  neigh.list = neigh.list[-which(sapply(neigh.list, is.null))]
  
  # fit model ---------------------------------------------------------------
  # without environmental covariates
  
  caracoles.fit <- cxr_pm_multifit(data = neigh.list,
                                   model_family = model_family,
                                   focal_column = names(neigh.list),
                                   covariates = NULL,
                                   optimization_method = optimization_method,
                                   alpha_form = alpha_form,
                                   lambda_cov_form = lambda_cov_form,
                                   alpha_cov_form = alpha_cov_form,
                                   initial_values = initial_values,
                                   lower_bounds = lower_bounds,
                                   upper_bounds = upper_bounds,
                                   fixed_terms = fixed_terms,
                                   bootstrap_samples = bootstrap_samples)
  
  
  # project abundances ------------------------------------------------------
  # at the subplot level
  
  proj.abund <- subset(base.abund,year == pred.years[i.year] & 
                 species %in% all.sp.training)
  proj.abund <- proj.abund[,c("year","plot","subplot","species","individuals")]
  names(proj.abund)[5] <- "observed"
  
  # proj.abund <- expand.grid(year = pred.years[i.year],
  #                           plot = plots,
  #                           subplot = subplots, 
  #                           species = all.sp.training)

  proj.abund$predicted <- NA_real_

    for(i.plot in 1:length(plots)){
      for(i.sub in 1:length(subplots)){
        
        # species and parameters in this subplot
        # 1 - abundance
        present.abund <- base.abund$species[base.abund$year == pred.years[i.year] &
                                              base.abund$plot == plots[i.plot] &
                                              base.abund$subplot == subplots[i.sub] &
                                              base.abund$individuals > 0]
        # 2 - focal
        present.focal <- unique(comp$focal[comp$year == pred.years[i.year] &
                                             comp$plot == plots[i.plot] &
                                             comp$subplot == subplots[i.sub]])
        present.sp <- sort(unique(intersect(present.abund,present.focal)))
        present.sp <- intersect(present.sp,names(neigh.list))
        
        if(length(present.sp)>1){
          
          sp.lambda <- caracoles.fit$lambda[present.sp]
          sp.alpha <- caracoles.fit$alpha_matrix[present.sp,present.sp]
          
          # CAREFUL!
          # NA to 0s, in case an interaction is not observed
          # so, assume that unobserved interactions
          # do not influence abundances
          sp.alpha[is.na(sp.alpha)] <- 0
          
          # abundance in the previous year as baseline
          prev.abund <- base.abund$individuals[base.abund$year == (pred.years[i.year] - 1) &
                                               base.abund$plot == plots[i.plot] &
                                               base.abund$subplot == subplots[i.sub] &
                                               base.abund$species %in% present.sp]
          names(prev.abund) <- present.sp
          
          # project abundances
          sub.abund <- abundance_projection(lambda = sp.lambda,
                                            alpha_matrix = sp.alpha,
                                            model_family = model_family,
                                            alpha_form = alpha_form,
                                            lambda_cov_form = lambda_cov_form,
                                            alpha_cov_form = alpha_cov_form,
                                            timesteps = timesteps,
                                            initial_abundances = prev.abund)
          
          # match obtained abundances to the results dataframe
          pos <- which(proj.abund$year == pred.years[i.year] &
                         proj.abund$plot == plots[i.plot] & 
                         proj.abund$subplot == subplots[i.sub])
          
          proj.abund$predicted[pos[which(proj.abund$species[pos] %in% present.sp)]] <- sub.abund[nrow(sub.abund),]
          
        }# if >1 sp
      }# i.subplot
    }# i.plot
  
  projected.abundances <- rbind(projected.abundances,proj.abund)
  
}# for pred.years

# save results ------------------------------------------------------------

write.csv2(projected.abundances,file = "./results/projected_abundances_TEST.csv",row.names = FALSE)




