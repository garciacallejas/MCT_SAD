# obtain alpha and lambda values for all species, 
# considering all years together except one

# read and prepare data ---------------------------------------------------
library(tidyverse)
abund <- read.csv2("../Caracoles/data/abundances.csv",header = TRUE,stringsAsFactors = FALSE)
comp <- read.csv2("../Caracoles/data/competition_wide.csv",header = TRUE,stringsAsFactors = FALSE)
sp.rates <- read.csv2("../Caracoles/data/plant_species_traits.csv",header = TRUE,stringsAsFactors = FALSE)
sp.valid <- sp.rates$species.code[which(!is.na(sp.rates$germination.rate))]

base.abund <- abund %>% 
  filter(species %in% sp.valid) 
# group_by(species) %>% 
# summarise(num = sum(individuals))
# year.off <- 2019

comp <- comp[,which(names(comp) %in% c("year","plot","subplot","focal","fruit","seed",sp.valid))]
comp <- subset(comp, seed > 0)
# comp <- subset(comp, year != year.off)

all.sp <- names(comp)[7:length(names(comp))]
focal.sp <- unique(comp$focal[which(comp$focal %in% all.sp)])

# create multifit list
neigh.list <- list()
for(i.sp in 1:length(focal.sp)){
  my.data <- subset(comp,focal == focal.sp[i.sp])
  # if(nrow(my.data)>0){
  my.data <- my.data[,6:length(my.data)]
  names(my.data)[1] <- "fitness"
  neigh.list[[i.sp]] <- my.data
  names(neigh.list)[i.sp] <- focal.sp[i.sp]
  # }# if observations
}

# names(neigh.list) <- all.sp

# initial values ----------------------------------------------------------

library(cxr)
source("R/removeNiches.R")

# load model
# this is an Annual Plant model, 
# using the Ricker model from Mayfield and Stouffer 2017
# We constrain alphas to be positive,

source("./R/AP_pm_alpha_pairwise_lambdacov_none_alphacov_none.R")
source("./R/AP_project_alpha_pairwise_lambdacov_none_alphacov_none.R")

model_family <- "AP"

optimization_method <- "bobyqa"
alpha_form <- "pairwise"
lambda_cov_form <- "none"
alpha_cov_form <- "none"

fixed_terms <- NULL
bootstrap_samples <- 0

initial_values <- list(lambda = 10, alpha_intra = 0.01, alpha_inter = 0.01)
# to remove niche diff, we need positive coefficients
lower_bounds <- list(lambda = 0, alpha_intra = 0, alpha_inter = 0)
upper_bounds <- list(lambda = 1e3, alpha_intra = 1, alpha_inter = 1)

timesteps <- 4

# fit model ---------------------------------------------------------------

# include environmental covariates? not for now
caracoles.fit <- cxr_pm_multifit(data = neigh.list,
                                 model_family = model_family,
                                 focal_column = all.sp,
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

# TODO: bind together all models in a single loop

years <- sort(unique(base.abund$year))
plots <- sort(unique(base.abund$plot))
subplots <- sort(unique(base.abund$subplot))

nf.abund <- expand.grid(base.year = years,
                        predicted.year = years+1,
                        plot = plots,
                        subplot = subplots, 
                        sp = all.sp)
nf.abund$individuals <- NA_real_
nf.abund$SAD.model <- "niche_fitness"

nf.abund <- subset(nf.abund,predicted.year > base.year & predicted.year < 2020)

# niche + fitness diff
for(i.year in 1:(length(years)-1)){
  for(i.plot in 1:length(plots)){
    for(i.sub in 1:length(subplots)){
      
      # species and parameters in this subplot
      # 1 - abundance
      present.abund <- base.abund$species[base.abund$year == years[i.year] &
                                            base.abund$plot == plots[i.plot] &
                                            base.abund$subplot == subplots[i.sub] &
                                            base.abund$individuals > 0]
      # 2 - focal
      present.focal <- unique(comp$focal[comp$year == years[i.year] &
                                           comp$plot == plots[i.plot] &
                                           comp$subplot == subplots[i.sub]])
      present.sp <- sort(unique(intersect(present.abund,present.focal)))
      
      if(length(present.sp)>1){
        
        sp.lambda <- caracoles.fit$lambda[present.sp]
        sp.alpha <- caracoles.fit$alpha_matrix[present.sp,present.sp]
        sp.abund <- base.abund$individuals[base.abund$year == years[i.year] &
                                             base.abund$plot == plots[i.plot] &
                                             base.abund$subplot == subplots[i.sub] &
                                             base.abund$species %in% present.sp]
        names(sp.abund) <- present.sp
        
        # project abundances
        sub.abund <- abundance_projection(lambda = sp.lambda,
                                          alpha_matrix = sp.alpha,
                                          model_family = model_family,
                                          alpha_form = alpha_form,
                                          lambda_cov_form = lambda_cov_form,
                                          alpha_cov_form = alpha_cov_form,
                                          timesteps = timesteps+1, # because t1 is the original
                                          initial_abundances = sp.abund)
        predicted.years <- years[i.year]:(years[i.year]+timesteps)
        row.names(sub.abund) <- predicted.years
        
        # since the first year on sub.abund is not actually a prediction,
        # remove it from the predicted.years vector. 
        # Remove years > 2020 as well
        predicted.years <- predicted.years[predicted.years < 2020]
        predicted.years <- predicted.years[predicted.years > years[i.year]]
        
        # match obtained abundances to the results dataframe
        for(i.pred in 1:length(predicted.years)){
          
          pos <- which(nf.abund$base.year == years[i.year] &
                         nf.abund$plot == plots[i.plot] &
                         nf.abund$subplot == subplots[i.sub] &
                         nf.abund$predicted.year == predicted.years[i.pred])
          nf.abund$individuals[pos[which(nf.abund$sp[pos] %in% present.sp)]] <- 
          sub.abund[which(row.names(sub.abund) == predicted.years[i.pred]),]
        }# for i.timestep
        
      }# if >1 sp
    }# i.subplot
  }# i.plot
}# i.year

# only fitness diff

of.abund <- expand.grid(year = years,
                        plot = plots,
                        subplot = subplots, 
                        sp = all.sp)
of.abund$individuals <- NA_real_
of.abund$SAD.model <- "only.fitness"

for(i.year in 1:length(years)){
  for(i.plot in 1:length(plots)){
    for(i.sub in 1:length(subplots)){
      
      # species and parameters in this subplot
      # 1 - abundance
      present.abund <- base.abund$species[base.abund$year == years[i.year] &
                                            base.abund$plot == plots[i.plot] &
                                            base.abund$subplot == subplots[i.sub] &
                                            base.abund$individuals > 0]
      # 2 - focal
      present.focal <- unique(comp$focal[comp$year == years[i.year] &
                                           comp$plot == plots[i.plot] &
                                           comp$subplot == subplots[i.sub]])
      present.sp <- sort(unique(intersect(present.abund,present.focal)))
      
      if(length(present.sp)>1){
        sp.lambda <- caracoles.fit$lambda[present.sp]
        sp.alpha <- caracoles.fit$alpha_matrix[present.sp,present.sp]
        
        # modify to nxn sp
        sp.alpha.f <- removeNiches_nsp(sp.alpha)
        
        sp.abund <- base.abund$individuals[base.abund$year == years[i.year] &
                                             base.abund$plot == plots[i.plot] &
                                             base.abund$subplot == subplots[i.sub] &
                                             base.abund$species %in% present.sp]
        names(sp.abund) <- present.sp
        
        # project abundances
        sub.abund <- abundance_projection(lambda = sp.lambda,
                                          alpha_matrix = sp.alpha.f,
                                          model_family = model_family,
                                          alpha_form = alpha_form,
                                          lambda_cov_form = lambda_cov_form,
                                          alpha_cov_form = alpha_cov_form,
                                          timesteps = timesteps,
                                          initial_abundances = sp.abund)
        
        # match obtained abundances to the results dataframe
        pos <- which(of.abund$year == years[i.year] &
                       of.abund$plot == plots[i.plot] & 
                       of.abund$subplot == subplots[i.sub])
        of.abund$individuals[pos[which(of.abund$sp[pos] %in% present.sp)]] <- sub.abund[nrow(sub.abund),]
      }# if >1 sp
    }# i.subplot
  }# i.plot
}# i.year

# only niche diff

on.abund <- expand.grid(year = years,
                        plot = plots,
                        subplot = subplots, 
                        sp = all.sp)
on.abund$individuals <- NA_real_
on.abund$SAD.model <- "only.niche"

for(i.year in 1:length(years)){
  for(i.plot in 1:length(plots)){
    for(i.sub in 1:length(subplots)){
      
      # species and parameters in this subplot
      # 1 - abundance
      present.abund <- base.abund$species[base.abund$year == years[i.year] &
                                            base.abund$plot == plots[i.plot] &
                                            base.abund$subplot == subplots[i.sub] &
                                            base.abund$individuals > 0]
      # 2 - focal
      present.focal <- unique(comp$focal[comp$year == years[i.year] &
                                           comp$plot == plots[i.plot] &
                                           comp$subplot == subplots[i.sub]])
      present.sp <- sort(unique(intersect(present.abund,present.focal)))
      
      if(length(present.sp)>1){
        
        sp.lambda <- rep(mean(caracoles.fit$lambda[present.sp]),length(caracoles.fit$lambda[present.sp]))
        names(sp.lambda) <- present.sp
        # caracoles.fit$lambda[present.sp]
        
        sp.alpha <- caracoles.fit$alpha_matrix[present.sp,present.sp]
        
        # modify to nxn sp
        # sp.alpha.f <- removeNiches_nsp(sp.alpha)
        
        sp.abund <- base.abund$individuals[base.abund$year == years[i.year] &
                                             base.abund$plot == plots[i.plot] &
                                             base.abund$subplot == subplots[i.sub] &
                                             base.abund$species %in% present.sp]
        names(sp.abund) <- present.sp
        
        # project abundances
        sub.abund <- abundance_projection(lambda = sp.lambda,
                                          alpha_matrix = sp.alpha,
                                          model_family = model_family,
                                          alpha_form = alpha_form,
                                          lambda_cov_form = lambda_cov_form,
                                          alpha_cov_form = alpha_cov_form,
                                          timesteps = timesteps,
                                          initial_abundances = sp.abund)
        
        # match obtained abundances to the results dataframe
        pos <- which(on.abund$year == years[i.year] &
                       on.abund$plot == plots[i.plot] & 
                       on.abund$subplot == subplots[i.sub])
        on.abund$individuals[pos[which(on.abund$sp[pos] %in% present.sp)]] <- sub.abund[nrow(sub.abund),]
      }# if >1 sp
    }# i.subplot
  }# i.plot
}# i.year

# increased dispersal -----------------------------------------------------
# relative to abundance in previous year

# maximum relative increase in fecundity due to dispersal
increase.rate <- 2

disp.abund <- expand.grid(year = years,
                          plot = plots,
                          subplot = subplots, 
                          sp = all.sp)
disp.abund$individuals <- NA_real_
disp.abund$SAD.model <- "increased.dispersal"

for(i.year in 1:length(years)){
  for(i.plot in 1:length(plots)){
    for(i.sub in 1:length(subplots)){
      
      # species and parameters in this subplot
      # 1 - abundance
      present.abund <- base.abund$species[base.abund$year == years[i.year] &
                                            base.abund$plot == plots[i.plot] &
                                            base.abund$subplot == subplots[i.sub] &
                                            base.abund$individuals > 0]
      # 2 - focal
      present.focal <- unique(comp$focal[comp$year == years[i.year] &
                                           comp$plot == plots[i.plot] &
                                           comp$subplot == subplots[i.sub]])
      present.sp <- sort(unique(intersect(present.abund,present.focal)))
      
      if(length(present.sp)>1){
        
        # base fecundity
        base.lambda <- caracoles.fit$lambda[present.sp]
        
        # previous abundance across all caracoles
        prev.full.abund <- base.abund %>%
          filter(year == years[i.year] & species %in% present.sp) %>%
          group_by(species) %>%
          summarise(sum.ind = sum(individuals))
        
        # the species with the highest abundance gets the maximum increase rate
        max.sp <- prev.full.abund$species[which(prev.full.abund$sum.ind == max(prev.full.abund$sum.ind))] 
        # from the maximum increase and previous abundance
        # I can extract the rate of increase in lambda per unit of abundance
        unit.increase <- (base.lambda[which(names(base.lambda) == max.sp)] * increase.rate)/prev.full.abund$sum.ind[prev.full.abund$species == max.sp]
        
        sp.lambda <- base.lambda + unit.increase * prev.full.abund$sum.ind
        
        sp.alpha <- caracoles.fit$alpha_matrix[present.sp,present.sp]
        
        # modify to nxn sp
        # sp.alpha.f <- removeNiches_nsp(sp.alpha)
        
        sp.abund <- base.abund$individuals[base.abund$year == years[i.year] &
                                             base.abund$plot == plots[i.plot] &
                                             base.abund$subplot == subplots[i.sub] &
                                             base.abund$species %in% present.sp]
        names(sp.abund) <- present.sp
        
        # project abundances
        sub.abund <- abundance_projection(lambda = sp.lambda2,
                                          alpha_matrix = sp.alpha,
                                          model_family = model_family,
                                          alpha_form = alpha_form,
                                          lambda_cov_form = lambda_cov_form,
                                          alpha_cov_form = alpha_cov_form,
                                          timesteps = timesteps,
                                          initial_abundances = sp.abund)
        
        # match obtained abundances to the results dataframe
        pos <- which(disp.abund$year == years[i.year] &
                       disp.abund$plot == plots[i.plot] & 
                       disp.abund$subplot == subplots[i.sub])
        disp.abund$individuals[pos[which(disp.abund$sp[pos] %in% present.sp)]] <- sub.abund[nrow(sub.abund),]
      }# if >1 sp
    }# i.subplot
  }# i.plot
}# i.year

# save results ------------------------------------------------------------

projected.abund <- dplyr::bind_rows(nf.abund,of.abund,on.abund,disp.abund)

write.csv2(projected.abund,file = "./results/projected_abundances_subplot.csv",row.names = FALSE)




