# obtain alpha and lambda values for all species, 
# considering all years together except one

# read and prepare data ---------------------------------------------------
library(tidyverse)
abund <- read.csv2("../Caracoles/data/abundances.csv",header = TRUE,stringsAsFactors = FALSE)
comp <- read.csv2("../Caracoles/data/competition_wide.csv",header = TRUE,stringsAsFactors = FALSE)
sp.rates <- read.csv2("../Caracoles/data/plant_species_traits.csv",header = TRUE,stringsAsFactors = FALSE)
sp.valid <- sp.rates$species.code[which(!is.na(sp.rates$germination.rate))]
precip <- read.csv2("../Caracoles/data/precipitation_per_year.csv",header = TRUE,stringsAsFactors = FALSE)

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
precip.list <- list()
for(i.sp in 1:length(focal.sp)){
  my.data <- subset(comp,focal == focal.sp[i.sp])
  precip.list[[i.sp]] <- data.frame(precipitaion = precip$prec[match(my.data$year,precip$hidrologic.year)])
  names(precip.list)[i.sp] <- focal.sp[i.sp]
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
# We constrain intraspecific alphas to be positive,
# but interspecific ones can take any value

source("./R/AP_pm_alpha_pairwise_lambdacov_global_alphacov_global.R")
source("./R/AP_project_alpha_pairwise_lambdacov_global_alphacov_global.R")

model_family <- "AP"

optimization_method <- "bobyqa"
alpha_form <- "pairwise"
lambda_cov_form <- "global"
alpha_cov_form <- "global"

fixed_terms <- NULL
bootstrap_samples <- 0

# to remove niche diff, we need positive coefficients

initial_values <- list(lambda = 10, 
                       alpha_intra = 0.01, 
                       alpha_inter = 0.01,
                       lambda_cov = 0, 
                       alpha_cov = 0)

lower_bounds <- list(lambda = 0, 
                     alpha_intra = 0,
                     alpha_inter = 0,
                     lambda_cov = 0,
                     alpha_cov = 0)

upper_bounds <- list(lambda = 1e3, 
                     alpha_intra = 1,
                     alpha_inter = 1,
                     lambda_cov = 1,
                     alpha_cov = 1)

# project only one timestep
timesteps <- 2

# fit model ---------------------------------------------------------------

# include environmental covariates? not for now
caracoles.fit <- cxr_pm_multifit(data = neigh.list,
                                 model_family = model_family,
                                 focal_column = all.sp,
                                 covariates = precip.list,
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

years <- sort(unique(base.abund$year))
plots <- sort(unique(base.abund$plot))
subplots <- sort(unique(base.abund$subplot))

nf.abund <- expand.grid(year = years,
                        plot = plots,
                        subplot = subplots, 
                        sp = all.sp)
nf.abund$individuals <- NA_real_
nf.abund$SAD.model <- "niche_fitness"

# niche + fitness diff
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
        sp.abund <- base.abund$individuals[base.abund$year == years[i.year] &
                                             base.abund$plot == plots[i.plot] &
                                             base.abund$subplot == subplots[i.sub] &
                                             base.abund$species %in% present.sp]
        names(sp.abund) <- present.sp
        
        # lambda_cov (named matrix with covariates in columns and taxa in rows)
        lambda_cov <- as.matrix(caracoles.fit$lambda_cov)
        lambda_cov <- matrix(lambda_cov[which(rownames(lambda_cov) %in% present.sp)],
                             nrow = length(present.sp),
                             dimnames = list(present.sp,"precipitation"))
        
        # alpha_cov (list of one element with a single alpha_cov value -- per species --)
        alpha_cov_data <- caracoles.fit$alpha_cov[[1]]
        alpha_cov_data <- alpha_cov_data[present.sp,present.sp]
        # alpha_cov_data <- matrix(alpha_cov_data,
        #                      nrow = length(present.sp),
        #                      dimnames = list(present.sp,"precipitation"))
        alpha_cov <- list(precipitation = alpha_cov_data)
        
        # covariates (matrix: covariates in columns, timesteps in rows)
        covariates <- matrix(precip$prec[precip$hidrologic.year %in% c(as.numeric(years[i.year]),as.numeric(years[i.year])+1)],nrow = 2)
        colnames(covariates) <- "precipitation"
        
        # project abundances
        sub.abund <- abundance_projection(lambda = sp.lambda,
                                          alpha_matrix = sp.alpha,
                                          model_family = model_family,
                                          alpha_form = alpha_form,
                                          lambda_cov_form = lambda_cov_form,
                                          alpha_cov_form = alpha_cov_form,
                                          lambda_cov = lambda_cov,
                                          alpha_cov = alpha_cov,
                                          covariates = covariates,
                                          timesteps = timesteps,
                                          initial_abundances = sp.abund)
        
        # match obtained abundances to the results dataframe
        pos <- which(nf.abund$year == years[i.year] &
                       nf.abund$plot == plots[i.plot] & 
                       nf.abund$subplot == subplots[i.sub])
        nf.abund$individuals[pos[which(nf.abund$sp[pos] %in% present.sp)]] <- sub.abund[nrow(sub.abund),]
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
        
        # lambda_cov (named matrix with covariates in columns and taxa in rows)
        lambda_cov <- as.matrix(caracoles.fit$lambda_cov)
        lambda_cov <- matrix(lambda_cov[which(rownames(lambda_cov) %in% present.sp)],
                             nrow = length(present.sp),
                             dimnames = list(present.sp,"precipitation"))
        
        # alpha_cov (list of one element with a single alpha_cov value -- per species --)
        alpha_cov_data <- caracoles.fit$alpha_cov[[1]]
        alpha_cov_data <- alpha_cov_data[present.sp,present.sp]
        # alpha_cov_data <- matrix(alpha_cov_data,
        #                      nrow = length(present.sp),
        #                      dimnames = list(present.sp,"precipitation"))
        alpha_cov <- list(precipitation = alpha_cov_data)
        
        # covariates (matrix: covariates in columns, timesteps in rows)
        covariates <- matrix(precip$prec[precip$hidrologic.year %in% c(as.numeric(years[i.year]),as.numeric(years[i.year])+1)],nrow = 2)
        colnames(covariates) <- "precipitation"
        
        # project abundances
        sub.abund <- abundance_projection(lambda = sp.lambda,
                                          alpha_matrix = sp.alpha.f,
                                          model_family = model_family,
                                          alpha_form = alpha_form,
                                          lambda_cov_form = lambda_cov_form,
                                          alpha_cov_form = alpha_cov_form,
                                          lambda_cov = lambda_cov,
                                          alpha_cov = alpha_cov,
                                          covariates = covariates,
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
        
        # lambda_cov (named matrix with covariates in columns and taxa in rows)
        lambda_cov <- as.matrix(caracoles.fit$lambda_cov)
        lambda_cov <- matrix(lambda_cov[which(rownames(lambda_cov) %in% present.sp)],
                             nrow = length(present.sp),
                             dimnames = list(present.sp,"precipitation"))
        
        # alpha_cov (list of one element with a single alpha_cov value -- per species --)
        alpha_cov_data <- caracoles.fit$alpha_cov[[1]]
        alpha_cov_data <- alpha_cov_data[present.sp,present.sp]
        # alpha_cov_data <- matrix(alpha_cov_data,
        #                      nrow = length(present.sp),
        #                      dimnames = list(present.sp,"precipitation"))
        alpha_cov <- list(precipitation = alpha_cov_data)
        
        # covariates (matrix: covariates in columns, timesteps in rows)
        covariates <- matrix(precip$prec[precip$hidrologic.year %in% c(as.numeric(years[i.year]),as.numeric(years[i.year])+1)],nrow = 2)
        colnames(covariates) <- "precipitation"
        
        # project abundances
        sub.abund <- abundance_projection(lambda = sp.lambda,
                                          alpha_matrix = sp.alpha,
                                          model_family = model_family,
                                          alpha_form = alpha_form,
                                          lambda_cov_form = lambda_cov_form,
                                          alpha_cov_form = alpha_cov_form,
                                          lambda_cov = lambda_cov,
                                          alpha_cov = alpha_cov,
                                          covariates = covariates,
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

# save results ------------------------------------------------------------

projected.abund <- dplyr::bind_rows(nf.abund,of.abund,on.abund)

write.csv2(projected.abund,file = "./results/projected_abundances_subplot_precip.csv",row.names = FALSE)




