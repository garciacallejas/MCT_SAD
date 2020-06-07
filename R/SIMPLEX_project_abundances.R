# obtain alpha and lambda values for all species, 
# with two possible training/test datasets

# 1 - all years but one as training, single year as test
# 2 - 80% spatial observations training, 20% test

# training.data.type <- "yearly"
training.data.type <- "spatial"

# also, using either no covariates or precipitation as covariate (model 4)

include.precipitation <- TRUE

# read and prepare data ---------------------------------------------------
library(tidyverse)
library(cxr)

abund <- read.csv2("../Caracoles/data/abundances.csv",header = TRUE,stringsAsFactors = FALSE)
comp <- read.csv2("../Caracoles/data/competition_wide.csv",header = TRUE,stringsAsFactors = FALSE)
sp.rates <- read.csv2("../Caracoles/data/plant_species_traits.csv",header = TRUE,stringsAsFactors = FALSE)
sp.valid <- sp.rates$species.code[which(!is.na(sp.rates$germination.rate))]
precip <- read.csv2("../Caracoles/data/precipitation_per_year.csv",header = TRUE,stringsAsFactors = FALSE)

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
focal.sp <- unique(comp$focal[which(comp$focal %in% all.sp)])

plots <- sort(unique(base.abund$plot))
subplots <- sort(unique(base.abund$subplot))

# initial values ----------------------------------------------------------

# load model
# this is an Annual Plant model, 
# using the Ricker model from Mayfield and Stouffer 2017

source("./R/AP_pm_alpha_pairwise_lambdacov_none_alphacov_none.R")
source("./R/AP_project_alpha_pairwise_lambdacov_none_alphacov_none.R")
source("./R/AP_pm_alpha_pairwise_lambdacov_global_alphacov_pairwise.R")
source("./R/AP_project_alpha_pairwise_lambdacov_global_alphacov_pairwise.R")
source("./R/AP_pm_alpha_pairwise_lambdacov_global_alphacov_global.R")
source("./R/AP_project_alpha_pairwise_lambdacov_global_alphacov_global.R")

model_family <- "AP"

optimization_method <- "bobyqa"
alpha_form <- "pairwise"

if(include.precipitation){
  lambda_cov_form <- "global"
  alpha_cov_form <- "global"
  
  initial_values <- list(lambda = 10, 
                         alpha_intra = 0.1, 
                         alpha_inter = 0,
                         lambda_cov = 0, 
                         alpha_cov = 0)
  
  lower_bounds <- list(lambda = 0, 
                       alpha_intra = 0,
                       alpha_inter = -1,
                       lambda_cov = 0,
                       alpha_cov = 0)
  
  upper_bounds <- list(lambda = 1e3, 
                       alpha_intra = 1,
                       alpha_inter = 1,
                       lambda_cov = 1,
                       alpha_cov = 1)
  
}else{
  lambda_cov_form <- "none"
  alpha_cov_form <- "none"
  
  initial_values <- list(lambda = 10, alpha_intra = 0.01, alpha_inter = 0.01)
  
  lower_bounds <- list(lambda = 0, alpha_intra = 0, alpha_inter = -1)
  upper_bounds <- list(lambda = 1e3, alpha_intra = 1, alpha_inter = 1)
}

fixed_terms <- NULL
bootstrap_samples <- 0 # WHEN OK, SET TO 100

timesteps <- 2

# project one year based on all others ------------------------------------
# except for 2015, for which we don't have initial abunances

if(training.data.type == "yearly"){
  
  pred.years <- 2016:2019
  
  projected.abundances <- NULL
  for(i.year in 1:length(pred.years)){
    
    comp.training <- subset(comp, year != pred.years[i.year])
    all.sp.training <- names(comp.training)[7:length(names(comp.training))]
    
    # create multifit list
    neigh.list <- list()
    precip.list <- list()
    
    for(i.sp in 1:length(focal.sp)){
      my.data <- subset(comp.training,focal == focal.sp[i.sp])
      
      if(include.precipitation){
        precip.list[[i.sp]] <- data.frame(precipitaion = precip$prec[match(my.data$year,precip$hidrologic.year)])
        names(precip.list)[i.sp] <- focal.sp[i.sp]
      }
      
      my.data <- my.data[,6:length(my.data)]
      names(my.data)[1] <- "fitness"
      neigh.list[[i.sp]] <- my.data
      names(neigh.list)[i.sp] <- focal.sp[i.sp]
    }
    
    # fit model ---------------------------------------------------------------
    
    # data = neigh.list
    # focal_column = all.sp.training
    if(include.precipitation){
      covariates <- precip.list
    }else{
      covariates <- NULL
    }
    
    caracoles.fit <- cxr_pm_multifit(data = neigh.list,
                                     model_family = model_family,
                                     focal_column = all.sp.training,
                                     covariates = covariates,
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
        
        if(length(present.sp)>1){
          
          sp.lambda <- caracoles.fit$lambda[present.sp]
          sp.alpha <- caracoles.fit$alpha_matrix[present.sp,present.sp]
          
          # abundance in the previous year as baseline
          prev.abund <- base.abund$individuals[base.abund$year == (pred.years[i.year] - 1) &
                                                 base.abund$plot == plots[i.plot] &
                                                 base.abund$subplot == subplots[i.sub] &
                                                 base.abund$species %in% present.sp]
          names(prev.abund) <- present.sp
          
          # project abundances
          
          if(include.precipitation){
            
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
            covariates <- matrix(precip$prec[precip$hidrologic.year %in% c(as.numeric(pred.years[i.year])-1,as.numeric(pred.years[i.year]))],nrow = 2)
            colnames(covariates) <- "precipitation"
            
            # lambda = sp.lambda
            # alpha_matrix = sp.alpha
            # initial_abundances = prev.abund
            
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
                                              initial_abundances = prev.abund)
          }else{
            sub.abund <- abundance_projection(lambda = sp.lambda,
                                              alpha_matrix = sp.alpha,
                                              model_family = model_family,
                                              alpha_form = alpha_form,
                                              lambda_cov_form = lambda_cov_form,
                                              alpha_cov_form = alpha_cov_form,
                                              timesteps = timesteps,
                                              initial_abundances = prev.abund)
          }
          
          
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
  
}else if(training.data.type == "spatial"){
  
  # i need to divide training and test data
  # making sure that test data does not come
  # from the first year (2015)
  # because I need a t-1 abundance to run the projection
  
  # so, include the whole 2015 data in the training set
  comp.not.first <- subset(comp, year != "2015")
  
  num.test.data <- round(0.2*nrow(comp))
  freq.test <- num.test.data/nrow(comp.not.first)
  test.rows <- sample(nrow(comp.not.first),num.test.data)
  
  test.data <- comp.not.first[test.rows,]
  training.data <- comp.not.first[-test.rows,]
  training.data <- rbind(comp[comp$year == 2015,],training.data)
  
  projected.abundances <- NULL
  
  all.sp.training <- names(training.data)[7:length(names(training.data))]
  
  # create multifit list
  neigh.list <- list()
  precip.list <- list()
  
  for(i.sp in 1:length(focal.sp)){
    my.data <- subset(training.data,focal == focal.sp[i.sp])
    
    if(include.precipitation){
      precip.list[[i.sp]] <- data.frame(precipitaion = precip$prec[match(my.data$year,precip$hidrologic.year)])
      names(precip.list)[i.sp] <- focal.sp[i.sp]
    }
    
    my.data <- my.data[,6:length(my.data)]
    names(my.data)[1] <- "fitness"
    neigh.list[[i.sp]] <- my.data
    names(neigh.list)[i.sp] <- focal.sp[i.sp]
  }
  
  #  fit model ---------------------------------------------------------------
  
  if(include.precipitation){
    covariates <- precip.list
  }else{
    covariates <- NULL
  }
  
  # TEST
  # data = neigh.list
  # focal_column = all.sp.training
  # source("../cxr/R/cxr_check_initial_values.R")
  # source("../cxr/R/cxr_check_input_data.R")
  # source("../cxr/R/cxr_check_method_boundaries.R")
  # source("../cxr/R/cxr_check_pm_input.R")
  # source("../cxr/R/cxr_get_init_params.R")
  # source("../cxr/R/cxr_get_model_bounds.R")
  # source("../cxr/R/cxr_pm_bootstrap.R")
  # source("../cxr/R/cxr_retrieve_params.R")
  # source("../cxr/R/cxr_return_init_length.R")
  # source("../cxr/R/cxr_sort_params.R")
  # source("../cxr/R/BH_pm_alpha_pairwise_lambdacov_none_alphacov_none.R")
  # source("../cxr/R/BH_pm_alpha_pairwise_lambdacov_global_alphacov_global.R")
  # source("../cxr/R/cxr_pm_fit.R")
  # source("../cxr/R/cxr_pm_multifit.R")
  
  caracoles.fit <- cxr_pm_multifit(data = neigh.list,
                                   model_family = model_family,
                                   focal_column = all.sp.training,
                                   covariates = covariates,
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
  
  proj.abund <- left_join(test.data[,1:3],base.abund[,c(1,4:7)])
  # proj.abund <- proj.abund[,c("year","plot","subplot","species","individuals")]
  names(proj.abund)[5] <- "observed"
  
  # proj.abund <- expand.grid(year = pred.years[i.year],
  #                           plot = plots,
  #                           subplot = subplots, 
  #                           species = all.sp.training)
  
  proj.abund$predicted <- NA_real_
  
  id.char <- unique(paste(proj.abund$year,"_",proj.abund$plot,"_",proj.abund$subplot,sep=""))
  
  for(i.id in 1:length(id.char)){
    i.year <- substr(id.char[i.id],1,4)
    i.plot <- substr(id.char[i.id],6,6)
    i.sub <- substr(id.char[i.id],8,nchar(id.char[i.id]))    
    
    # species and parameters in this subplot
    # 1 - abundance
    present.abund <- base.abund$species[base.abund$year == i.year &
                                          base.abund$plot == i.plot &
                                          base.abund$subplot == i.sub &
                                          base.abund$individuals > 0]
    # 2 - focal
    present.focal <- unique(comp$focal[comp$year == i.year &
                                         comp$plot == i.plot &
                                         comp$subplot == i.sub])
    
    present.sp <- sort(unique(intersect(present.abund,present.focal)))
    
    if(length(present.sp)>1){
      
      sp.lambda <- caracoles.fit$lambda[present.sp]
      sp.alpha <- caracoles.fit$alpha_matrix[present.sp,present.sp]
      
      # abundance in the previous year as baseline
      prev.abund <- base.abund$individuals[base.abund$year == (as.integer(i.year) - 1) &
                                             base.abund$plot == i.plot &
                                             base.abund$subplot == i.sub &
                                             base.abund$species %in% present.sp]
      names(prev.abund) <- present.sp
      
      # project abundances
      
      if(include.precipitation){
        
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
        covariates <- matrix(precip$prec[precip$hidrologic.year %in% c(as.numeric(i.year)-1,as.numeric(i.year))],nrow = 2)
        colnames(covariates) <- "precipitation"
        
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
                                          initial_abundances = prev.abund)
      }else{
        sub.abund <- abundance_projection(lambda = sp.lambda,
                                          alpha_matrix = sp.alpha,
                                          model_family = model_family,
                                          alpha_form = alpha_form,
                                          lambda_cov_form = lambda_cov_form,
                                          alpha_cov_form = alpha_cov_form,
                                          timesteps = timesteps,
                                          initial_abundances = prev.abund)
      }
      

      
      # match obtained abundances to the results dataframe
      pos <- which(proj.abund$year == i.year &
                     proj.abund$plot == i.plot & 
                     proj.abund$subplot == i.sub)
      
      proj.abund$predicted[pos[which(proj.abund$species[pos] %in% present.sp)]] <- sub.abund[nrow(sub.abund),]
      
    }# if >1 sp
    
  }# for i.id
  
  projected.abundances <- rbind(projected.abundances,proj.abund)
  
}# if-else test data type

# save results ------------------------------------------------------------

precip.char <- ifelse(include.precipitation,"precip","no_precip")

write.csv2(projected.abundances,file = paste("./results/projected_abundances_",training.data.type,"_",precip.char,".csv",sep=""),row.names = FALSE)




