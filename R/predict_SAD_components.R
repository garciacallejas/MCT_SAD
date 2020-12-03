
library(tidyverse)
library(cxr)
source("R/hill_diversity.R")

# read data ---------------------------------------------------------------

load("results/communities.Rdata")

# some constants ----------------------------------------------------------

years <- names(communities)

# all years are predicted except the first
years.predicted <- years[2:length(years)]

# likewise, start by projecting abundances from all years except the most recent
# so that I can compare projected and observed
# in addition, note that all projections are from t to t+1, for now
initial.years <- years[1:(length(years)-1)]

plots <- 1:length(communities[[1]])
types <- names(communities[[1]][[1]])
steps <- length(communities[[1]][[1]][["nd"]][["alpha"]])

metrics <- c("richness","abundance","evenness")

# also, load the model
source("./R/AP_pm_alpha_pairwise_lambdacov_none_alphacov_none.R")
source("./R/AP_project_alpha_pairwise_lambdacov_none_alphacov_none.R")

model_family <- "AP"
optimization_method <- "bobyqa"
alpha_form <- "pairwise"
lambda_cov_form <- "none"
alpha_cov_form <- "none"

# results data structure --------------------------------------------------

obs.sad <- expand.grid(year.predicted = years.predicted,plot = plots,
                       metric = metrics,type = c("observed"),
                       value = NA)

pert.sad <- expand.grid(year.predicted = years.predicted,plot = plots,
                        metric = metrics,type = types, intensity = 1:steps,
                        value = NA)

# project abundances for every community ----------------------------------

for(i.year in 1:length(initial.years)){
  for(i.plot in 1:length(plots)){
    for(i.type in 1:length(types)){
      
      if(types[i.type] == "obs"){
        
        lambda.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[types[i.type]]][["lambda"]]
        sp.alpha <- communities[[initial.years[i.year]]][[plots[i.plot]]][[types[i.type]]][["alpha"]]
        abund.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[types[i.type]]][["abundances"]]
        
        if(sum(sum(is.na(lambda.df)),sum(is.na(sp.alpha)),sum(is.na(abund.df))) == 0){
          
          sp.lambda <- lambda.df$lambda
          names(sp.lambda) <- lambda.df$sp
          
          sp.abund <- abund.df$abundance
          names(sp.abund) <- abund.df$species
          
          sub.abund <- cxr::abundance_projection(lambda = sp.lambda,
                                                 alpha_matrix = sp.alpha,
                                                 model_family = model_family,
                                                 alpha_form = alpha_form,
                                                 lambda_cov_form = lambda_cov_form,
                                                 alpha_cov_form = alpha_cov_form,
                                                 timesteps = 2, # because t1 is the original
                                                 initial_abundances = sp.abund)
          projected.abund <- sub.abund[2,]
          projected.richness <- sum(projected.abund > 0)
          projected.abundance <- sum(projected.abund)
          projected.evenness <- hill.diversity(projected.abund)
          
          pos <- which(pert.sad$year.predicted == (as.numeric(initial.years[i.year]) + 1) &
                         pert.sad$plot == plots[i.plot] & 
                         pert.sad$type == "obs" )
          
          pert.sad$value[pos[which(pert.sad$metric[pos] == "richness")]] <- 
            projected.richness
          pert.sad$value[pos[which(pert.sad$metric[pos] == "abundance")]] <- 
            projected.abundance
          pert.sad$value[pos[which(pert.sad$metric[pos] == "evenness")]] <- 
            projected.evenness
        }# if !na
      }else if(types[i.type] == "nd"){
        
        for(i.step in 1:steps){
          
          lambda.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[types[i.type]]][["lambda"]]
          if(sum(is.na(lambda.df)) == 0){
            sp.alpha <- communities[[initial.years[i.year]]][[plots[i.plot]]][[types[i.type]]][["alpha"]][[i.step]]
          }else{
            sp.alpha <- NA
          }
          abund.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[types[i.type]]][["abundances"]]
          
          if(sum(sum(is.na(lambda.df)),sum(is.na(sp.alpha)),sum(is.na(abund.df))) == 0){
            
            sp.lambda <- lambda.df$lambda
            names(sp.lambda) <- lambda.df$sp
            
            sp.abund <- abund.df$abundance
            names(sp.abund) <- abund.df$species
            
            sub.abund <- cxr::abundance_projection(lambda = sp.lambda,
                                                   alpha_matrix = sp.alpha,
                                                   model_family = model_family,
                                                   alpha_form = alpha_form,
                                                   lambda_cov_form = lambda_cov_form,
                                                   alpha_cov_form = alpha_cov_form,
                                                   timesteps = 2, # because t1 is the original
                                                   initial_abundances = sp.abund)
            projected.abund <- sub.abund[2,]
            projected.richness <- sum(projected.abund > 0)
            projected.abundance <- sum(projected.abund)
            projected.evenness <- hill.diversity(projected.abund)
            
            pos <- which(pert.sad$year.predicted == (as.numeric(initial.years[i.year]) + 1) &
                           pert.sad$plot == plots[i.plot] & 
                           pert.sad$type == types[i.type] &
                           pert.sad$intensity == i.step)
            
            pert.sad$value[pos[which(pert.sad$metric[pos] == "richness")]] <- 
              projected.richness
            pert.sad$value[pos[which(pert.sad$metric[pos] == "abundance")]] <- 
              projected.abundance
            pert.sad$value[pos[which(pert.sad$metric[pos] == "evenness")]] <- 
              projected.evenness
          }# if !na
        }# for i.step
      }else if(types[i.type] == "fd"){
        
        for(i.step in 1:steps){
          
          sp.alpha <- communities[[initial.years[i.year]]][[plots[i.plot]]][[types[i.type]]][["alpha"]]
          abund.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[types[i.type]]][["abundances"]]
          if(sum(is.na(sp.alpha)) == 0){
            lambda.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[types[i.type]]][["lambda"]][[i.step]]
          }else{
            sp.alpha <- NA
          }
          
          if(sum(sum(is.na(lambda.df)),sum(is.na(sp.alpha)),sum(is.na(abund.df))) == 0){
            
            sp.lambda <- lambda.df$lambda
            names(sp.lambda) <- lambda.df$sp
            
            sp.abund <- abund.df$abundance
            names(sp.abund) <- abund.df$species
            
            sub.abund <- cxr::abundance_projection(lambda = sp.lambda,
                                                   alpha_matrix = sp.alpha,
                                                   model_family = model_family,
                                                   alpha_form = alpha_form,
                                                   lambda_cov_form = lambda_cov_form,
                                                   alpha_cov_form = alpha_cov_form,
                                                   timesteps = 2, # because t1 is the original
                                                   initial_abundances = sp.abund)
            projected.abund <- sub.abund[2,]
            projected.richness <- sum(projected.abund > 0)
            projected.abundance <- sum(projected.abund)
            projected.evenness <- hill.diversity(projected.abund)
            
            pos <- which(pert.sad$year.predicted == (as.numeric(initial.years[i.year]) + 1) &
                           pert.sad$plot == plots[i.plot] & 
                           pert.sad$type == types[i.type] &
                           pert.sad$intensity == i.step)
            
            pert.sad$value[pos[which(pert.sad$metric[pos] == "richness")]] <- 
              projected.richness
            pert.sad$value[pos[which(pert.sad$metric[pos] == "abundance")]] <- 
              projected.abundance
            pert.sad$value[pos[which(pert.sad$metric[pos] == "evenness")]] <- 
              projected.evenness
          }# if !na
        }# for i.step
        
      }else if(types[i.type] == "ia"){
        
        for(i.step in 1:steps){
          
          lambda.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[types[i.type]]][["lambda"]]
          if(sum(is.na(lambda.df)) == 0){
            sp.alpha <- communities[[initial.years[i.year]]][[plots[i.plot]]][[types[i.type]]][["alpha"]][[i.step]]
          }else{
            sp.alpha <- NA
          }
          abund.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[types[i.type]]][["abundances"]]
          
          if(sum(sum(is.na(lambda.df)),sum(is.na(sp.alpha)),sum(is.na(abund.df))) == 0){
            
            sp.lambda <- lambda.df$lambda
            names(sp.lambda) <- lambda.df$sp
            
            sp.abund <- abund.df$abundance
            names(sp.abund) <- abund.df$species
            
            sub.abund <- cxr::abundance_projection(lambda = sp.lambda,
                                                   alpha_matrix = sp.alpha,
                                                   model_family = model_family,
                                                   alpha_form = alpha_form,
                                                   lambda_cov_form = lambda_cov_form,
                                                   alpha_cov_form = alpha_cov_form,
                                                   timesteps = 2, # because t1 is the original
                                                   initial_abundances = sp.abund)
            projected.abund <- sub.abund[2,]
            projected.richness <- sum(projected.abund > 0)
            projected.abundance <- sum(projected.abund)
            projected.evenness <- hill.diversity(projected.abund)
            
            pos <- which(pert.sad$year.predicted == (as.numeric(initial.years[i.year]) + 1) &
                           pert.sad$plot == plots[i.plot] & 
                           pert.sad$type == types[i.type] &
                           pert.sad$intensity == i.step)
            
            pert.sad$value[pos[which(pert.sad$metric[pos] == "richness")]] <- 
              projected.richness
            pert.sad$value[pos[which(pert.sad$metric[pos] == "abundance")]] <- 
              projected.abundance
            pert.sad$value[pos[which(pert.sad$metric[pos] == "evenness")]] <- 
              projected.evenness
            
          }# if !na
        }# for i.step
        
      }else if(types[i.type] == "id"){
        
        for(i.step in 1:steps){
          
          lambda.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[types[i.type]]][["lambda"]]
          if(sum(is.na(lambda.df)) == 0){
            sp.alpha <- communities[[initial.years[i.year]]][[plots[i.plot]]][[types[i.type]]][["alpha"]][[i.step]]
          }else{
            sp.alpha <- NA
          }
          abund.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[types[i.type]]][["abundances"]]
          
          if(sum(sum(is.na(lambda.df)),sum(is.na(sp.alpha)),sum(is.na(abund.df))) == 0){
            
            sp.lambda <- lambda.df$lambda
            names(sp.lambda) <- lambda.df$sp
            
            sp.abund <- abund.df$abundance
            names(sp.abund) <- abund.df$species
            
            sub.abund <- cxr::abundance_projection(lambda = sp.lambda,
                                                   alpha_matrix = sp.alpha,
                                                   model_family = model_family,
                                                   alpha_form = alpha_form,
                                                   lambda_cov_form = lambda_cov_form,
                                                   alpha_cov_form = alpha_cov_form,
                                                   timesteps = 2, # because t1 is the original
                                                   initial_abundances = sp.abund)
            projected.abund <- sub.abund[2,]
            projected.richness <- sum(projected.abund > 0)
            projected.abundance <- sum(projected.abund)
            projected.evenness <- hill.diversity(projected.abund)
            
            pos <- which(pert.sad$year.predicted == (as.numeric(initial.years[i.year]) + 1) &
                           pert.sad$plot == plots[i.plot] & 
                           pert.sad$type == types[i.type] &
                           pert.sad$intensity == i.step)
            
            pert.sad$value[pos[which(pert.sad$metric[pos] == "richness")]] <- 
              projected.richness
            pert.sad$value[pos[which(pert.sad$metric[pos] == "abundance")]] <- 
              projected.abundance
            pert.sad$value[pos[which(pert.sad$metric[pos] == "evenness")]] <- 
              projected.evenness
          }# if !na
        }# for i.step
        
      }# if-else type
      
    }# for i.type
  }# for i.plot
}# for i.year


# store results -----------------------------------------------------------

write.csv2(pert.sad,file = "results/predicted_SAD_components.csv",
           row.names = FALSE)

