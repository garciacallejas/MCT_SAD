
library(tidyverse)
library(cxr)
source("R/hill_diversity.R")

# read data ---------------------------------------------------------------

load("results/communities_subplot.Rdata")

# species names
lambda.orig <- read.csv2("results/lambda_RK.csv",stringsAsFactors = FALSE)

# some constants ----------------------------------------------------------

timesteps <- 20
persistence.threshold <- 1e-3

years <- names(communities)

# I can "only" predict 15->16, 16->17; I lack 2018 to predict 2019.
years.predicted <- c("2016","2017")

# likewise, start by projecting abundances from 2015 and 2016
initial.years <- c("2015","2016")

species <- sort(unique(lambda.orig$sp))
plots <- 1:length(communities[[1]])

# types <- names(communities[[1]][[1]][[1]])
types <- c("obs","ia","id","dd")

subplots <- names(communities[[1]][[1]])
steps <- length(communities[[1]][[1]][[1]][["nd"]][["alpha"]])

metrics <- c("richness","abundance","evenness")

# load the model
# source("./R/AP_pm_alpha_pairwise_lambdacov_none_alphacov_none.R")
# source("./R/AP_project_alpha_pairwise_lambdacov_none_alphacov_none.R")

model_family <- "RK"
optimization_method <- "bobyqa"
alpha_form <- "pairwise"
lambda_cov_form <- "none"
alpha_cov_form <- "none"

# results data structures -------------------------------------------------

# predicted abundance of each sp at the plot level
# pred.plot <- expand.grid(year.predicted = years.predicted,
#                          plot = plots,
#                          type = types, 
#                          intensity = 1:steps,
#                          timesteps = 1:timesteps,
#                          species = species,
#                          abund = 0)

# list keeping the metrics for each plot/year
all.plots.list <- list()
# count <- 1
timesteps.to.keep <- c(2,5,10,15,20)

# project abundances for every community ----------------------------------
# and store them at the plot level

for(i.year in 1:length(initial.years)){
  for(i.plot in 1:length(plots)){
    
    # list keeping all raw abundance projections
    plot.list <- list()
    count <- 1
    
    for(i.sub in 1:length(subplots)){
    for(i.type in 1:length(types)){
      
      if(types[i.type] == "obs"){
        
        lambda.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["lambda"]]
        sp.alpha <- communities[[initial.years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["alpha"]]
        abund.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["abundances"]]
        
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
                                                 timesteps = timesteps,
                                                 initial_abundances = sp.abund)
          
          sub.abund.df <- as.data.frame(sub.abund)
          sub.abund.df$year.predicted <- (as.numeric(initial.years[i.year]) + 1)
          sub.abund.df$plot <- plots[i.plot]
          sub.abund.df$type <- types[i.type]
          sub.abund.df$timestep <- 1:nrow(sub.abund.df)
          
          sub.abund.df <- pivot_longer(sub.abund.df,1:ncol(sub.abund),
                                       names_to = "species",
                                       values_to = "abund")
          
          for(i.int in 1:steps){
            sub.abund.df$intensity <- i.int
            sub.abund.df <- sub.abund.df[,c("year.predicted","plot","type",
                                            "intensity","timestep","species","abund")]
            plot.list[[count]] <- sub.abund.df
            count <- count + 1
          }
          
        }# if !na
      }else if(types[i.type] == "nd"){
        
        for(i.step in 1:steps){
          
          lambda.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["lambda"]]
          if(sum(is.na(lambda.df)) == 0){
            sp.alpha <- communities[[initial.years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["alpha"]][[i.step]]
          }else{
            sp.alpha <- NA
          }
          abund.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["abundances"]]
          
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
                                                   timesteps = timesteps, # because t1 is the original
                                                   initial_abundances = sp.abund)
            sub.abund.df <- as.data.frame(sub.abund)
            sub.abund.df$year.predicted <- (as.numeric(initial.years[i.year]) + 1)
            sub.abund.df$plot <- plots[i.plot]
            sub.abund.df$type <- types[i.type]
            sub.abund.df$intensity <- i.step
            sub.abund.df$timestep <- 1:nrow(sub.abund.df)
            
            sub.abund.df <- pivot_longer(sub.abund.df,1:ncol(sub.abund),
                                         names_to = "species",
                                         values_to = "abund")
            
            plot.list[[count]] <- sub.abund.df
            count <- count + 1
            
          }# if !na
        }# for i.step
        
      }else if(types[i.type] == "fd"){
        
        for(i.step in 1:steps){
          
          sp.alpha <- communities[[initial.years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["alpha"]]
          abund.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["abundances"]]
          if(sum(is.na(sp.alpha)) == 0){
            lambda.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["lambda"]][[i.step]]
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
                                                   timesteps = timesteps, # because t1 is the original
                                                   initial_abundances = sp.abund)
            sub.abund.df <- as.data.frame(sub.abund)
            sub.abund.df$year.predicted <- (as.numeric(initial.years[i.year]) + 1)
            sub.abund.df$plot <- plots[i.plot]
            sub.abund.df$type <- types[i.type]
            sub.abund.df$intensity <- i.step
            sub.abund.df$timestep <- 1:nrow(sub.abund.df)
            
            sub.abund.df <- pivot_longer(sub.abund.df,1:ncol(sub.abund),
                                         names_to = "species",
                                         values_to = "abund")
            
            plot.list[[count]] <- sub.abund.df
            count <- count + 1
            
          }# if !na
        }# for i.step
        
      }else if(types[i.type] == "ia"){
        
        for(i.step in 1:steps){
          
          lambda.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["lambda"]]
          if(sum(is.na(lambda.df)) == 0){
            sp.alpha <- communities[[initial.years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["alpha"]][[i.step]]
          }else{
            sp.alpha <- NA
          }
          abund.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["abundances"]]
          
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
                                                   timesteps = timesteps, # because t1 is the original
                                                   initial_abundances = sp.abund)
            
            sub.abund.df <- as.data.frame(sub.abund)
            sub.abund.df$year.predicted <- (as.numeric(initial.years[i.year]) + 1)
            sub.abund.df$plot <- plots[i.plot]
            sub.abund.df$type <- types[i.type]
            sub.abund.df$intensity <- i.step
            sub.abund.df$timestep <- 1:nrow(sub.abund.df)
            
            sub.abund.df <- pivot_longer(sub.abund.df,1:ncol(sub.abund),
                                         names_to = "species",
                                         values_to = "abund")
            
            plot.list[[count]] <- sub.abund.df
            count <- count + 1
            
          }# if !na
        }# for i.step
        
      }else if(types[i.type] == "id"){
        
        for(i.step in 1:steps){
          
          lambda.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["lambda"]]
          if(sum(is.na(lambda.df)) == 0){
            sp.alpha <- communities[[initial.years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["alpha"]][[i.step]]
          }else{
            sp.alpha <- NA
          }
          abund.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["abundances"]]
          
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
                                                   timesteps = timesteps, # because t1 is the original
                                                   initial_abundances = sp.abund)
            
            sub.abund.df <- as.data.frame(sub.abund)
            sub.abund.df$year.predicted <- (as.numeric(initial.years[i.year]) + 1)
            sub.abund.df$plot <- plots[i.plot]
            sub.abund.df$type <- types[i.type]
            sub.abund.df$intensity <- i.step
            sub.abund.df$timestep <- 1:nrow(sub.abund.df)
            
            sub.abund.df <- pivot_longer(sub.abund.df,1:ncol(sub.abund),
                                         names_to = "species",
                                         values_to = "abund")
            
            plot.list[[count]] <- sub.abund.df
            count <- count + 1
            
          }# if !na
        }# for i.step
        
      }else if(types[i.type] == "dd"){
        
        for(i.step in 1:steps){
          
          lambda.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["lambda"]]
          if(sum(is.na(lambda.df)) == 0){
            sp.alpha <- communities[[initial.years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["alpha"]][[i.step]]
          }else{
            sp.alpha <- NA
          }
          abund.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["abundances"]]
          
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
                                                   timesteps = timesteps, # because t1 is the original
                                                   initial_abundances = sp.abund)
            
            sub.abund.df <- as.data.frame(sub.abund)
            sub.abund.df$year.predicted <- (as.numeric(initial.years[i.year]) + 1)
            sub.abund.df$plot <- plots[i.plot]
            sub.abund.df$type <- types[i.type]
            sub.abund.df$intensity <- i.step
            sub.abund.df$timestep <- 1:nrow(sub.abund.df)
            
            sub.abund.df <- pivot_longer(sub.abund.df,1:ncol(sub.abund),
                                         names_to = "species",
                                         values_to = "abund")
            
            plot.list[[count]] <- sub.abund.df
            count <- count + 1
            
          }# if !na
        }# for i.step
        
      }# if-else type
      
    }# for i.type
    }# for i.sub
    
    if(length(plot.list)>0){
      plot.data <- bind_rows(plot.list) %>% 
        filter(abund > persistence.threshold) %>%
        filter(timestep %in% timesteps.to.keep) %>%
        group_by(year.predicted,plot,type,intensity,timestep,species) %>%
        summarise(abund.plot = sum(abund))
      
      plot.metrics <- plot.data %>%
        group_by(year.predicted,plot,type,intensity,timestep) %>%
        summarise(abundance = sum(abund.plot),
                  richness = n(),
                  evenness = hill.diversity(abund.plot))
      
      all.plots.list[[length(all.plots.list)+1]] <- plot.metrics
    }# if plot ok
    
  }# for i.plot
}# for i.year

# metrics at the plot level -----------------------------------------------

pred.plot <- bind_rows(all.plots.list)

# store results -----------------------------------------------------------

pert.sad <- pivot_longer(pred.plot,cols = abundance:evenness,
                         names_to = "metric",
                         values_to = "value")

write.csv2(pert.sad,file = "results/predicted_SAD_components_subplot_aggregated.csv",
           row.names = FALSE)

