
library(tidyverse)
library(cxr)
source("R/hill_diversity.R")

# read data ---------------------------------------------------------------

load("results/communities_subplot.Rdata")

# species names
# lambda.orig <- read.csv2("results/lambda.csv",stringsAsFactors = FALSE)

# load("results/SIMPLEX_fit_co3_AP.Rdata")
load("results/SIMPLEX_fit_AP.Rdata")

salinity <- read.csv2("../Caracoles/data/salinity.csv")

# some constants ----------------------------------------------------------

timesteps <- 3
persistence.threshold <- 1e-3

years <- names(communities)

# likewise, start by projecting abundances from all years except the most recent
# so that I can compare projected and observed
# in addition, note that all projections are from t to t+1, for now
initial.years <- c("2015","2016","2018")#years[1:(length(years)-1)]

species <- sort(unique(pred[[1]]$sp))
plots <- 1:length(communities[[1]])

# types <- names(communities[[1]][[1]][[1]])
types <- c("obs")

subplots <- names(communities[[1]][[1]])
# steps <- length(communities[[1]][[1]][[1]][["nd"]][["alpha"]])

# metrics <- c("richness","abundance","evenness")

# load the model
source("./R/AP_pm_alpha_pairwise_lambdacov_none_alphacov_none.R")
source("./R/AP_project_alpha_pairwise_lambdacov_none_alphacov_none.R")
source("./R/AP_pm_alpha_pairwise_lambdacov_global_alphacov_global.R")
source("./R/AP_project_alpha_pairwise_lambdacov_global_alphacov_global.R")
source("./R/AP_pm_alpha_pairwise_lambdacov_global_alphacov_pairwise.R")
source("./R/AP_project_alpha_pairwise_lambdacov_global_alphacov_pairwise.R")

model_family <- "AP"
optimization_method <- "bobyqa"
alpha_form <- "pairwise"
lambda_cov_form <- "none"
alpha_cov_form <- "none"
# lambda_cov_form <- "global"
# alpha_cov_form <- "global"

# results data structures -------------------------------------------------

abund.projected <- list()

# list keeping the metrics for each plot/year
all.plots.list <- list()
# count <- 1
timesteps.to.keep <- c(2)

# project abundances for every community ----------------------------------
# and store them at the plot level

# i.year <- i.plot <- i.sub <- i.type <- 1

for(i.year in 1:length(initial.years)){
  for(i.plot in 1:length(plots)){
    
    # list keeping all raw abundance projections
    plot.list <- list()
    count <- 1
    
    for(i.sub in 1:length(subplots)){
    for(i.type in 1:length(types)){
      
      if(types[i.type] == "obs"){
        
        abund.df <- communities[[initial.years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["abundances"]]
        
        if(is.data.frame(abund.df)){
          lambda.df <- arrange(pred[[1]],sp)
          lambda.df <- subset(lambda.df,sp %in% abund.df$species)
          
          sp.alpha <- pred[[2]][abund.df$species,abund.df$species]
          
          # uncomment if necessary
          # lambda_cov <- matrix(data = pred[[3]][abund.df$species,"co3"],
          #                      ncol = 1, 
          #                      dimnames = list(abund.df$species,"co3"))
          # 
          # alpha_cov <- list(co3 = pred[[4]][abund.df$species,abund.df$species])
          
          # values for the projected abundances
          # subplot.salinity <- subset(salinity, plot == plots[i.plot] &
          #                              year %in% c((as.numeric(initial.years[i.year])),
          #                                          (as.numeric(initial.years[i.year]))+1,
          #                                          (as.numeric(initial.years[i.year]))+2) &
          #                              subplot == subplots[i.sub])
          # subplot.co3 <- data.frame(co3 = subplot.salinity[,c("co3")])
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
                                                 covariates = NULL,
                                                 # lambda_cov = lambda_cov,
                                                 # alpha_cov = alpha_cov,
                                                 # covariates = subplot.co3,
                                                 timesteps = timesteps,
                                                 initial_abundances = sp.abund)
          
          sub.abund.df <- data.frame(species = names(sp.abund), 
                                     projected.abundance = sub.abund[2,])
          
          # add observed abundance
          ny <- as.character(as.numeric(initial.years[i.year])+1)
          obs.abund <- communities[[ny]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["abundances"]]
          
          if(is.data.frame(obs.abund)){
            sub.abund.df <- left_join(sub.abund.df,obs.abund)
            names(sub.abund.df)[3] <- "observed.abundance"
            sub.abund.df$observed.abundance[which(is.na(sub.abund.df$observed.abundance))] <- 0
          }else{
            sub.abund.df$observed.abundance <- 0
          }
          
          sub.abund.df$year.predicted <- (as.numeric(initial.years[i.year]) + 1)
          sub.abund.df$plot <- plots[i.plot]
          sub.abund.df$subplot <- subplots[i.sub]
          # sub.abund.df$type <- types[i.type]
          # sub.abund.df$timestep <- 1
          
          abund.projected[[length(abund.projected)+1]] <- sub.abund.df
          
          # sub.abund.df <- pivot_longer(sub.abund.df,1:ncol(sub.abund),
          #                              names_to = "species",
          #                              values_to = "abund")
          
          # for(i.int in 1:steps){
          #   sub.abund.df$intensity <- i.int
          #   sub.abund.df <- sub.abund.df[,c("year.predicted","plot","type",
          #                                   "intensity","timestep","species","abund")]
          #   plot.list[[count]] <- sub.abund.df
          #   count <- count + 1
          # }
          
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
    
    # if(length(plot.list)>0){
    #   plot.data <- bind_rows(plot.list) %>% 
    #     filter(abund > persistence.threshold) %>%
    #     filter(timestep %in% timesteps.to.keep) %>%
    #     group_by(year.predicted,plot,type,intensity,timestep,species) %>%
    #     summarise(abund.plot = sum(abund))
    #   
    #   plot.metrics <- plot.data %>%
    #     group_by(year.predicted,plot,type,intensity,timestep) %>%
    #     summarise(abundance = sum(abund.plot),
    #               richness = n(),
    #               evenness = hill.diversity(abund.plot))
    #   
    #   all.plots.list[[length(plot.list)+1]] <- plot.metrics
    # }# if plot ok
    
  }# for i.plot
}# for i.year

# temp - for simplex
pred.subplot <- bind_rows(abund.projected)

write.csv2(pred.subplot,file = "results/simplex/predicted_abundances_AP.csv",row.names = FALSE)

# calculate error ---------------------------------------------------------

pred.no.cov <- read.csv2("results/simplex/predicted_abundances_AP.csv")
pred.cov <- read.csv2("results/simplex/predicted_abundances_AP_co3.csv")

pred.no.cov$covariates <- "none"
pred.cov$covariates <- "co3"

pred.all <- bind_rows(pred.no.cov,pred.cov)

rse <- function(obs,pred){
  num <- sum((pred-obs)^2)
  den <- sum((mean(obs)-obs)^2)
  num/den
}

rmse <- function(obs,pred){
  num <- sum((pred-obs)^2)
  return(sqrt(num/length(obs)))
}

all.sp <- sort(unique(pred.all$species))

pred.clean <- pred.all %>% filter(!is.na(observed.abundance)) %>%
  group_by(species, covariates, year.predicted) %>% 
  summarise(rse = rse(obs = observed.abundance,pred = projected.abundance),
            rmse = rmse(obs = observed.abundance,pred = projected.abundance))

write.csv2(pred.clean,file = "results/simplex/prediction_errors_per_sp_year.csv")

# check
rse.plot <- pred.clean %>% filter(!is.infinite(rse) & !is.na(rse)) %>%
  group_by(species,covariates) %>% summarise(mean.rse = mean(rse),
                                             mean.rmse = mean(rmse)) %>%
  ggplot(aes(y = log(mean.rse), x = species)) +
  geom_point() +
  facet_grid(covariates~., scales = "free_y")+
  NULL

#rse table
pred.clean %>% filter(!is.infinite(rse) & !is.na(rse)) %>%
  group_by(covariates) %>%
  summarise(median.rse = median(rse),
            median.rmse = median(rmse))


ggsave("results/simplex/rse_plot.pdf",plot = rse.plot,device = cairo_pdf,
              width = 6,height = 4,dpi = 300)

# metrics at the plot level -----------------------------------------------

# pred.plot <- bind_rows(all.plots.list)

# store results -----------------------------------------------------------

# pert.sad <- pivot_longer(pred.plot,cols = abundance:evenness,
#                          names_to = "metric",
#                          values_to = "value")

# write.csv2(pert.sad,file = "results/predicted_SAD_components_subplot_aggregated_v2.csv",
#            row.names = FALSE)

