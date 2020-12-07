
# a toy script for seeing the effects of removing niche/fitness diff
# in a given community (matrix + lambda)

library(tidyverse)

source("R/gradient_asymmetry.R")
source("R/gradient_fitness_diff.R")
source("R/gradient_niche_diff.R")
source("R/gradient_strength_dist.R")
source("R/hill_diversity.R")

# some constants ----------------------------------------------------------

metrics <- c("richness","abundance","evenness")

# load the model
source("./R/AP_pm_alpha_pairwise_lambdacov_none_alphacov_none.R")
source("./R/AP_project_alpha_pairwise_lambdacov_none_alphacov_none.R")

model_family <- "AP"
optimization_method <- "bobyqa"
alpha_form <- "pairwise"
lambda_cov_form <- "none"
alpha_cov_form <- "none"

# read data ----------------------------------------------------

lambda <- read.csv2("./results/lambda.csv",stringsAsFactors = FALSE)
alpha.df <- read.csv2("results/alpha.csv",stringsAsFactors = FALSE,
                      row.names = 1)
alpha.matrix <- as.matrix(alpha.df)
alpha.matrix[is.na(alpha.matrix)] <- 0
non.zero <- alpha.matrix[which(alpha.matrix != 0)]

abund.orig <- read.csv2("../Caracoles/data/abundances.csv",
                        header = TRUE,stringsAsFactors = FALSE)
sp.rates <- read.csv2("../Caracoles/data/plant_species_traits.csv",
                      header = TRUE,stringsAsFactors = FALSE)
sp.valid <- sp.rates$species.code[which(!is.na(sp.rates$germination.rate))]

base.abund <- abund.orig %>% 
  filter(species %in% sp.valid & species %in% rownames(alpha.matrix) & 
           individuals > 0) 

# matrix replicates
replicates <- 1 #TODO test
sp <- c("CETE","LYTR","PAIN","POMA","POMO","PUPA","SASO","SPRU")
lambda.obs <- lambda$lambda[lambda$sp %in% sp]
abund.obs <- c(3,3,6,27,1,3,2,1)
alpha.obs <- alpha.matrix[sp,sp]

# gradient steps
steps <- 10
# abundance projection timesteps
timesteps <- 2000
# number of initial species (up to 12)
nsp <- 12
# persistence threshold
persistence.threshold <- 1e-2

# types
types <- c("obs","nd","fd","ia","id")

communities <- list()

# results data structure --------------------------------------------------

pert.matrices <- expand.grid(rep = 1:replicates,
                             type = types,
                             intensity = 1:steps,
                             metric = metrics,
                             value = NA_real_)

# build communities -------------------------------------------------------

for(i.rep in 1:replicates){
  
  communities[[i.rep]] <- list()
  
# build community ---------------------------------------------------------

  # my.names <- paste("sp",1:nsp,sep="")
  # 
  # alpha.obs <- matrix(data = sample(non.zero,length(non.zero)),nrow = nsp,
  #                     dimnames = list(my.names,my.names))
  # lambda.obs <- sample(lambda$lambda,nsp)
  # 
  # # abundances also randomly sampled
  # abund.obs <- sample(base.abund$individuals,nsp)

# perturb commun ----------------------------------------------------------
  
  for(i.type in 1:length(types)){
    
    communities[[i.rep]][[i.type]] <- list()
    
    if(types[i.type] == "obs"){
      
        communities[[i.rep]][[i.type]][[1]] <- abund.obs
        communities[[i.rep]][[i.type]][[2]] <- lambda.obs
        communities[[i.rep]][[i.type]][[3]] <- alpha.obs
 
    }else if(types[i.type] == "nd"){
      
        alpha.nd <- gradient_niche_diff(A = alpha.obs,steps = steps)
        
        communities[[i.rep]][[i.type]][[1]] <- abund.obs
        communities[[i.rep]][[i.type]][[2]] <- lambda.obs
        communities[[i.rep]][[i.type]][[3]] <- alpha.nd        
    
    }else if(types[i.type] == "fd"){
      
        lambda.fd <- gradient_fitness_diff(lambda = lambda.obs,
                                           steps = steps)
        
        communities[[i.rep]][[i.type]][[1]] <- abund.obs
        communities[[i.rep]][[i.type]][[2]] <- lambda.fd
        communities[[i.rep]][[i.type]][[3]] <- alpha.obs
        
    }else if(types[i.type] == "ia"){

        alpha.ia <- gradient_asymmetry(A = alpha.obs,steps = steps)
        
        communities[[i.rep]][[i.type]][[1]] <- abund.obs
        communities[[i.rep]][[i.type]][[2]] <- lambda.obs
        communities[[i.rep]][[i.type]][[3]] <- alpha.ia
      
    }else if(types[i.type] == "id"){
      
        alpha.id <- gradient_strength_dist(A = alpha.obs,steps = steps)
        
        communities[[i.rep]][[i.type]][[1]] <- abund.obs
        communities[[i.rep]][[i.type]][[2]] <- lambda.obs
        communities[[i.rep]][[i.type]][[3]] <- alpha.id
    }
    
    names(communities[[i.rep]][[i.type]]) <- c("abundances",
                                               "lambda","alpha")
    
  }# for each type
  names(communities[[i.rep]]) <- types
}# for i.rep

# extract SAD components --------------------------------------------------

for(i.rep in 1:replicates){
      for(i.type in 1:length(types)){
        
        if(types[i.type] == "obs"){
          
          lambda.df <- communities[[i.rep]][[types[i.type]]][["lambda"]]
          sp.alpha <- communities[[i.rep]][[types[i.type]]][["alpha"]]
          abund.df <- communities[[i.rep]][[types[i.type]]][["abundances"]]
          
          if(sum(sum(is.na(lambda.df)),sum(is.na(sp.alpha)),sum(is.na(abund.df))) == 0){
            
            sp.lambda <- lambda.df
            names(sp.lambda) <- colnames(sp.alpha)
            
            sp.abund <- abund.df
            names(sp.abund) <- colnames(sp.alpha)
            
            sub.abund <- cxr::abundance_projection(lambda = sp.lambda,
                                                   alpha_matrix = sp.alpha,
                                                   model_family = model_family,
                                                   alpha_form = alpha_form,
                                                   lambda_cov_form = lambda_cov_form,
                                                   alpha_cov_form = alpha_cov_form,
                                                   timesteps = timesteps,
                                                   initial_abundances = sp.abund)
            projected.abund <- sub.abund[timesteps,]
            projected.richness <- sum(projected.abund > persistence.threshold)
            projected.abundance <- sum(projected.abund)
            projected.evenness <- hill.diversity(projected.abund)
            
            pos <- which(pert.matrices$rep == i.rep &
                           pert.matrices$type == "obs" )
            
            pert.matrices$value[pos[which(pert.matrices$metric[pos] == "richness")]] <- 
              projected.richness
            pert.matrices$value[pos[which(pert.matrices$metric[pos] == "abundance")]] <- 
              projected.abundance
            pert.matrices$value[pos[which(pert.matrices$metric[pos] == "evenness")]] <- 
              projected.evenness
          }# if !na
        }else if(types[i.type] == "nd"){
          
          for(i.step in 1:steps){
            
            lambda.df <- communities[[i.rep]][[types[i.type]]][["lambda"]]
            if(sum(is.na(lambda.df)) == 0){
              sp.alpha <- communities[[i.rep]][[types[i.type]]][["alpha"]][[i.step]]
            }else{
              sp.alpha <- NA
            }
            abund.df <- communities[[i.rep]][[types[i.type]]][["abundances"]]
            
            if(sum(sum(is.na(lambda.df)),sum(is.na(sp.alpha)),sum(is.na(abund.df))) == 0){
              
              sp.lambda <- lambda.df
              names(sp.lambda) <- colnames(sp.alpha)
              
              sp.abund <- abund.df
              names(sp.abund) <- colnames(sp.alpha)
              
              sub.abund <- cxr::abundance_projection(lambda = sp.lambda,
                                                     alpha_matrix = sp.alpha,
                                                     model_family = model_family,
                                                     alpha_form = alpha_form,
                                                     lambda_cov_form = lambda_cov_form,
                                                     alpha_cov_form = alpha_cov_form,
                                                     timesteps = timesteps, # because t1 is the original
                                                     initial_abundances = sp.abund)
              projected.abund <- sub.abund[timesteps,]
              projected.richness <- sum(projected.abund > persistence.threshold)
              projected.abundance <- sum(projected.abund)
              projected.evenness <- hill.diversity(projected.abund)
              
              pos <- which(pert.matrices$rep == i.rep &
                             pert.matrices$type == types[i.type] &
                             pert.matrices$intensity == i.step)
              
              pert.matrices$value[pos[which(pert.matrices$metric[pos] == "richness")]] <- 
                projected.richness
              pert.matrices$value[pos[which(pert.matrices$metric[pos] == "abundance")]] <- 
                projected.abundance
              pert.matrices$value[pos[which(pert.matrices$metric[pos] == "evenness")]] <- 
                projected.evenness
            }# if !na
          }# for i.step
        }else if(types[i.type] == "fd"){
          
          for(i.step in 1:steps){
            
            sp.alpha <- communities[[i.rep]][[types[i.type]]][["alpha"]]
            abund.df <- communities[[i.rep]][[types[i.type]]][["abundances"]]
            if(sum(is.na(sp.alpha)) == 0){
              lambda.df <- communities[[i.rep]][[types[i.type]]][["lambda"]][[i.step]]
            }else{
              lambda.df <- NA
            }
            
            if(sum(sum(is.na(lambda.df)),sum(is.na(sp.alpha)),sum(is.na(abund.df))) == 0){
              
              sp.lambda <- lambda.df
              names(sp.lambda) <- colnames(sp.alpha)
              
              sp.abund <- abund.df
              names(sp.abund) <- colnames(sp.alpha)
              
              sub.abund <- cxr::abundance_projection(lambda = sp.lambda,
                                                     alpha_matrix = sp.alpha,
                                                     model_family = model_family,
                                                     alpha_form = alpha_form,
                                                     lambda_cov_form = lambda_cov_form,
                                                     alpha_cov_form = alpha_cov_form,
                                                     timesteps = timesteps, # because t1 is the original
                                                     initial_abundances = sp.abund)
              projected.abund <- sub.abund[timesteps,]
              projected.richness <- sum(projected.abund > persistence.threshold)
              projected.abundance <- sum(projected.abund)
              projected.evenness <- hill.diversity(projected.abund)
              
              pos <- which(pert.matrices$rep == i.rep &
                             pert.matrices$type == types[i.type] &
                             pert.matrices$intensity == i.step)
              
              pert.matrices$value[pos[which(pert.matrices$metric[pos] == "richness")]] <- 
                projected.richness
              pert.matrices$value[pos[which(pert.matrices$metric[pos] == "abundance")]] <- 
                projected.abundance
              pert.matrices$value[pos[which(pert.matrices$metric[pos] == "evenness")]] <- 
                projected.evenness
            }# if !na
          }# for i.step
          
        }else if(types[i.type] == "ia"){
          
          for(i.step in 1:steps){
            
            lambda.df <- communities[[i.rep]][[types[i.type]]][["lambda"]]
            if(sum(is.na(lambda.df)) == 0){
              sp.alpha <- communities[[i.rep]][[types[i.type]]][["alpha"]][[i.step]]
            }else{
              sp.alpha <- NA
            }
            abund.df <- communities[[i.rep]][[types[i.type]]][["abundances"]]
            
            if(sum(sum(is.na(lambda.df)),sum(is.na(sp.alpha)),sum(is.na(abund.df))) == 0){
              
              sp.lambda <- lambda.df
              names(sp.lambda) <- colnames(sp.alpha)
              
              sp.abund <- abund.df
              names(sp.abund) <- colnames(sp.alpha)
              
              sub.abund <- cxr::abundance_projection(lambda = sp.lambda,
                                                     alpha_matrix = sp.alpha,
                                                     model_family = model_family,
                                                     alpha_form = alpha_form,
                                                     lambda_cov_form = lambda_cov_form,
                                                     alpha_cov_form = alpha_cov_form,
                                                     timesteps = timesteps, # because t1 is the original
                                                     initial_abundances = sp.abund)
              projected.abund <- sub.abund[timesteps,]
              projected.richness <- sum(projected.abund > persistence.threshold)
              projected.abundance <- sum(projected.abund)
              projected.evenness <- hill.diversity(projected.abund)
              
              pos <- which(pert.matrices$rep == i.rep &
                             pert.matrices$type == types[i.type] &
                             pert.matrices$intensity == i.step)
              
              pert.matrices$value[pos[which(pert.matrices$metric[pos] == "richness")]] <- 
                projected.richness
              pert.matrices$value[pos[which(pert.matrices$metric[pos] == "abundance")]] <- 
                projected.abundance
              pert.matrices$value[pos[which(pert.matrices$metric[pos] == "evenness")]] <- 
                projected.evenness
              
            }# if !na
          }# for i.step
          
        }else if(types[i.type] == "id"){
          
          for(i.step in 1:steps){
            
            lambda.df <- communities[[i.rep]][[types[i.type]]][["lambda"]]
            if(sum(is.na(lambda.df)) == 0){
              sp.alpha <- communities[[i.rep]][[types[i.type]]][["alpha"]][[i.step]]
            }else{
              sp.alpha <- NA
            }
            abund.df <- communities[[i.rep]][[types[i.type]]][["abundances"]]
            
            if(sum(sum(is.na(lambda.df)),sum(is.na(sp.alpha)),sum(is.na(abund.df))) == 0){
              
              sp.lambda <- lambda.df
              names(sp.lambda) <- colnames(sp.alpha)
              
              sp.abund <- abund.df
              names(sp.abund) <- colnames(sp.alpha)
              
              sub.abund <- cxr::abundance_projection(lambda = sp.lambda,
                                                     alpha_matrix = sp.alpha,
                                                     model_family = model_family,
                                                     alpha_form = alpha_form,
                                                     lambda_cov_form = lambda_cov_form,
                                                     alpha_cov_form = alpha_cov_form,
                                                     timesteps = timesteps, # because t1 is the original
                                                     initial_abundances = sp.abund)
              projected.abund <- sub.abund[timesteps,]
              projected.richness <- sum(projected.abund > persistence.threshold)
              projected.abundance <- sum(projected.abund)
              projected.evenness <- hill.diversity(projected.abund)
              
              pos <- which(pert.matrices$rep == i.rep &
                             pert.matrices$type == types[i.type] &
                             pert.matrices$intensity == i.step)
              
              pert.matrices$value[pos[which(pert.matrices$metric[pos] == "richness")]] <- 
                projected.richness
              pert.matrices$value[pos[which(pert.matrices$metric[pos] == "abundance")]] <- 
                projected.abundance
              pert.matrices$value[pos[which(pert.matrices$metric[pos] == "evenness")]] <- 
                projected.evenness
            }# if !na
          }# for i.step
          
        }# if-else type
        
      }# for i.type
}# for i.rep

# plot --------------------------------------------------------------------

obs.pred.plot <- ggplot(pert.matrices, aes(x = type, y = value)) + 
  geom_point(aes(fill = type), size = 3, shape = 21, 
             position = position_jitterdodge()) +
  # geom_boxplot(aes(fill = type),alpha = 0.60) +
  scale_fill_OkabeIto() +
  facet_grid(metric~intensity,scales = "free_y") +
  labs(x = NULL,y = NULL)+
  theme_bw()+
  scale_x_discrete(breaks=NULL)+
  theme(strip.background = element_blank())+
  NULL
 # obs.pred.plot

# library(vegan)
# tt <- c(1,2,3,1,1,1)
# H <- diversity(tt)
# S <- specnumber(tt) ## rowSums(BCI > 0) does the same...
# J <- H/log(S)
# J
# hill.diversity(tt)
