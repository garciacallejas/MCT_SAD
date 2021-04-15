
# perturb the stabilized communities

library(tidyverse)

source("R/gradient_asymmetry.R")
source("R/gradient_fitness_diff.R")
source("R/gradient_niche_diff.R")
source("R/gradient_strength_dist.R")
source("R/gradient_diag_dominance.R")

# some constants ----------------------------------------------------------

# model_family <- c("BH","LW","RK")
model_family <- "LV"
steps <- 10
types <- c("obs","ia","id","dd")
communities <- list()

# read data ---------------------------------------------------------------

load("results/communities_subplot_stabilized.Rdata")

years <- names(stable_communities)
plots <- 1:length(stable_communities[[1]])
subplots <- names(stable_communities[[1]][[1]])

# generate perturbed values -----------------------------------------------

for(i.year in 1:length(years)){
  
  communities[[i.year]] <- list()
  
  for(i.plot in 1:length(plots)){
    
    communities[[i.year]][[i.plot]] <- list()
    
    for(i.sub in 1:length(subplots)){
      
      if(!is.na(stable_communities[[years[i.year]]][[i.plot]][[subplots[i.sub]]][["corrected_alpha"]])){
        
        abund.obs <- stable_communities[[years[i.year]]][[i.plot]][[subplots[i.sub]]][["abundances"]]
        corrected_r <- stable_communities[[years[i.year]]][[i.plot]][[subplots[i.sub]]][["corrected_r"]]
        corrected_alpha <- stable_communities[[years[i.year]]][[i.plot]][[subplots[i.sub]]][["corrected_alpha"]]
        
        communities[[i.year]][[i.plot]][[i.sub]] <- list()
        
        for(i.type in 1:length(types)){
          
          communities[[i.year]][[i.plot]][[i.sub]][[i.type]] <- list()
          
          if(types[i.type] == "obs"){
            
            communities[[i.year]][[i.plot]][[i.sub]][[i.type]][[1]] <- abund.obs
            communities[[i.year]][[i.plot]][[i.sub]][[i.type]][[2]] <- corrected_r
            communities[[i.year]][[i.plot]][[i.sub]][[i.type]][[3]] <- corrected_alpha
            
          }else if(types[i.type] == "nd"){
            
            alpha.nd <- gradient_niche_diff(A = corrected_alpha,steps = steps)
            
            communities[[i.year]][[i.plot]][[i.sub]][[i.type]][[1]] <- abund.obs
            communities[[i.year]][[i.plot]][[i.sub]][[i.type]][[2]] <- corrected_r
            communities[[i.year]][[i.plot]][[i.sub]][[i.type]][[3]] <- alpha.nd        
            
          }else if(types[i.type] == "fd"){
            
            r.fd <- gradient_fitness_diff(lambda = corrected_r$rfit,
                                          steps = steps)
            
            r.fd.list <- list()
            for(i.step in 1:length(lambda.fd)){
              r.fd.list[[i.step]] <- data.frame(sp = corrected_r$sp,
                                                r = r.fd[[i.step]])
            }
            
            communities[[i.year]][[i.plot]][[i.sub]][[i.type]][[1]] <- abund.obs
            communities[[i.year]][[i.plot]][[i.sub]][[i.type]][[2]] <- r.fd.list
            communities[[i.year]][[i.plot]][[i.sub]][[i.type]][[3]] <- corrected_alpha
            
          }else if(types[i.type] == "ia"){
            
            alpha.ia <- gradient_asymmetry(A = corrected_alpha,steps = steps)
            
            communities[[i.year]][[i.plot]][[i.sub]][[i.type]][[1]] <- abund.obs
            communities[[i.year]][[i.plot]][[i.sub]][[i.type]][[2]] <- corrected_r
            communities[[i.year]][[i.plot]][[i.sub]][[i.type]][[3]] <- alpha.ia
            
          }else if(types[i.type] == "id"){
            
            alpha.id <- gradient_strength_dist(A = corrected_alpha,steps = steps)
            
            communities[[i.year]][[i.plot]][[i.sub]][[i.type]][[1]] <- abund.obs
            communities[[i.year]][[i.plot]][[i.sub]][[i.type]][[2]] <- corrected_r
            communities[[i.year]][[i.plot]][[i.sub]][[i.type]][[3]] <- alpha.id
            
          }else if(types[i.type] == "dd"){
            
            alpha.id <- gradient_diag_dominance(A = corrected_alpha,steps = steps)
            
            communities[[i.year]][[i.plot]][[i.sub]][[i.type]][[1]] <- abund.obs
            communities[[i.year]][[i.plot]][[i.sub]][[i.type]][[2]] <- corrected_r
            communities[[i.year]][[i.plot]][[i.sub]][[i.type]][[3]] <- alpha.id
            
          }# if-else types
          
          names(communities[[i.year]][[i.plot]][[i.sub]][[i.type]]) <- c("abundances",
                                                                         "r","alpha")
          
        }# for each type
        names(communities[[i.year]][[i.plot]][[i.sub]]) <- types
        
      }else{
        
        communities[[i.year]][[i.plot]][[i.sub]] <- NA
        
      }# if-else valid corrected community
    }# for i.sub
    names(communities[[i.year]][[i.plot]]) <- subplots
  }# for i.plot
}# for i.year
names(communities) <- years

# write to disk -----------------------------------------------------------

save(communities,file = "results/communities_subplot_perturbed.Rdata")

