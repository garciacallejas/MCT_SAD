
# generate perturbed alphas/lambdas

library(tidyverse)

source("R/gradient_asymmetry.R")
source("R/gradient_fitness_diff.R")
source("R/gradient_niche_diff.R")
source("R/gradient_strength_dist.R")
source("R/gradient_diag_dominance.R")

# some constants ----------------------------------------------------------

model_family <- c("BH","LW","RK")
steps <- 10
types <- c("obs","nd","fd","ia","id","dd")
communities <- list()

# read data ---------------------------------------------------------------

abund <- read.csv2("results/abund_filtered.csv",
                   header = TRUE,stringsAsFactors = FALSE)

years <- sort(unique(abund$year))
plots <- sort(unique(abund$plot))
subplots <- sort(unique(abund$subplot))

# generate perturbed values -----------------------------------------------

for(i.model in 1:length(model_family)){
  
  # specific model coefs
  lambda <- read.csv2(paste("./results/lambda_",model_family[i.model],".csv",sep=""),stringsAsFactors = FALSE)
  alpha.df <- read.csv2(paste("results/alpha_",model_family[i.model],".csv",sep=""),stringsAsFactors = FALSE,
                        row.names = 1)
  alpha.matrix <- as.matrix(alpha.df)
  
  # it should be equal for all models, but just in case
  base.abund <- abund %>% 
    filter(species %in% rownames(alpha.matrix)) 
  
  communities[[i.model]] <- list()
  
  for(i.year in 1:length(years)){
    
    communities[[i.model]][[i.year]] <- list()
    
    for(i.plot in 1:length(plots)){
      
      communities[[i.model]][[i.year]][[i.plot]] <- list()
      
      for(i.sub in 1:length(subplots)){
        
        communities[[i.model]][[i.year]][[i.plot]][[i.sub]] <- list()
        
        # subset present species
        present.sp <- sort(unique(base.abund$species[base.abund$year == years[i.year] &
                                                       base.abund$plot == plots[i.plot] &
                                                       base.abund$subplot == subplots[i.sub] &
                                                       base.abund$individuals > 0]))
        
        # 3 - sum observed abundances
        # surely not needed at the subplot level...
        abund.obs <- base.abund %>% 
          filter(year == years[i.year] &
                   plot == plots[i.plot] &
                   subplot == subplots[i.sub] &
                   species %in% present.sp) %>%
          group_by(species) %>%
          summarise(abundance = sum(individuals))
        
        lambda.obs <- lambda[lambda$sp %in% present.sp,]
        lambda.obs <- arrange(lambda.obs,sp)
        
        alpha.obs <- alpha.matrix[present.sp,present.sp]
        alpha.obs[which(is.na(alpha.obs))] <- 0
        
        for(i.type in 1:length(types)){
          
          communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]] <- list()
          
          if(types[i.type] == "obs"){
            
            if(length(present.sp)>1){
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[1]] <- abund.obs
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[2]] <- lambda.obs
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[3]] <- alpha.obs
            }else{
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[1]] <- NA
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[2]] <- NA
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[3]] <- NA
            }
          }else if(types[i.type] == "nd"){
            
            if(length(present.sp)>1){
              alpha.nd <- gradient_niche_diff(A = alpha.obs,steps = steps)
              
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[1]] <- abund.obs
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[2]] <- lambda.obs
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[3]] <- alpha.nd        
            }else{
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[1]] <- NA
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[2]] <- NA
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[3]] <- NA
            }
            
          }else if(types[i.type] == "fd"){
            
            if(length(present.sp)>1){
              lambda.fd <- gradient_fitness_diff(lambda = lambda.obs$lambda,
                                                 steps = steps)
              
              lambda.fd.list <- list()
              for(i.step in 1:length(lambda.fd)){
                lambda.fd.list[[i.step]] <- data.frame(sp = lambda.obs$sp,
                                                       lambda = lambda.fd[[i.step]])
              }
              
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[1]] <- abund.obs
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[2]] <- lambda.fd.list
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[3]] <- alpha.obs
              
            }else{
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[1]] <- NA
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[2]] <- NA
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[3]] <- NA
            }
            
          }else if(types[i.type] == "ia"){
            
            if(length(present.sp)>1){
              alpha.ia <- gradient_asymmetry(A = alpha.obs,steps = steps)
              
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[1]] <- abund.obs
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[2]] <- lambda.obs
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[3]] <- alpha.ia
              
            }else{
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[1]] <- NA
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[2]] <- NA
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[3]] <- NA
            }
            
          }else if(types[i.type] == "id"){
            
            if(length(present.sp)>1){
              alpha.id <- gradient_strength_dist(A = alpha.obs,steps = steps)
              
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[1]] <- abund.obs
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[2]] <- lambda.obs
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[3]] <- alpha.id
              
            }else{
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[1]] <- NA
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[2]] <- NA
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[3]] <- NA
            }
            
          }else if(types[i.type] == "dd"){
            
            if(length(present.sp)>1){
              alpha.id <- gradient_diag_dominance(A = alpha.obs,steps = steps)
              
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[1]] <- abund.obs
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[2]] <- lambda.obs
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[3]] <- alpha.id
              
            }else{
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[1]] <- NA
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[2]] <- NA
              communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]][[3]] <- NA
            }
            
          }
          
          names(communities[[i.model]][[i.year]][[i.plot]][[i.sub]][[i.type]]) <- c("abundances",
                                                                                    "lambda","alpha")
          
        }# for each type
        names(communities[[i.model]][[i.year]][[i.plot]][[i.sub]]) <- types
      }# for i.sub
      names(communities[[i.model]][[i.year]][[i.plot]]) <- subplots
    }# for i.plot
  }# for i.year
  names(communities[[i.model]]) <- years
  
}# for i.model

names(communities) <- model_family

# store results -----------------------------------------------------------

save(communities,file = "results/communities_subplot.Rdata")
