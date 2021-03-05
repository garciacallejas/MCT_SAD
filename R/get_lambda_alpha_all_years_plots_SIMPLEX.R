
library(tidyverse)
library(cxr)

source("./R/AP_pm_alpha_pairwise_lambdacov_none_alphacov_none.R")
source("./R/AP_project_alpha_pairwise_lambdacov_none_alphacov_none.R")
source("./R/AP_pm_alpha_pairwise_lambdacov_global_alphacov_global.R")
source("./R/AP_project_alpha_pairwise_lambdacov_global_alphacov_global.R")
source("./R/AP_pm_alpha_pairwise_lambdacov_global_alphacov_pairwise.R")
source("./R/AP_project_alpha_pairwise_lambdacov_global_alphacov_pairwise.R")

# read abundance data per plot and year -----------------------------------

abund <- read.csv2("../Caracoles/data/abundances.csv",
                   header = TRUE,stringsAsFactors = FALSE)
comp <- read.csv2("../Caracoles/data/competition_wide.csv",
                  header = TRUE,stringsAsFactors = FALSE)
sp.rates <- read.csv2("../Caracoles/data/plant_species_traits.csv",
                      header = TRUE,stringsAsFactors = FALSE)
sp.valid <- sp.rates$species.code[which(!is.na(sp.rates$germination.rate))]

# covariate - precipitation
prec <- read.csv2("../Caracoles/data/precipitation_per_year.csv")
# covariate: salinity
salinity <- read.csv2("../Caracoles/data/salinity.csv")

base.abund <- abund %>% 
  filter(species %in% sp.valid) 
# group_by(species) %>% 
# summarise(num = sum(individuals))
# year.off <- 2019

comp <- comp[,which(names(comp) %in% c("year","plot","subplot",
                                       "focal","fruit","seed",sp.valid))]
comp <- subset(comp, seed > 0)
# comp <- subset(comp, year != year.off)

all.sp <- names(comp)[7:length(names(comp))]
focal.sp <- unique(comp$focal[which(comp$focal %in% all.sp)])

# create multifit list
neigh.list <- list()
# precipitation covariate
sal.precip.list <- list()

for(i.sp in 1:length(focal.sp)){
  my.data <- subset(comp,focal == focal.sp[i.sp])
  my.data2 <- left_join(my.data,salinity)
  # if(nrow(my.data)>0){
  my.data <- my.data[,6:length(my.data)]
  names(my.data)[1] <- "fitness"
  neigh.list[[i.sp]] <- my.data
  names(neigh.list)[i.sp] <- focal.sp[i.sp]
  
  sal.precip.list[[i.sp]] <- data.frame(co3 = my.data2[,"co3"])
                                        # prec = my.data2[,"precip"])
  # }# if observations
}

# names(neigh.list) <- all.sp

# initial values ----------------------------------------------------------

# load model
# this is an Annual Plant model, 
# using the Ricker model from Mayfield and Stouffer 2017
# We constrain alphas to be positive,


# TEST standard Ricker model
# the AP model includes g,s, which are species-specific, 
# and therefore would contribute to fitness diff. If I want 
# to have complete control, I need a model only with lambda/alphas
model_family <- "AP"
# model_family <- "RK"

optimization_method <- "bobyqa"
alpha_form <- "pairwise"
lambda_cov_form <- "none"
alpha_cov_form <- "none"
# lambda_cov_form <- "global"
# alpha_cov_form <- "global"

fixed_terms <- NULL
bootstrap_samples <- 0

# initial_values <- list(lambda = 10, alpha_intra = 0.01, alpha_inter = 0.01)
# # to remove niche diff, we need positive coefficients
# lower_bounds <- list(lambda = 0, alpha_intra = 0, alpha_inter = 0)
# upper_bounds <- list(lambda = 1e3, alpha_intra = 1, alpha_inter = 1)

initial_values <- list(lambda = 10, 
                       alpha_intra = 0.01, 
                       alpha_inter = 0.01)
                       # lambda_cov = 0.01,
                       # alpha_cov = 0.01)
# to remove niche diff, we need positive coefficients
lower_bounds <- list(lambda = 0, 
                     alpha_intra = 0, 
                     alpha_inter = 0) 
                     # lambda_cov = 0,
                     # alpha_cov = 0)

upper_bounds <- list(lambda = 1e3, 
                     alpha_intra = 1, 
                     alpha_inter = 1) 
                     # lambda_cov = 1,
                     # alpha_cov = 1)

# fit model ---------------------------------------------------------------

# include environmental covariates? not for now
caracoles.fit <- cxr_pm_multifit(data = neigh.list,
                                 model_family = model_family,
                                 focal_column = all.sp,
                                 # covariates = sal.precip.list,
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

# store lambdas and alphas ------------------------------------------------

sp.lambdas <- data.frame(sp = focal.sp, lambda = caracoles.fit$lambda)
alpha.matrix <- caracoles.fit$alpha_matrix
# lambda_cov <- data.frame(sp = focal.sp, lambda = caracoles.fit$lambda_cov)
# alpha_cov <- caracoles.fit$alpha_cov[["co3"]]

pred <- list(lambda = sp.lambdas, alpha = alpha.matrix) 
             # lambda_cov = lambda_cov, alpha_cov = alpha_cov)
save(pred,file = "./results/SIMPLEX_fit_AP.Rdata")

# 
# write.csv2(sp.lambdas, file = "./results/lambda_salinity_AP.csv",row.names = FALSE)
# write.csv2(alpha.matrix,file = "./results/alpha_salinity_AP.csv")
# write.csv2(lambda_cov, file = "./results/lambda_cov_salinity_AP.csv",row.names = FALSE)
# write.csv2(alpha_cov,file = "./results/alpha_cov_salinity_AP.csv",row.names = TRUE)
