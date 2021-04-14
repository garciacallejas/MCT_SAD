
library(tidyverse)
library(cxr)

# read abundance data per plot and year -----------------------------------

abund <- read.csv2("results/abund_filtered.csv",
                   header = TRUE,stringsAsFactors = FALSE)
comp <- read.csv2("results/neigh_filtered.csv",
                  header = TRUE,stringsAsFactors = FALSE)

comp <- subset(comp, seed > 0)

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

# load models

# Test all standard models in cxr
<<<<<<< HEAD
model_family <- c("LV")
=======
# model_family <- c("BH","LW","RK")

# try a High Density dependence (HD) model,
# based on the ricker formulation
model_family <- "HD"
source("R/HD_pm_alpha_pairwise_lambdacov_none_alphacov_none.R")
>>>>>>> 4355279d30e67477a2a968de83466284d6a57b70

optimization_method <- "bobyqa"
alpha_form <- "pairwise"
lambda_cov_form <- "none"
alpha_cov_form <- "none"

fixed_terms <- NULL
bootstrap_samples <- 0

# initial_values <- list(lambda = 10, alpha_intra = 0.01, alpha_inter = 0.01)
# 
# lower_bounds <- list(lambda = 0, alpha_intra = 0, alpha_inter = 0)
# upper_bounds <- list(lambda = 1e3, alpha_intra = 1, alpha_inter = 1)

initial_values <- list(lambda = 1, alpha_intra = 0, alpha_inter = 0)

lower_bounds <- list(lambda = 0, alpha_intra = 0, alpha_inter = 0)
upper_bounds <- list(lambda = 10, alpha_intra = 1, alpha_inter = 0.5)


# fit models --------------------------------------------------------------

for(i.model in 1:length(model_family)){
  
  # include environmental covariates? not for now
  caracoles.fit <- cxr_pm_multifit(data = neigh.list,
                                   model_family = model_family[i.model],
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
  
  # store lambdas and alphas ------------------------------------------------
  
  sp.lambdas <- data.frame(sp = focal.sp, lambda = caracoles.fit$lambda)
  alpha.matrix <- caracoles.fit$alpha_matrix
  
  write.csv2(sp.lambdas, file = paste("./results/lambda_",model_family[i.model],".csv",sep=""),row.names = FALSE)
  write.csv2(alpha.matrix,file = paste("./results/alpha_",model_family[i.model],".csv",sep=""))
  
}
