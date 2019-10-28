# generate SAD from a community given model parameters
# generate SADs for 
# 1) the complete model
# 2) model without niche differences
# 3) model without fitness demographic ratio
# it is unclear -yet- if competitive response ratio can be knocked off
# use caracoles data

# devtools::install_github("ibartomeus/cxr")
library(cxr)
library(tidyverse)

# models to generate
models <- c("complete", 
            "no_niche", 
            "no_fitness_dem")

# read model parameters from
# previous simulations
# lambda in sorted numeric vector
lambda.orig <- 1
# alpha in sorted numeric matrix
alpha.orig <- 1

# initial values and params for abundance projection
abundance.model <- model_abundBH3
timesteps <- 2

init.abund <- expand.grid(1:num.obs,1:num.sp)
names(init.abund) <- c("site","species")
init.abund$abundance <- rnorm(nrow(init.abund),10,2)

# results dataframe
abundance.projections <- NULL

# for each model, generate SAD
# using predictAbundances function

for(i.model in 1:length(models)){
  if(models(i.model) == "complete"){
    lambda.model <- lambda.orig
    germ.rate.model <- germ.rate.orig
    surv.rate.model <- surv.rate.orig
    alpha.model <- alpha.orig
  }else if(models(i.model) == "no_fitness"){
    lambda.model <- mean(lambda.orig)
    germ.rate.model <- mean(germ.rate.orig)
    surv.rate.model <- mean(surv.rate.orig)
    alpha.model <- alpha.orig
  }else if(models(i.model) == "no_fitness"){
    lambda.model <- mean(lambda.orig)
    germ.rate.model <- mean(germ.rate.orig)
    surv.rate.model <- mean(surv.rate.orig)
    alpha.model <- removeNiches(alpha.orig)
  }
  
  sp.par <- data.frame(species = 1:num.sp,
                       lambda = lambda.model,
                       germ.rate = germ.rate.model, 
                       survival.rate = surv.rate.model)
  
  par <- list(sp.par = sp.par, 
              initial.values = init.abund, 
              covariates = 0,
              other.par = list(alpha.matrix = alpha.model))
  
  predicted.abundances <- PredictAbundances(par = par,
                                            timesteps = timesteps,
                                            abundance.model = abundance.model,
                                            return.seeds = TRUE)
  # include a "model" column
  abundance.projections <- rbind(predicted.abundances,abundance.projections)
  
}# for i.model





