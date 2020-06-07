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
source("R/removeNiches.R")
source("R/model_abund_negbin.R")

store.results <- F

##########################
# from which method are parameters obtained
method <- "negbin"
# method <- "bobyqa"

##########################
# models to generate
models <- c("complete", 
            "no_niche" 
            # "no_fitness"
            )
##########################

# abundances and vital rates
abundances <- read.csv2(file = "../Caracoles/data/abundances.csv",header = T,stringsAsFactors = F)
sp.list <- read.csv2(file = "../Caracoles/data/plant_species_traits.csv",header = T,stringsAsFactors = F)

valid.sp <- sp.list$species.code[which(!is.na(sp.list$germination.rate))]

abund.plot.year <- abundances %>% 
  filter(species %in% valid.sp) %>% 
  group_by(year,plot,species) %>% 
  summarise(abund = sum(individuals))

# germination and survival rates
germ.rate.orig <- sp.list$germination.rate[match(sort(valid.sp),sp.list$species.code)]
surv.rate.orig <- sp.list$seed.survival[match(sort(valid.sp),sp.list$species.code)]

# model parameters
if(method == "negbin"){
  alpha.orig <- read.csv2(file = "../Spatial_coexistence/results/alpha_negbin_multilevel_heterogeneous_both_clean.csv",header = T,stringsAsFactors = F)
  lambda.orig <- read.csv2(file = "../Spatial_coexistence/results/lambda_negbin_multilevel_heterogeneous_both_clean.csv",header = T,stringsAsFactors = F)
  
  # lambda should be in its original scale
  lambda.orig$lambda <- log(lambda.orig$lambda)
  
  abundance.model <- model_abund_negbin
}else{
  alpha.orig <- read.csv2(file = paste("../Spatial_coexistence/results/cxr_alpha_plot_per_year_BH3_",optim.method,".csv",sep=""),header = T,stringsAsFactors = F)
  lambda.orig <- read.csv2(file = paste("../Spatial_coexistence/results/cxr_lambda_plot_per_year_BH3_",optim.method,".csv",sep=""),header = T,stringsAsFactors = F)
  abundance.model <- cxr::model_abundBH3
}
alpha.orig <- subset(alpha.orig,focal %in% valid.sp & competitor %in% valid.sp)
lambda.orig <- subset(lambda.orig,sp %in% valid.sp)

# initial values and params for abundance projection
timesteps <- 2

init.abund <- abund.plot.year[abund.plot.year$year == 2019,c("plot","species","abund")]
names(init.abund) <- c("site","species","abundance")

##########
# stick with plot 1 for now
# abundance
init.abund <- subset(init.abund, site == 2)

# lambda
lambda.orig <- subset(lambda.orig, year == 2019 & plot == 2)

# lambda for non-present species? 
# TODO for now, zero
absent.sp <- lambda.orig$sp[which(is.na(lambda.orig$lambda))]
lambda.orig$lambda[which(is.na(lambda.orig$lambda))] <- 0#sp.list$lambda[match(absent.sp,sp.list$species.code)]
# lambda.orig$lambda <- as.numeric(lambda.orig$lambda)

# alpha
alpha.orig <- subset(alpha.orig,year == 2019 & plot == 1)
alpha.orig.m <- tidyr::spread(alpha.orig,competitor,magnitude)
alpha.orig.mat <- as.matrix(alpha.orig.m[,4:ncol(alpha.orig.m)])
rownames(alpha.orig.mat) <- alpha.orig.m$focal
alpha.orig.mat[is.na(alpha.orig.mat)] <- 0
##########

# results dataframe
abundance.projections <- NULL

# for each model, generate SAD
# using predictAbundances function

for(i.model in 1:length(models)){
  if(models[i.model] == "complete"){
    lambda.model <- lambda.orig$lambda
    germ.rate.model <- germ.rate.orig
    surv.rate.model <- surv.rate.orig
    alpha.model <- alpha.orig.mat
  }else if(models[i.model] == "no_fitness"){
    lambda.model <- rep(mean(lambda.orig$lambda),length(lambda.orig$lambda))
    germ.rate.model <- rep(mean(germ.rate.orig),length(germ.rate.orig))
    surv.rate.model <- rep(mean(surv.rate.orig),length(surv.rate.orig))
    alpha.model <- alpha.orig.mat
  }else if(models[i.model] == "no_niche"){
    lambda.model <- lambda.orig$lambda
    germ.rate.model <- germ.rate.orig
    surv.rate.model <- surv.rate.orig
    alpha.model <- removeNiches(alpha.orig.mat)
    alpha.model[which(is.na(alpha.model))] <- 0
  }
  
  sp.par <- data.frame(species = 1:length(valid.sp),
                       lambda = lambda.model,
                       germ.rate = germ.rate.model, 
                       survival.rate = surv.rate.model,stringsAsFactors = F)
  
  par <- list(sp.par = sp.par, 
              initial.values = init.abund, 
              covariates = 0,
              other.par = list(alpha.matrix = alpha.model))
  
  predicted.abundances <- PredictAbundances(par = par,
                                            timesteps = timesteps,
                                            abundance.model = abundance.model,
                                            return.seeds = FALSE)
  # include a "model" column
  predicted.abundances$model <- models[i.model]
  
  # append results
  abundance.projections <- rbind(predicted.abundances,abundance.projections)
  
}# for i.model

if(store.results){
  write.csv2(abundance.projections,paste("./results/abundance_projections_",optim.method,".csv",sep=""))
}

####################
# rank-abundance curves
rank.abundances <- abundance.projections %>% group_by(timestep,site,model) %>% mutate(relative.abund = abundance/sum(abundance))
rank.abundances <- rank.abundances %>% group_by(site,timestep,model) %>% mutate(species.rank = rank(-relative.abund,ties.method = "first"))

my.palette <- c("darkgreen","#009E73","#E69F00","#D55E00")
rad.plot <- ggplot(rank.abundances,aes(x = species.rank,y = relative.abund, color = model, group = model)) +
  geom_point() +
  geom_line() +
  facet_grid(model~timestep,scales = "free")+
  scale_color_manual(values = my.palette)+
  # DGC::theme_Publication()+
  # theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  # guides(color=FALSE)+
  NULL
rad.plot


# density plot
density.plot <- ggplot(abundance.projections)+
  geom_density(aes(x = log(abundance),group = model, fill = model), alpha = .5)+
  facet_grid(timestep~.)+
  # xlim(0,1000)+
  NULL
density.plot

