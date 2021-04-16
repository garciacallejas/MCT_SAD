
# obtain steady-state abundances of the perturbed communities
# and their SAD components

source("R/GLV_functions.R")

# read data ---------------------------------------------------------------

# load the perturbed communities
load("results/communities_subplot_perturbed.Rdata")

years <- names(communities)
plots <- 1:length(communities[[1]])
subplots <- names(communities[[1]][[1]])
types <- names(communities[[1]][[1]][[1]])
steps <- length(communities[[1]][[1]][[1]][["ia"]][["alpha"]])

metrics <- c("richness","abundance","evenness")

# obtain metrics ----------------------------------------------------------

# list keeping the metrics for each plot/year
all.plots.list <- list()

# loop through years and plots
# i.year <- i.plot <- i.sub <- i.type <- 1
for(i.year in 1:length(initial.years)){
  for(i.plot in 1:length(plots)){
    
    # list keeping all raw abundance projections
    plot.list <- list()
    count <- 1
    
    # loop through subplots
    for(i.sub in 1:length(subplots)){
      
      # if this subplot has a valid corrected set of parameters
      if(inherits(communities[[years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]],"list")){
        
      # obtain SAD components from every perturbation type, including the original SAD
      for(i.type in 1:length(types)){
        
        if(types[i.type] == "obs"){
          
          r <- communities[[years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["r"]]
          alpha.matrix <- communities[[years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["alpha"]]
          abund <- communities[[years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["abundances"]]
          
          if(sum(sum(is.na(r)),sum(is.na(alpha)),sum(is.na(abund))) == 0){
            
            # alpha.sign <- alpha.matrix * -1
            
            sub.abund <- integrate_GLV(r = r$rfit,
                                       A = alpha.matrix,
                                       x0 = abund$abundance)
            # these should be stable, but they are not
            
          }# if data ok
        }# if type
      }# for i.type
      }# if subplot ok
    }# for i.sub
  }# for i.plot
}# for i.year