
# obtain steady-state abundances of the perturbed communities
# and their SAD components
source("R/hill_diversity.R")
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


# some parameters ---------------------------------------------------------

# how many timesteps for obtaining average abundances
# the idea is to look at the dynamics before it explodes, as it may do so
# just because of numerical roundings even for stable communities
num.stable.timesteps <- 20

# steps for perturbation intensity
# this should be the same as in "generate_perturbed_communities".
steps <- 10

# obtain metrics ----------------------------------------------------------

# list keeping the metrics for each subplot/plot/year
all.subplots.list <- list()

# loop through years and plots
# i.year <- i.plot <- i.sub <- i.type <- 1
for(i.year in 1:length(years)){
  for(i.plot in 1:length(plots)){

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
          
          if(sum(sum(is.na(r)),sum(is.na(alpha.matrix)),sum(is.na(abund))) == 0){
            
            # alpha.sign <- alpha.matrix * -1
            
            sub.abund <- integrate_GLV(r = r$rfit,
                                       A = alpha.matrix,
                                       x0 = abund$abundance)
                                       # maxtime = 100,
                                       # steptime = .5)
            # test
            # out <- as.data.frame(sub.abund)
            # colnames(out) <- c("time", paste("sp", 1:(ncol(out) -1), sep = "_"))
            # out <- as_tibble(out) %>% gather(species, density, -time)
            # pl <- ggplot(data = out) +
            #     aes(x = time, y = density, colour = species) + ylim(0,1e5) +
            #     geom_line()
            sub.abund.df <- as.data.frame(sub.abund)
            colnames(sub.abund.df) <- c("time",r$sp)
            
            # check where are the last timesteps before the dynamics explode
            stable.timesteps <- apply(sub.abund[,2:ncol(sub.abund)],1,FUN = function(x)sum(x > 100))
            last.stable <- max(which(stable.timesteps == 0))
            # if there are at least X timesteps with controlled dynamics, use these
            # and obtain the mean abundances in them
            # if not?
            if(last.stable > num.stable.timesteps){
              my.abund.df <- sub.abund.df[(last.stable-num.stable.timesteps):last.stable,]
            }else{
              my.abund.df <- NA
            }
            
            if(inherits(my.abund.df,"data.frame")){
              my.abund.long <- gather(my.abund.df,species, density, -time)
              mean.abund <- my.abund.long %>%
                group_by(species) %>%
                summarise(abundance = mean(density))
              
              my.richness <- sum(mean.abund$abundance>0)
              my.abund <- round(sum(mean.abund$abundance))
              my.evenness <- hill.diversity(mean.abund$abundance)
            }else{
              my.richness <- NA_real_
              my.abund <- NA_real_
              my.evenness <- NA_real_
            }
            
            # in the "observed" type, I add "intensity" for standardizing
            # the dataframe, but there is no "intensity" here.
            my.data <- data.frame(year = as.numeric(years[i.year]),
                                  plot = i.plot,
                                  subplot = subplots[i.sub],
                                  type = types[i.type],
                                  intensity = 1:steps,
                                  abundance = my.abund,
                                  richness = my.richness,
                                  evenness = my.evenness)
            
            all.subplots.list[[length(all.subplots.list)+1]] <- my.data
          }# if data ok
        }else if(types[i.type] %in% c("dd","id","ia")){
          
          r <- communities[[years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["r"]]
          abund <- communities[[years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["abundances"]]
          
          for(i.step in 1:steps){
    
            alpha.matrix <- communities[[years[i.year]]][[plots[i.plot]]][[subplots[i.sub]]][[types[i.type]]][["alpha"]][[i.step]]

            if(sum(sum(is.na(r)),sum(is.na(alpha.matrix)),sum(is.na(abund))) == 0){
              
              sub.abund <- integrate_GLV(r = r$rfit,
                                         A = alpha.matrix,
                                         x0 = abund$abundance)

              sub.abund.df <- as.data.frame(sub.abund)
              colnames(sub.abund.df) <- c("time",r$sp)
              
              # check where are the last timesteps before the dynamics explode
              stable.timesteps <- apply(sub.abund[,2:ncol(sub.abund)],1,FUN = function(x)sum(x > 100))
              last.stable <- max(which(stable.timesteps == 0))
              # if there are at least X timesteps with controlled dynamics, use these
              # and obtain the mean abundances in them
              # if not?
              if(last.stable > num.stable.timesteps){
                my.abund.df <- sub.abund.df[(last.stable-num.stable.timesteps):last.stable,]
              }else{
                my.abund.df <- NA
              }
              
              if(inherits(my.abund.df,"data.frame")){
                my.abund.long <- gather(my.abund.df,species, density, -time)
                mean.abund <- my.abund.long %>%
                  group_by(species) %>%
                  summarise(abundance = mean(density))
                
                my.richness <- sum(mean.abund$abundance>0)
                my.abund <- round(sum(mean.abund$abundance))
                my.evenness <- hill.diversity(mean.abund$abundance)
              }else{
                my.richness <- NA_real_
                my.abund <- NA_real_
                my.evenness <- NA_real_
              }
              
              # single row
              my.data <- data.frame(year = as.numeric(years[i.year]),
                                    plot = i.plot,
                                    subplot = subplots[i.sub],
                                    type = types[i.type],
                                    intensity = i.step,
                                    abundance = my.abund,
                                    richness = my.richness,
                                    evenness = my.evenness)
              
              all.subplots.list[[length(all.subplots.list)+1]] <- my.data
              
            }# if !na
          }# for i.step
          
       }# if-else type
        cat(years[i.year],i.plot,subplots[i.sub],types[i.type],"\n",sep="-")
      }# for i.type
      }# if subplot ok
    }# for i.sub
  }# for i.plot
}# for i.year

# clean up results --------------------------------------------------------

pred.plot <- bind_rows(all.subplots.list)

# store results -----------------------------------------------------------

pert.sad <- pivot_longer(pred.plot,cols = abundance:evenness,
                         names_to = "metric",
                         values_to = "value")

write.csv2(pert.sad,file = "results/predicted_SAD_components_subplot_v2.csv",
           row.names = FALSE)
