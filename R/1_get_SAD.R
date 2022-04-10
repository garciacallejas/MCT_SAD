
library(tidyverse)
source("R/GLV_functions.R")
source("R/hill_diversity.R")

# -------------------------------------------------------------------------

# 1 - solve Lotka-Volterra systems associated to the matrices generated in silico
# 2 - for those with > N surviving species, recover their SAD
# 3 - store, for each system selected, their full SAD, 
# alongside the information of the interaction matrix (degree of each "treatment")

# is it ok to assume r = 1? so far, go with it
# is it ok to assume equal abundances at the beginning of the simulation? likewise

load("results/community_matrices_silico.RData")
gradient.df <- read.csv2("results/community_matrices_silico_parameters.csv")

# -------------------------------------------------------------------------

SAD.metrics.df <- gradient.df
SAD.metrics.df$steady.state.richness <- NA
SAD.metrics.df$steady.state.abundance <- NA
SAD.metrics.df$steady.state.evenness <- NA

richness <- unique(gradient.df$richness)
sp.names <- paste("sp",1:richness,sep="")

SAD.list <- list()

for(i in 1:nrow(SAD.metrics.df)){
    A <- matrix.list[[i]]
    r <- rep(1,gradient.df$richness[i])
    x0 <- rep(.1,gradient.df$richness[i])
    
    lv.dynamics <- integrate_GLV(r = r,
                                 A = A,
                                 x0 = x0,
                                 positive.competition = TRUE)
    
    # pl <- ggplot(data = lv.dynamics) +
    #     aes(x = time, y = density, colour = species) + #ylim(0,1e5) +
    #     geom_line()
    
    steady.state.abundances <- subset(lv.dynamics, time == max(time))$density
    steady.state.abundances[steady.state.abundances < 1e-5] <- 0
    
    my.abundances <- expand_grid(richness = SAD.metrics.df$richness[i],
                                 connectance = SAD.metrics.df$connectance[i],
                                 diagonal.dominance = SAD.metrics.df$diagonal.dominance[i],
                                 mean.strength = SAD.metrics.df$mean.strength[i],
                                 # tau = SAD.metrics.df$tau[i],
                                 replicate = SAD.metrics.df$replicate[i],
                                 matrix.code = SAD.metrics.df$matrix_code[i],
                                 species = sp.names,
                                 abundance = NA)
    my.abundances$abundance <- steady.state.abundances
    
    SAD.metrics.df$steady.state.richness[i] <- sum(steady.state.abundances>0)
    SAD.metrics.df$steady.state.abundance[i] <- round(sum(steady.state.abundances))
    # this equation returns evenness = 1 for single-species communities, so correct it
    # and assign NA instead
    SAD.metrics.df$steady.state.evenness[i] <- ifelse(sum(steady.state.abundances>0)>1,
                                                      hill.diversity(steady.state.abundances)/sum(steady.state.abundances>0),
                                                      NA_real_)
    
    SAD.list[[i]] <- my.abundances
    
}

SAD.abundances <- bind_rows(SAD.list)

# -------------------------------------------------------------------------

write.csv2(SAD.metrics.df,"results/SAD_metrics_silico.csv",row.names = FALSE)
write.csv2(SAD.abundances, "results/SAD_abundances_silico.csv", row.names = FALSE)

