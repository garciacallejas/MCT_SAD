
library(tidyverse)
source("R/horizontal_community_matrix.R")

# -------------------------------------------------------------------------

# TODO higher tau variability, richness = 50

# generate interaction matrices with constraints on
# 1 - diagonal dominance
# 2 - connectance
# 3 - heterogeneity of the interaction strength distribution

# richness
S <- 25

# replicates of each combination
replicates <- 10

# gradient in connectance
connectance.gradient <- seq(from = 0.05,to = 0.75, length.out = 10)
# diagonal dominance
diag.dom.gradient <- seq(from = 0, to = 1, length.out = 10)
# kurtosis: 
# low tau values mean high kurtosis and viceversa
# as per preliminary tests, 1.5 is approximate to a normal dist 
# (using rSHASHO to generate random samples)
kurtosis.tau.gradient <- seq(from = 3, to = 1, length.out = 10)
# species richness
mean.strength.gradient <- seq(from = 0, to = 2, length.out = 10)
# -------------------------------------------------------------------------

# this dataframe will hold information about each matrix
gradient.df <- expand_grid(richness = S,
                           connectance = connectance.gradient,
                           diagonal.dominance = diag.dom.gradient,
                           # tau = kurtosis.tau.gradient,
                           mean.strength = mean.strength.gradient,
                           replicate = 1:replicates) %>%
    mutate(matrix_code = paste("c",round(connectance,2),"_d",round(diagonal.dominance,2),
                               "_s",round(mean.strength,2),"_r",replicate,sep=""))

# this list will hold the actual matrices
matrix.list <- list()

# -------------------------------------------------------------------------

for(i in 1:nrow(gradient.df)){
    matrix.list[[i]] <- horizontal_community_matrix(S = S,
                                                    c = gradient.df$connectance[i],
                                                    # tau = gradient.df$tau[i],
                                                    min.diag.dom = gradient.df$diagonal.dominance[i],
                                                    int.mean = gradient.df$mean.strength[i],
                                                    restricted.positive = TRUE)
}

names(matrix.list) <- gradient.df$matrix_code

# -------------------------------------------------------------------------

save(matrix.list,file = "results/community_matrices_silico.RData")
write.csv2(gradient.df,"results/community_matrices_silico_parameters.csv",row.names = FALSE)
