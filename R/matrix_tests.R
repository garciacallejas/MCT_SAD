# simple test
library(cxr)
source("R/removeNiches.R")

tm <- matrix(c(1,2,3,4),nrow = 2)

rn <- removeNiches(tm)
# should be 1
NicheOverlap(rn)
# NicheOverlap(tm)

rf <- removeFitness(tm)
AvgFitnessRatio(lambda = c(1.2,1.2),pair.matrix = rf)
AvgFitnessRatio(lambda = c(1.2,1.2),pair.matrix = tm)

pair.matrix <- rf
pair.matrix <- tm
pair.matrix <- rn
sqrt((pair.matrix[1,2]/pair.matrix[2,2])*(pair.matrix[1,1]/pair.matrix[2,1]))
