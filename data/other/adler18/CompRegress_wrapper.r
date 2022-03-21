
# this script calls other scripts to compete all stages of analysis

rm(list=ls())

# make the directory containing this script the working directory
setwd("C:/Users/adler/Box Sync/activeProjects/CompetitionReview/datapluscode")

# required packages
req_libs <- c("lme4","tidyr", "dplyr", "ggplot2","ggthemes","viridis","gridExtra",
              "grid","gtable","cowplot","texreg")

# use this function to check if each package is on the local machine
# if a package is installed, it will be loaded
# if any are not, the missing package(s) will be installed and loaded
package_check <- lapply(req_libs, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# import data
rawD <- read.csv("CompetitionRegressionData071517.csv")

# clean up data
source("CompRegress_cleanData.r")  # this returns a "cleanD" data frame
rm(rawD)  # use cleanD from here on

# scale the competition coefficients
source("CompRegress_scaleAlphas.r")  # this returns "effect" and "response" data frames

# analyze pairwise comparisons (effects and responses)
source("CompRegress_analyzePairs.r")

# analyze rho (niche differences)
source("CompRegress_analyzeRho.r")

# make figures
source("CompRegress_figures.r")

# simulation of rho, competitive effects and competitive responses
source("CompRegress_simulation.r")

# The Discussion includes some quantitative results from
# other recent studies that were not included in our original search.
# To see where those numbers come from, open "other_studies.R",
# download and save the data from the listed url's, and then
# execute the code.
