# tidy empirical network data from Kilnock 2019 AmNat
# NOTE 
# I will not be using these networks in the end
# so writing to disk is not yet implemented


# INPUTS
# - table with competition coefficients, supp mat of Kilnock 19
# https://www.journals.uchicago.edu/doi/abs/10.1086/705293
# - dataframe with coefficients from Caracoles (one network per year)
# "data/alpha_heterogeneous_time.csv"

# OUTPUTS
# - list of network matrices
# "results/plant_competition_networks.RData"
# - dataframe of network info
# "results/plant_competition_networks_info.csv"

# -------------------------------------------------------------------------

library(tidyverse)

# -------------------------------------------------------------------------

k19 <- read.csv("data/other/kilnock19/Data_PlantInteractionNetworks.csv")
caracoles <- read.csv2("data/alpha_heterogeneous_time.csv")

caracoles.lat <- 37.06708333
caracoles.lon <- -6.32116667

k19.1 <- k19[,c("Habitat",
                "Latitude","Longitude",
                "Species","WoodyHerbaceous","PlantPerformanceMetric",
                "TypeOfControl","NetworkMonoCtrl","NetworkTrueCtrl")]

# -------------------------------------------------------------------------

plant.competition.network.info <- NULL
plant.competition.networks <- NULL

# -------------------------------------------------------------------------
# clean caracoles

caracoles.years <- sort(unique(caracoles$year))

# subset weird years (2018)
caracoles.years <- caracoles.years[-4]

caracoles.id <- paste("c",caracoles.years,sep="")

for(ic in 1:length(caracoles.years)){
    my.data <- subset(caracoles, year == caracoles.years[ic] & 
                          !is.na(magnitude))
    
    # my.sp <- sort(unique(my.data$focal))
    
    # interactions not recorded are given values of 0
    
    my.data.clean <- my.data %>%
        filter(!is.na(magnitude)) %>%
        select(focal,neighbour,magnitude) %>%
        distinct() %>%
        subset(neighbour %in% unique(focal)) %>%
        subset(focal %in% unique(neighbour)) %>%
        complete(focal,neighbour,fill = list(magnitude = 0)) %>%
        pivot_wider(names_from = neighbour,values_from = magnitude) 
    
    my.matrix <- as.matrix(my.data.clean[,2:ncol(my.data.clean)])
    rownames(my.matrix) <- colnames(my.matrix)
    
    plant.competition.networks[[length(plant.competition.networks)+1]] <- my.matrix
    
    # -------------------------------------------------------------------------
    # info
    
    my.info <- data.frame(
    net_ID = caracoles.id[ic],
    richness = nrow(my.matrix),
    habitat = "Grassland",
    latitude = caracoles.lat,
    longitude = caracoles.lon,
    plant_morphology = "Herbaceous",
    fitness_metric = "VialbleSeeds",
    reference = "Garcia-Callejas_2021"
    )
    
    plant.competition.network.info <- bind_rows(plant.competition.network.info,
                                                my.info)
}

# -------------------------------------------------------------------------
# clean kinlock 19

kinlock.id <- paste("k",sprintf("%02d", 1:nrow(k19.1)),sep="")

for(ik in 1:length(kinlock.id)){
    
    # -------------------------------------------------------------------------
    # matrix
    
    my.sp <- unlist(str_split(k19.1$Species[ik],", "))
    if(k19.1$TypeOfControl[ik] == "True ctrl"){
        my.coefs <- as.numeric(unlist(str_split(k19.1$NetworkTrueCtrl[ik],",")))
    }else{
        my.coefs <- as.numeric(unlist(str_split(k19.1$NetworkMonoCtrl[ik],",")))
    }
    
    my.matrix <- matrix(my.coefs,
                        nrow = length(my.sp),
                        ncol = length(my.sp),
                        byrow = TRUE,
                        dimnames = list(my.sp,my.sp))
    
    plant.competition.networks[[length(plant.competition.networks)+1]] <- my.matrix
    
    # -------------------------------------------------------------------------
    # info
    
    my.info <- data.frame(
        net_ID = kinlock.id[ik],
        richness = nrow(my.matrix),
        habitat = k19.1$Habitat[ik],
        latitude = k19.1$Latitude[ik],
        longitude = k19.1$Longitude[ik],
        plant_morphology = k19.1$WoodyHerbaceous[ik],
        fitness_metric = k19.1$PlantPerformanceMetric[ik],
        reference = "Kilnock_2019"
    )
    
    plant.competition.network.info <- bind_rows(plant.competition.network.info,
                                                my.info)
    }

# TODO write to disk


