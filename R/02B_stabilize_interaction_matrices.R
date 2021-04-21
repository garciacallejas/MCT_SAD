
library(tidyverse)
source("other_code/reconciling_coexistence-master/reconciling_coexistence-master/code/functions_algae.R")
source("R/GLV_functions.R")

# read data ---------------------------------------------------------------

model_family <- "LV"

lambda <- read.csv2(paste("./results/lambda_",model_family,".csv",sep=""),stringsAsFactors = FALSE)
alpha.df <- read.csv2(paste("results/alpha_",model_family,".csv",sep=""),stringsAsFactors = FALSE,
                      row.names = 1)
alpha.matrix <- as.matrix(alpha.df)

abund <- read.csv2("results/abund_filtered.csv",
                   header = TRUE,stringsAsFactors = FALSE)
base.abund <- abund %>% 
  filter(species %in% rownames(alpha.matrix)) 

years <- sort(unique(abund$year))
plots <- sort(unique(abund$plot))
subplots <- sort(unique(abund$subplot))

# -------------------------------------------------------------------------
# list with observed abundances, corrected interaction matrices, and growth rates
# per subplot and year

stable_communities <- list()

# i.year <- i.plot <- i.sub <- 1
for(i.year in 1:length(years)){
  
  stable_communities[[i.year]] <- list()
  
  for(i.plot in 1:length(plots)){
    
    stable_communities[[i.year]][[i.plot]] <- list()
    
    for(i.sub in 1:length(subplots)){
      
      stable_communities[[i.year]][[i.plot]][[i.sub]] <- list()
      
      # subset present species
      present.sp <- sort(unique(base.abund$species[base.abund$year == years[i.year] &
                                                     base.abund$plot == plots[i.plot] &
                                                     base.abund$subplot == subplots[i.sub] &
                                                     base.abund$individuals > 0]))
      
      if(length(present.sp)>1){
        
        # sum observed abundances
        # surely not needed at the subplot level...
        abund.obs <- base.abund %>% 
          filter(year == years[i.year] &
                   plot == plots[i.plot] &
                   subplot == subplots[i.sub] &
                   species %in% present.sp) %>%
          group_by(species) %>%
          summarise(abundance = sum(individuals))
        
        lambda.obs <- lambda[lambda$sp %in% present.sp,]
        lambda.obs <- arrange(lambda.obs,sp)
        
        # the quadratic programming function accepts negative coefficients
        # for competition
        alpha.obs <- alpha.matrix[present.sp,present.sp] * -1
        alpha.obs[which(is.na(alpha.obs))] <- 0
        
        corrected_A <- fit_qp_LV(A=alpha.obs,
                                 r=lambda.obs$lambda,
                                 x_obs=abund.obs$abundance,
                                 tol=1000)
      
        if(!inherits(corrected_A,"warning")){
          nspp <- ncol(alpha.obs)
          Afit <- t(matrix(corrected_A$X[1:nspp^2], nspp,nspp,
                           dimnames = list(rownames(alpha.obs),colnames(alpha.obs))))
          rfit <- data.frame(sp = lambda.obs$sp, rfit = corrected_A$X[(nspp^2+1):(nspp^2+nspp)])
          # round(Afit,2)
          # x_obs <- abund.obs$abundance
          # # exact solution?
          # round(x_obs*(rfit$rfit+Afit%*%x_obs),12)
          # # steady state SAD
          # test <- integrate_GLV(r = rfit$rfit,
          #                       A = Afit,
          #                       x0 = x_obs)
          # out <- as.data.frame(test)
          # colnames(out) <- c("time", paste("sp", 1:(ncol(out) -1), sep = "_"))
          # out <- as_tibble(out) %>% gather(species, density, -time)
          # pl <- ggplot(data = out) +
          #     aes(x = time, y = density, colour = species) +
          #     geom_line()
          
          stable_communities[[i.year]][[i.plot]][[i.sub]][[1]] <- abund.obs
          stable_communities[[i.year]][[i.plot]][[i.sub]][[2]] <- rfit
          stable_communities[[i.year]][[i.plot]][[i.sub]][[3]] <- Afit
          # cat("\n",years[i.year],i.plot,subplots[i.sub],"OK",sep="-")
        }else{
          stable_communities[[i.year]][[i.plot]][[i.sub]][[1]] <- NA
          stable_communities[[i.year]][[i.plot]][[i.sub]][[2]] <- NA
          stable_communities[[i.year]][[i.plot]][[i.sub]][[3]] <- NA
          # cat("\n",years[i.year],i.plot,subplots[i.sub],"--- NOT OK ---",sep="-")
        }
      }else{
        stable_communities[[i.year]][[i.plot]][[i.sub]][[1]] <- NA
        stable_communities[[i.year]][[i.plot]][[i.sub]][[2]] <- NA
        stable_communities[[i.year]][[i.plot]][[i.sub]][[3]] <- NA
        # cat("\n-1sp",years[i.year],i.plot,subplots[i.sub],"--- NOT OK ---",sep="-")
      }
      
      names(stable_communities[[i.year]][[i.plot]][[i.sub]]) <- c("abundances",
                                                                  "corrected_r",
                                                                  "corrected_alpha")
      
    }# for i.sub
    names(stable_communities[[i.year]][[i.plot]]) <- subplots
  }# for i.plot
}# for i.year
names(stable_communities) <- years

# write to disk -----------------------------------------------------------

save(stable_communities,file = "results/communities_subplot_stabilized.Rdata")
