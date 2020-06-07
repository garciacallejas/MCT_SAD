
library(tidyverse)
library(patchwork)


# read different datasets -------------------------------------------------

training.data.type <- c("yearly","spatial")
covariates <- c("no_precip","precip")

caracoles.projections <- NULL
for(i.tr in 1:length(training.data.type)){
  for(i.cov in 1:length(covariates)){
    my.data <- read.csv2(file = paste("./results/projected_abundances_",training.data.type[i.tr],"_",covariates[i.cov],".csv",sep=""),
                         header = TRUE,stringsAsFactors = FALSE)
    my.data$training.type <- training.data.type[i.tr]
    my.data$covariates <- covariates[i.cov]
    caracoles.projections <- rbind(caracoles.projections,my.data)
  }
}

# RMSE --------------------------------------------------------------------

# rmse is extremely high in all cases, but it is orders of magnitude higher
# without covariates, in particular because it produces some extreme outliers
# with predicted abundances of ~1e7. 
# Even discarding these outliers, rmse is one order of magnitude higher
# without covariates

rmse <- caracoles.projections %>%
  filter(!is.na(predicted)) %>%
  group_by(training.type,covariates) %>% 
  summarise(rmse = sqrt(sum((predicted - observed)^2)/n()))

# plot --------------------------------------------------------------------

abund.plot <- ggplot(caracoles.projections,aes(x = log(observed),y = log(predicted),group = species)) + 
  geom_point(aes(color = species)) + 
  geom_abline(slope = 1,linetype = "dashed",color = "grey") +
  facet_grid(training.type~covariates) + 
  theme_bw() +
  # xlim(0,50)+ylim(0,50)+
  ggtitle("log predicted-observed")+
  NULL
# abund.plot

oabund.plot <- ggplot(caracoles.projections,aes(x = observed, y = predicted,group = species)) + 
  geom_point(aes(color = species)) + 
  geom_abline(slope = 1,linetype = "dashed",color = "grey") +
  facet_grid(training.type~covariates) + 
  theme_bw() +
  guides(color = FALSE) +
  xlim(0,50)+ylim(0,50)+
  ggtitle("predicted-observed")+
  NULL
# oabund.plot

tt <- oabund.plot + abund.plot

ggsave(filename = paste("results/images/SIMPLEX_mechanistic_projected_abundances.pdf",sep=""),plot = tt,device = cairo_pdf,width = 12,height = 6,dpi = 600)
write.csv2(caracoles.projections,file = "results/projected_abundances_mechanistic_model.csv",row.names = FALSE)

