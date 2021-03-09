# plot SAD components
library(tidyverse)
library(colorblindr)
library(patchwork)

# read data ---------------------------------------------------------------

obs <- read.csv2(file = "results/observed_SAD_components_subplot.csv",
                 stringsAsFactors = FALSE)
pred <- read.csv2(file = "results/predicted_SAD_components_subplot_aggregated.csv",
                  stringsAsFactors = FALSE)

# pred <- subset(pred, timestep == 2)

# plots -------------------------------------------------------------------

# -------------------------------------------------------------------------
# 1 - temporal trends for all models and types of matrix perturbations

# select an intensity value for this plot. Note that for the unperturbed matrix
# all intensity values are equivalent, they are put there for comparison.
pred.2 <- pred %>% filter(intensity == 5) %>%
  group_by(model,type,timestep,metric) %>%
  summarise(min.value = min(value),
            mean.value = mean(value),
            max.value = max(value),
            sd.value = sd(value))

# add log(abundance) category
pred.log.abund <- pred.2 %>% filter(metric == "abundance")
pred.log.abund$min.value <- log(pred.log.abund$min.value)
pred.log.abund$mean.value <- log(pred.log.abund$mean.value)
pred.log.abund$sd.value <- log(pred.log.abund$sd.value)
pred.log.abund$max.value <- log(pred.log.abund$max.value)
pred.log.abund$metric <- "log(abundance)"

pred.3 <- rbind(pred.2,pred.log.abund) %>% filter(metric != "abundance")

# plot
pd <- 0.5
ew <- 0.5
dk <- 0.2

f1 <- ggplot(pred.3,aes(x = timestep, y = mean.value, group = model)) + 
  geom_errorbar(aes(ymin = mean.value-sd.value, ymax = mean.value+sd.value,
                    color = model),
                position = position_dodge(pd),
                width = ew)+
  geom_line(aes(color = model),position = position_dodge(pd)) +
  geom_point(aes(fill = model), 
             shape = 21,
             position = position_dodge(pd),
             size = 2) +
  facet_grid(metric~type, scales = "free_y") +
  theme_bw() +
  theme(strip.background = element_blank())+
  scale_fill_OkabeIto(darken = dk) +
  scale_color_OkabeIto(darken = dk) +
  ylab("value")+
  NULL
# f1

ggsave(filename = "results/images/timestep_trends.pdf",
       plot = f1,
       device = cairo_pdf,
       width = 10,height = 10,dpi = 300)

# -------------------------------------------------------------------------
# 2 - variations with respect to the standard matrix
# save a plot for each model and timestep

model_family <- unique(pred$model)
timesteps <- max(pred$timestep)

pred.wide <- pred %>% 
  filter(!is.na(value)) %>%
  spread(key = "type", value = "value")
# pred.wide$delta.fd <- pred.wide$fd - pred.wide$obs
# pred.wide$delta.nd <- pred.wide$nd - pred.wide$obs
pred.wide$delta.ia <- pred.wide$ia - pred.wide$obs
pred.wide$delta.id <- pred.wide$id - pred.wide$obs
pred.wide$delta.dd <- pred.wide$dd - pred.wide$obs

pred.wide.delta <- pred.wide[,c("model",
                                "year.predicted",
                                "plot",
                                "metric",
                                "intensity",
                                "timestep",
                                # "delta.fd",
                                # "delta.nd",
                                "delta.ia","delta.id","delta.dd")]
names(pred.wide.delta) <- c("model",
                            "year.predicted",
                            "plot",
                            "metric",
                            "intensity",
                            "timestep",
                            # "fd","nd",
                            "ia","id","dd")

pred.delta <- pred.wide.delta %>% gather(key = "type", value = "value",ia:dd)

# plot each combination
pd <- 0.2
ew <- 0.2
dk <- 0.2

for(i.model in 1:length(model_family)){
  for(i.timestep in 1:timesteps){
    
    pred.delta.avg <- pred.delta %>% 
      filter(model == model_family[i.model] & 
               timestep == i.timestep) %>%
      group_by(metric, intensity, type) %>%
      summarise(mean.value = mean(value,na.rm = TRUE), 
                sd.value = sd(value,na.rm = TRUE))
    
    pred.dif.plot <- ggplot(pred.delta.avg,aes(x = intensity, y = mean.value)) + 
      geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
      geom_errorbar(aes(ymin = mean.value-sd.value, ymax = mean.value+sd.value,
                        color = type),
                    position = position_dodge(pd),
                    width = ew)+
      geom_line(aes(color = type), position = position_dodge(pd))+
      geom_point(aes(fill = type), shape = 21, position = position_dodge(pd)) + 
      scale_fill_OkabeIto(darken = dk) + 
      scale_color_OkabeIto(darken = dk) + 
      facet_grid(metric~.,scales = "free_y") +
      theme_bw()+
      # scale_x_discrete(breaks=NULL)+
      theme(strip.background = element_blank())+
      ggtitle(paste(model_family[i.model]," - t",i.timestep,sep="")) +
      ylab("value")+
      NULL
    # pred.dif.plot
    
    ggsave(filename = paste("results/images/variations_",
                            model_family[i.model],"_t",
                            i.timestep,".pdf",sep=""),
           plot = pred.dif.plot,
           device = cairo_pdf,
           width = 8,height = 4,dpi = 300)
  }# for i.timestep
}# for i.model

# -------------------------------------------------------------------------
# 3 - observed-predicted

# abundance is summed up across subplots, whereas richness and evenness are avg

obs.subset <- subset(obs, year %in% unique(pred$year.predicted))
obs.wide <- pivot_wider(obs.subset,names_from = metric,values_from = value)
obs.wide.ag <- obs.wide %>% group_by(year,plot,type) %>%
  summarise(all.abund = sum(abundance),
            mean.rich = mean(richness),
            mean.ev = mean(evenness))
names(obs.wide.ag)[4:6] <- c("abundance","richness","evenness")
obs.ag <- pivot_longer(obs.wide.ag,cols = abundance:evenness,
                       names_to = "metric",values_to = "value")

obs.ag <- obs.ag[,c("year","plot","metric","value","type")]

pred.subset <- subset(pred,type == "obs")
pred.subset <- pred.subset %>%
  filter(!is.na(value),timestep == 2) %>%
  group_by(model,year.predicted,plot,metric) %>%
  summarise(mean.value = mean(value))
names(pred.subset) <- c("type","year","plot","metric","value")
pred.subset <- pred.subset[,c("year","plot","metric","value","type")]

obs.pred <- bind_rows(obs.ag,pred.subset)
 
obs.pred.plot <- ggplot(obs.pred, aes(x = type, y = value)) +
  geom_point(aes(fill = type), size = 1.5, shape = 21,
             position = position_jitterdodge()) +
  geom_boxplot(aes(fill = type),alpha = 0.60) +
  scale_fill_OkabeIto(darken = dk) +
  facet_grid(metric~.,scales = "free_y") +
  labs(x = NULL,y = NULL)+
  theme_bw()+
  scale_x_discrete(breaks=NULL)+
  theme(strip.background = element_blank())+
  NULL
obs.pred.plot

ggsave(filename = "results/images/observed_predicted_plot.pdf",
       plot = obs.pred.plot,
       device = cairo_pdf,
       width = 6,height = 4,dpi = 300)



