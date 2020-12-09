# plot SAD components
library(tidyverse)
library(colorblindr)
library(patchwork)

# read data ---------------------------------------------------------------

obs <- read.csv2(file = "results/observed_SAD_components_subplot.csv",
                 stringsAsFactors = FALSE)
pred <- read.csv2(file = "results/predicted_SAD_components_subplot_aggregated.csv",
                  stringsAsFactors = FALSE)

# plot --------------------------------------------------------------------
# 1 - observed-predicted

obs.subset <- subset(obs, year > 2015)

pred.subset <- subset(pred,type == "obs")
pred.subset <- pred.subset %>% 
  filter(!is.na(value)) %>%
  group_by(year.predicted,plot,subplot,metric) %>% 
  summarise(mean.value = mean(value))
names(pred.subset) <- c("year","plot","subplot","metric","value")
pred.subset$type <- "projected"

obs.pred <- bind_rows(obs.subset,pred.subset)

obs.pred.plot <- ggplot(obs.pred, aes(x = type, y = value)) + 
  geom_point(aes(fill = type), size = 1.5, shape = 21, 
             position = position_jitterdodge()) +
  geom_boxplot(aes(fill = type),alpha = 0.60) +
  scale_fill_OkabeIto() +
  facet_grid(metric~.,scales = "free_y") +
  labs(x = NULL,y = NULL)+
  theme_bw()+
  scale_x_discrete(breaks=NULL)+
  theme(strip.background = element_blank())+
  NULL
# obs.pred.plot

# 2 - variation against standard predicted
pred.wide <- pred %>% 
  filter(!is.na(value)) %>%
  spread(key = "type", value = "value")
pred.wide$delta.fd <- pred.wide$fd - pred.wide$obs
pred.wide$delta.nd <- pred.wide$nd - pred.wide$obs
pred.wide$delta.ia <- pred.wide$ia - pred.wide$obs
pred.wide$delta.id <- pred.wide$id - pred.wide$obs
pred.wide$delta.dd <- pred.wide$dd - pred.wide$obs

pred.wide.delta <- pred.wide[,c("year.predicted","plot","metric","intensity","delta.fd","delta.nd","delta.ia","delta.id","delta.dd")]
names(pred.wide.delta) <- c("year","plot","metric","intensity","fd","nd","ia","id","dd")

pred.delta <- pred.wide.delta %>% gather(key = "type", value = "value",fd:dd)

pred.delta.avg <- pred.delta %>% group_by(metric, intensity, type) %>%
  summarise(mean.value = mean(value,na.rm = TRUE), 
            sd.value = sd(value,na.rm = TRUE))

pred.dif.plot <- ggplot(pred.delta.avg,aes(x = intensity, y = mean.value)) + 
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_errorbar(aes(ymin = mean.value-sd.value, ymax = mean.value+sd.value,
                    color = type),
                position = position_dodge(.2),
                width = 0.2)+
  geom_line(aes(color = type), position = position_dodge(.2))+
  geom_point(aes(fill = type), shape = 21, position = position_dodge(.2)) + 
  scale_fill_OkabeIto(darken = 0.2) + 
  scale_color_OkabeIto(darken = 0.2) + 
  facet_grid(metric~.,scales = "free_y") +
  theme_bw()+
  # scale_x_discrete(breaks=NULL)+
  theme(strip.background = element_blank())+
  NULL
# pred.dif.plot

# store plots -------------------------------------------------------------

# ggsave(filename = "results/images/observed_predicted_subplot.pdf",
#        plot = obs.pred.plot,
#        device = cairo_pdf,
#        width = 6,height = 4,dpi = 300)
# 
# ggsave(filename = "results/images/predicted_perturbations_subplot.pdf",
#        plot = pred.dif.plot,
#        device = cairo_pdf,
#        width = 8,height = 4,dpi = 300)


