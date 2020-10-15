# plot rad components
library(tidyverse)
library(colorblindr)

rad.components <- read.csv2("./results/rad_components.csv",header = TRUE,stringsAsFactors = FALSE)

rad.components$SAD.model <- factor(rad.components$SAD.model,levels = c("observed","niche_fitness","increased.dispersal","only.niche","only.fitness"))

components.plot <- ggplot(rad.components,aes(x = SAD.model, y = value)) +
  geom_point(aes(fill = SAD.model), size = 1.5, shape = 21, position = position_jitterdodge()) +
  geom_boxplot(aes(fill = SAD.model),alpha = 0.60) +
  scale_fill_OkabeIto() +
  facet_grid(component~.,scales = "free_y") +
  labs(x = NULL,y = NULL)+
  theme_bw()+
  scale_x_discrete(breaks=NULL)+
  theme(strip.background = element_blank())+
  NULL
components.plot

ggsave("./results/images/SAD_components.pdf",components.plot,width = 5,height = 8,dpi = 600)


