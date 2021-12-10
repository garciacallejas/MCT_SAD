# plot SAD components
library(tidyverse)
library(colorblindr)
library(patchwork)

# read data ---------------------------------------------------------------

pred <- read.csv2(file = "results/predicted_SAD_components_subplot_v2.csv",
                  stringsAsFactors = FALSE)

pred$type <- recode(pred$type, obs = "observed", dd = "diagonal_dominance",
                    ia = "interaction_asymmetry",
                    id = "interaction_strength_distribution")

# TODO check discrepancies between dist.plot and pred.dif.plot

# plots -------------------------------------------------------------------

save.plots <- FALSE

# -------------------------------------------------------------------------
# 1 - overall distributions

my.intensities <- c(2,5,10)

pred.wide <- pivot_wider(pred,names_from = metric,values_from = value)
pred.wide$log.abundance <- log(pred.wide$abundance)
pred.2 <- pivot_longer(pred.wide,abundance:log.abundance,
                       names_to = "metric",
                       values_to = "value")
pred.3 <- filter(pred.2,intensity %in% my.intensities & !is.na(value)) %>%
  filter(metric != "abundance")
  

dist.plot <- ggplot(pred.3,aes(x = type, y = value)) +
  # geom_jitter(aes(fill = type), shape = 21, alpha = .2, size = .2) +
  geom_boxplot(aes(fill = type), alpha = .5) +
  facet_grid(metric~intensity, scales = "free_y") +
  theme_bw()+
  scale_fill_OkabeIto(darken = 0.2) +
  scale_x_discrete(breaks=NULL)+
  theme(strip.background = element_blank())+
  NULL
# dist.plot

# -------------------------------------------------------------------------
# frequency of stable communities

pred.count <- pred %>% 
  filter(metric == "abundance") %>%
  group_by(type,intensity) %>%
  summarise(num.stable = sum(!is.na(value)))

st.plot <- ggplot(pred.count,aes(x = intensity, y = num.stable, group = type)) + 
  geom_line(aes(color = type)) +
  geom_point(aes(fill = type), shape = 21, size = 2) +
  scale_color_OkabeIto(darken = 0.2) +
  scale_fill_OkabeIto(darken = 0.2) +
  theme_bw()+
  NULL
# st.plot

# -------------------------------------------------------------------------
# 2 - variations with respect to the standard matrix

pred.wide <- pred %>%
  filter(!is.na(value)) %>%
  pivot_wider(names_from = type,values_from = value)

# pred.wide$delta.fd <- pred.wide$fd - pred.wide$obs
# pred.wide$delta.nd <- pred.wide$nd - pred.wide$obs
pred.wide$delta.ia <- pred.wide$interaction_asymmetry - pred.wide$observed
pred.wide$delta.id <- pred.wide$interaction_strength_distribution - pred.wide$observed
pred.wide$delta.dd <- pred.wide$diagonal_dominance - pred.wide$observed

pred.wide.delta <- pred.wide[,c("year",
                                "plot",
                                "subplot",
                                "metric",
                                "intensity",
                                # "delta.fd",
                                # "delta.nd",
                                "delta.ia","delta.id","delta.dd")]
names(pred.wide.delta) <- c("year",
                            "plot",
                            "subplot",
                            "metric",
                            "intensity",
                            # "delta.fd",
                            # "delta.nd",
                            "interaction_asymmetry",
                            "interaction_strength_distribution",
                            "diagonal_dominance")

pred.delta <- pred.wide.delta %>% gather(key = "type", 
                                         value = "value",
                                         interaction_asymmetry:diagonal_dominance)

# plot 
pd <- 0.2
ew <- 0.2
dk <- 0.2
    
pred.delta.avg <- pred.delta %>% 
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
  ylab("value")+
  NULL
# pred.dif.plot
    
# save plots --------------------------------------------------------------

if(save.plots){
  ggsave(filename = "results/images/SAD_components_distributions.pdf",
         plot = dist.plot,
         device = cairo_pdf,
         width = 10,height = 10,dpi = 300)
  ggsave(filename = "results/images/number_stable_communities.pdf",
         plot = st.plot,
         device = cairo_pdf,
         width = 6,height = 4,dpi = 300)
  ggsave(filename = "results/images/variation_in_SAD_components.pdf",
         plot = pred.dif.plot,
         device = cairo_pdf,
         width = 8,height = 10,dpi = 300)
}

