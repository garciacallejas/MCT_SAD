# plot rank-abundance curves
library(tidyverse)

abund.proj <- read.csv2("./results/projected_abundances_subplot.csv",header = TRUE,stringsAsFactors = FALSE)
abund.obs <- read.csv2("../Caracoles/data/abundances.csv",header = TRUE,stringsAsFactors = FALSE)

# generate rank-abundance data
# 1 - observed
rank.obs <- abund.obs %>%
  filter(!is.na(individuals) &
           individuals > 0) %>%
  group_by(year,plot,species) %>%
  summarise(mean.abund = mean(individuals)) %>%
  group_by(year,plot) %>%
  mutate(plot.abund = sum(mean.abund)) %>%
  group_by(year,plot,species) %>%
  summarise(rel.abund = mean.abund/plot.abund) %>%
  mutate(sp.rank = rank(-rel.abund,ties.method = "first"))

rank.obs <- arrange(rank.obs,year,plot,desc(rel.abund))
rank.obs$SAD.model <- "observed"
names(rank.obs)[which(names(rank.obs) == "species")] <- "sp"
# write.csv2(rank.obs,file = "./results/rank_abundances_observed.csv",row.names = FALSE)

# 2 - projected
rank.abundances <- abund.proj %>% 
  filter(!is.na(individuals)) %>%
  group_by(year,plot,SAD.model,sp) %>% 
  summarise(mean.abund = mean(individuals)) %>%
  group_by(year,plot,SAD.model) %>%
  mutate(plot.abund = sum(mean.abund)) %>%
  group_by(year,plot,SAD.model,sp) %>%
  summarise(rel.abund = mean.abund/plot.abund) %>%
  mutate(sp.rank = rank(-rel.abund,ties.method = "first"))

rank.abundances <- arrange(rank.abundances,SAD.model,year,plot,desc(rel.abund))
# write.csv2(rank.abundances,file = "./results/rank_abundances_predicted.csv",row.names = FALSE)

# plot
rank.all <- bind_rows(rank.obs,rank.abundances)
my.palette <- c("gray60","#009E73","#E69F00","#D55E00")

rad.plot <- ggplot(rank.all,aes(x = sp.rank,y = rel.abund, group = SAD.model)) +
  geom_line(aes(color = SAD.model), size = 1.1) +#(aes(linetype = site.ID)) +
  # geom_line(aes(color = niche.apport), size = 1.1) +#(aes(linetype = site.ID)) +
  # geom_point(aes(fill = niche.apport),shape = 21, size = 1.5) +#(aes(shape = site.ID)) + 
  geom_point(aes(fill = SAD.model, shape = SAD.model), size = 2) +#(aes(shape = site.ID)) + 
  # facet_grid(trophic.guild~connectance+richness,scales = "free")+
  # facet_grid(fct_rev(trophic.guild)~connectance+resource.distribution,scales = "free")+
  facet_grid(plot~year,scales = "free")+
  scale_shape_manual(values = c(21,22,23,24))+
  scale_color_manual(values = my.palette)+
  scale_fill_manual(values = my.palette)+
  # xlim(0,20)+
  # scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
  xlab("species rank") + ylab("relative abundance") +
  DGC::theme_Publication()+
  theme(legend.title=element_blank())+
  theme(strip.background = element_blank())+#,strip.text.x = element_blank()) +
  #guides(color=FALSE)+#, fill = FALSE)+
  NULL
# rad.plot

# ggsave(filename = paste("results/images/SAD_plot_year_negbin.pdf",sep=""),plot = rad.plot,device = cairo_pdf,width = 12,height = 12,dpi = 600)

# are species ranks well predicted?
# rank.obs.wide <- spread(rank.obs,key = SAD.model,value = sp.rank)
# names(rank.obs.wide)[5] <- "observed_rank"
# names(rank.obs.wide)[4] <- "observed_rel_abund"
# 
# rank.proj <- rank.abundances
# names(rank.proj)[5] <- "predicted_rel_abund"
# names(rank.proj)[6] <- "predicted_rank"
# 
# rank.data <- left_join(rank.proj,rank.obs.wide)
# rank.data <- arrange(rank.data,year,plot,sp,SAD.model)
# 
# rank.plot <- ggplot(rank.data, aes(x = observed_rank,y = predicted_rank,group = sp)) +
#   geom_point(aes(color = sp)) +
#   geom_abline(slope = 1,linetype = "dashed",color = "grey")+
#   facet_grid(year~SAD.model)+
#   NULL
# rank.plot
# 
# rel.abund.plot <- ggplot(rank.data, aes(x = observed_rel_abund,y = predicted_rel_abund,group = sp)) +
#   geom_point(aes(color = sp)) +
#   geom_abline(slope = 1,linetype = "dashed",color = "grey")+
#   facet_grid(year~SAD.model)+
#   NULL
# rel.abund.plot

