# plot rank-abundance curves
library(tidyverse)
library(colorblindr)

# from "extract_RAD_components.R"
rank.all <- read.csv2("./results/rank_abundances_complete.csv",header = TRUE,stringsAsFactors = FALSE)
rank.all$SAD.model <- factor(rank.all$SAD.model,levels = c("observed","niche_fitness","only.niche","only.fitness","increased.dispersal"))

# TEST
# rank.all <- subset(rank.all, SAD.model %in% c("niche_fitness","increased.dispersal"))

my.palette <- c("gray60","#009E73","#E69F00","#D55E00","darkblue")

rad.plot <- ggplot(rank.all,aes(x = sp.rank,y = rel.abund, group = SAD.model)) +
  geom_line(aes(color = SAD.model), size = 1.1) +#(aes(linetype = site.ID)) +
  # geom_line(aes(color = niche.apport), size = 1.1) +#(aes(linetype = site.ID)) +
  # geom_point(aes(fill = niche.apport),shape = 21, size = 1.5) +#(aes(shape = site.ID)) + 
  geom_point(aes(fill = SAD.model, shape = SAD.model), size = 2) +#(aes(shape = site.ID)) + 
  # facet_grid(trophic.guild~connectance+richness,scales = "free")+
  # facet_grid(fct_rev(trophic.guild)~connectance+resource.distribution,scales = "free")+
  facet_grid(plot~year,scales = "free")+
  scale_shape_manual(values = c(21,22,23,24,25))+
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


# abundances bar plot -----------------------------------------------------

abund.clean <- abund.obs[,c("year","plot","subplot","species","individuals")]
abund.clean$SAD.model <- "observed"
names(abund.clean)[4] <- "sp"

abund.all <- rbind(abund.clean,abund.proj)

plot.abund <- abund.all %>% 
  filter(!is.na(individuals)) %>%
  group_by(year,plot,sp,SAD.model) %>%
  summarise(sum.ind = sum(individuals))

# plot.list <- list()

years <- sort(unique(plot.abund$year))
plots <- sort(unique(plot.abund$plot))

for(i.year in 1:length(years)){
  for(i.plot in 1:length(plots)){
    my.abund <- subset(plot.abund,year == years[i.year] & plot == plots[i.plot])
    
    # only sp with observed + at least 3 predicted 
    my.sp <- sort(unique(my.abund$sp))
    to.remove <- NULL
    for(i.sp in 1:length(my.sp)){
      if(length(which(my.abund$sp == my.sp[i.sp])) <= 3){
        to.remove <- c(my.sp[i.sp],to.remove)
      }
    }
    
    my.abund <- subset(my.abund,!(sp %in% to.remove))
    
    my.plot <- ggplot(my.abund,aes(x = sp,y = sum.ind)) + 
      geom_col(aes(fill = SAD.model),position = position_dodge()) + 
      # facet_grid(plot~year,scales = "free")+
      xlab("") + ylab("abundance") +
      theme_bw()+
      scale_fill_manual(values = my.palette)+
      # guides(fill = FALSE)+
      theme(legend.title=element_blank())+
      theme(legend.justification=c(.95,.95),legend.position=c(.95,.95))+
      labs(x = NULL,
           y = "abundance",
           title = paste("year:",years[i.year],"plot:",plots[i.plot]))+
      # theme(plot.title.position = "plot")+
      theme(axis.text.x  = element_text(angle=90, vjust=0.5))+
      NULL
    
    ggsave(filename = paste("results/images/projected_abundances_",years[i.year],"_",plots[i.plot],".pdf",sep=""),
           plot = my.plot,device = cairo_pdf,width = 8,height = 4,dpi = 600)
  }
}

# all together
# bar.plot <- ggplot(plot.abund,aes(x = sp,y = sum.ind)) + 
#   geom_col(aes(fill = SAD.model),position = position_dodge()) + 
#   facet_grid(plot~year,scales = "free")+
#   xlab("species rank") + ylab("relative abundance") +
#   DGC::theme_Publication()+
#   scale_fill_manual(values = my.palette)+
#   theme(legend.title=element_blank())+
#   theme(strip.background = element_blank())+
#   NULL
# bar.plot

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

