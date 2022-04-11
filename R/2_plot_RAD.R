# plot SAD metrics from simulations and/or empirical data

library(tidyverse)
library(patchwork)

# -------------------------------------------------------------------------

sad.abundances <- read.csv2("results/SAD_abundances_silico.csv")
sad.metrics <- read.csv2("results/SAD_metrics_silico.csv")

# remove outliers in steady state abundance
# sad.results <- subset(sad.results, steady.state.abundance < 1e3)

# -------------------------------------------------------------------------
connectance.mean <- unique(sad.metrics$connectance)[5]
diag.dom.mean <- unique(sad.metrics$diagonal.dominance)[5]
strength.mean <- unique(sad.metrics$mean.strength)[5]

# metrics <- c("steady.state.richness","steady.state.evenness","steady.state.abundance")

# plot.list <- list()

rank.abundances <- sad.abundances %>% 
    select(-matrix.code) %>%
    group_by(connectance,diagonal.dominance,mean.strength,replicate) %>% 
    mutate(relative.abund = abundance/sum(abundance)) %>%
    group_by(connectance,diagonal.dominance,mean.strength,replicate) %>% 
    mutate(species.rank = rank(-relative.abund,ties.method = "first"))

# rank.abundances.long <- rank.abundances %>%
#     pivot_longer(connectance:tau,names_to = "metric",values_to = "value")

# -------------------------------------------------------------------------
plot.list <- list()

for(i.metric in c("connectance","diagonal.dominance","mean.strength")){
    if(i.metric == "connectance"){
        my.metric <- subset(rank.abundances,
                            mean.strength == strength.mean &
                                diagonal.dominance == diag.dom.mean)
        my.metric$metric <- as.factor(round(my.metric$connectance,2))
    }else if(i.metric == "diagonal.dominance"){
        my.metric <- subset(rank.abundances,
                            mean.strength == strength.mean &
                                connectance == connectance.mean)
        my.metric$metric <- as.factor(round(my.metric$diagonal.dominance,2))
    }else if(i.metric == "mean.strength"){
        my.metric <- subset(rank.abundances,
                            connectance == connectance.mean &
                                diagonal.dominance == diag.dom.mean)
        my.metric$metric <- as.factor(round(my.metric$mean.strength,2))
    }
    
    plot.list[[length(plot.list)+1]] <- ggplot(my.metric, aes(x = species.rank, 
                                                              y = relative.abund,
                                                              color = metric,
                                                              fill = metric,
                                                              group = interaction(metric,replicate))) +
        # geom_point(shape = 21) +
        geom_line(alpha = .5) +
        ggtitle(i.metric) +
        theme_bw() +
        theme(legend.title=element_blank()) +
        NULL
}

patchwork::wrap_plots(plot.list)


# -------------------------------------------------------------------------
# SCRIPT TO PLOT RADS from the Ph.D study

my.palette <- c("darkgreen","#009E73","#E69F00","#D55E00")


# abundances data
abundances <- readr::read_delim(file = "./data/all_abundances.csv",delim = ";",col_types = list(col_character(),col_character(),col_double()))
abundances <- subset(abundances, !is.na(abundance))

# species data
species.data <- readr::read_delim(file = "./data/all_species.csv",
                                  delim = ";",
                                  col_types = list(col_character(),
                                                   col_character(),
                                                   col_character(),
                                                   col_character(),
                                                   col_character()))

abundances$trophic.guild <- species.data$trophic.guild[match(abundances$species.ID,species.data$species.ID)]

# some cleaning
abundances <- subset(abundances,!is.na(trophic.guild))

abundances <- subset(abundances, abundance > 0)
abundances$trophic.guild <- factor(abundances$trophic.guild, levels = c("plants","herbivores","omnivores","carnivores"))


# -------------------------------------------------------------------------
# rank-abundance curves
rank.abund.empirical <- abundances[,c("site.ID","trophic.guild","abundance")]
rank.abund.empirical <- rank.abund.empirical %>% group_by(site.ID,trophic.guild) %>% mutate(relative.abund = abundance/sum(abundance))
rank.abund.empirical <- rank.abund.empirical %>% group_by(site.ID,trophic.guild) %>% mutate(species.rank = rank(-relative.abund,ties.method = "first"))

rank.abund.theor <- rank.abundances
rank.abund.theor$metric <- paste(rank.abund.theor$connectance,"_",
                                 rank.abund.theor$diagonal.dominance,"_",
                                 rank.abund.theor$mean.strength,sep="")
# all sites
rad.plot <- ggplot(rank.abund.empirical) +
  # geom_point() +
    geom_line(data = rank.abund.theor, aes(x = species.rank,y = relative.abund,
                                           group = interaction(metric,replicate)),
              color = "lightgrey",size = .7, alpha = .4) +
    geom_line(aes(x = species.rank,y = relative.abund, color = trophic.guild, group = site.ID),
              alpha = .5, size = .7) +
    xlim(c(0,25)) +
  # facet_grid(.~trophic.guild,scales = "free")+
  scale_color_manual(values = my.palette)+
  # theme_Publication()+
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  guides(color="none")+
  NULL
rad.plot

# # select some sites
# # trophic.guilds.per.site <- rank.abundances %>% count(site.ID,trophic.guild)
# # trophic.guilds.per.site <- trophic.guilds.per.site %>% group_by(site.ID) %>% summarise(num.guilds = n())
# # diverse.sites <- trophic.guilds.per.site$site.ID[trophic.guilds.per.site$num.guilds > 2]
# # 
# # # are the guilds well balanced?
# # diverse.sites.abund <- subset(abundances,site.ID %in% diverse.sites)
# # diverse.sites.abund <- diverse.sites.abund %>% group_by(site.ID,trophic.guild) %>% summarise(num.sp = n())
# # 
# # diverse.sites.plot <- ggplot(diverse.sites.abund, aes(x = trophic.guild, y = num.sp)) +
# #   geom_point() + 
# #   facet_wrap(~site.ID) +
# #   NULL
# # diverse.sites.plot
# # 
# # # sites with all guilds with >1 species
# # one.sp.sites <- unique(diverse.sites.abund$site.ID[diverse.sites.abund$num.sp == 1])
# # diverse.sites.abund <- subset(diverse.sites.abund,!(site.ID %in% one.sp.sites))
# # 
# # diverse.sites.plot <- ggplot(diverse.sites.abund, aes(x = trophic.guild, y = num.sp)) +
# #   geom_point() + 
# #   facet_wrap(~site.ID) + ylim(c(3,20)) +
# #   NULL
# # diverse.sites.plot
# 
# selected.sites <- c("ALTERDOC",1200,1271,1282,1822)
# 
# selected.abundances <- subset(rank.abundances,site.ID %in% selected.sites)
# selected.abundances$site.ID <- factor(selected.abundances$site.ID, levels = c("ALTERDOC",1200,1271,1282,1822))
# 
# selected.rad.plot <- ggplot(selected.abundances,aes(x = species.rank,y = relative.abund, group = trophic.guild)) +
#     geom_line(aes(color = trophic.guild), size = 1.1) +#(aes(linetype = site.ID)) +
#     geom_point(aes(fill = trophic.guild),shape = 21, size = 1.5) +#(aes(shape = site.ID)) + 
#     facet_grid(.~site.ID,scales = "free")+
#     scale_color_manual(values = my.palette)+#, labels = c("plants", "herbivores", "omnivores", "carnivores (inv)", "carnivores (vert)"))+
#     scale_fill_manual(values = my.palette)+#, labels = c("plants", "foli/granivores", "frugi/nectarivores", "omnivores", "carnivores (inv)", "carnivores (vert)"))+
#     #c("primary\nproducers", "plant/seed\neaters", "frugivores", "omnivores", "carnivores\n(invertebrates)", "carnivores\n(vertebrates)"))+
#     scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
#     xlab("species rank") + ylab("relative abundance") +
#     #DGC::theme_Publication()+
#     theme(strip.background = element_blank(),strip.text.x = element_blank()) +
#     #guides(color=FALSE)+#, fill = FALSE)+
#     NULL
# 
# tiff("./results/images/empirical_RAD_sites.tiff", res=600, compression = "lzw", width = 5500, height = 2500, units = "px")
# selected.rad.plot
# dev.off()

# -------------------------------------------------------------------------


# -------------------------------------------------------------------------


