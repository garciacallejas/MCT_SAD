# plot SAD metrics from simulations

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
tau.mean <- unique(sad.metrics$tau)[5]

# metrics <- c("steady.state.richness","steady.state.evenness","steady.state.abundance")

# plot.list <- list()

rank.abundances <- sad.abundances %>% 
    select(-matrix.code) %>%
    group_by(connectance,diagonal.dominance,tau,replicate) %>% 
    mutate(relative.abund = abundance/sum(abundance)) %>%
    group_by(connectance,diagonal.dominance,tau,replicate) %>% 
    mutate(species.rank = rank(-relative.abund,ties.method = "first"))

# rank.abundances.long <- rank.abundances %>%
#     pivot_longer(connectance:tau,names_to = "metric",values_to = "value")

# -------------------------------------------------------------------------
plot.list <- list()

for(i.metric in c("connectance","diagonal.dominance","tau")){
    if(i.metric == "connectance"){
        my.metric <- subset(rank.abundances,
                            tau == tau.mean &
                                diagonal.dominance == diag.dom.mean)
        my.metric$metric <- as.factor(round(my.metric$connectance,2))
    }else if(i.metric == "diagonal.dominance"){
        my.metric <- subset(rank.abundances,
                            tau == tau.mean &
                                connectance == connectance.mean)
        my.metric$metric <- as.factor(round(my.metric$diagonal.dominance,2))
    }else if(i.metric == "tau"){
        my.metric <- subset(rank.abundances,
                            connectance == connectance.mean &
                                diagonal.dominance == diag.dom.mean)
        my.metric$metric <- as.factor(round(my.metric$tau,2))
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


ggsave(filename = "results/images/timestep_trends.pdf",
       plot = f1,
       device = cairo_pdf,
       width = 10,height = 10,dpi = 300)

# all sites
# rad.plot <- ggplot(rank.abundances,aes(x = species.rank,y = relative.abund, color = trophic.guild, group = site.ID)) +
#   geom_point() + 
#   geom_line() +
#   facet_grid(.~trophic.guild,scales = "free")+
#   scale_color_manual(values = my.palette)+
#   theme_Publication()+
#   theme(strip.background = element_blank(),strip.text.x = element_blank()) +
#   guides(color=FALSE)+
#   NULL
# rad.plot


# for(i.metric in 1:length(metrics)){
#     for(i.factor in 1:length(factors)){
#         my.data <- sad.results[,c(factors[i.factor],metrics[i.metric])]
#         names(my.data) <- c("x","y")
#         plot.list[[length(plot.list)+1]] <- ggplot(my.data,aes(x = x, y = y)) + 
#             geom_point() +
#             # geom_smooth() +
#             xlab(factors[i.factor]) + 
#             ylab(metrics[i.metric]) +
#             theme_bw() +
#             NULL
#     }# for i.metric
# }# for i.factor
# 
# patchwork::wrap_plots(plot.list,ncol = 3)




