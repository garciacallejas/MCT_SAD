# plot SAD metrics from simulations

library(tidyverse)
library(patchwork)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# -------------------------------------------------------------------------

sad.results <- read.csv2("results/SAD_metrics_silico.csv")

# remove outliers in steady state abundance
sad.results <- subset(sad.results, steady.state.abundance < 1e3)

# -------------------------------------------------------------------------

factors <- c("connectance","diagonal.dominance","mean.strength")
metrics <- c("steady.state.richness","steady.state.evenness","steady.state.abundance")

scatterplot.list <- list()

for(i.metric in 1:length(metrics)){
    for(i.factor in 1:length(factors)){
        my.data <- sad.results[,c(factors[i.factor],metrics[i.metric])]
        names(my.data) <- c("x","y")
        scatterplot.list[[length(scatterplot.list)+1]] <- ggplot(my.data,aes(x = x, y = y)) + 
            geom_jitter(width = .1, size = .5, alpha = .5) +
            # geom_point() +
            geom_smooth() +
            xlab(factors[i.factor]) + 
            ylab(metrics[i.metric]) +
            theme_bw() +
            NULL
    }# for i.metric
}# for i.factor

patchwork::wrap_plots(scatterplot.list,ncol = 3)




