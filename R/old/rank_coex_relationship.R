# relate SAD rank to species coexistence
library(tidyverse)

# read data ---------------------------------------------------------------

rank.data <- read.csv2("./results/rank_abundances_predicted.csv",
                       header = TRUE,
                       stringsAsFactors = FALSE)
obs.rank.data <- read.csv2("./results/rank_abundances_observed.csv",
                           header = TRUE,
                           stringsAsFactors = FALSE)
# species roles from spatial_coexistence study
# I can use them because we use the same Ricker model to infer
# interaction matrices (although, if picky, I should not have
# used the annual plant version to derive SADs here)

coex.data <- read.csv2("../Spatial_coexistence/results/observed_species_roles.csv",
                       header = TRUE,
                       stringsAsFactors = FALSE)

# TODO infer species roles for only niche and only fitness models
# for now, stick only with "niche_fitness" model
rank.data <- subset(rank.data, SAD.model == "niche_fitness")

# join rank and coex data -------------------------------------------------

# coex.rank <- left_join(rank.data,coex.data)
coex.rank <- left_join(obs.rank.data,coex.data)

# make pairwise-indirect dichotomic?
coex.rank$indirect <- ifelse(coex.rank$pairwise == TRUE,FALSE,coex.rank$indirect)

coex.rank.long <- gather(coex.rank,key = "coexistence_type",value = "value",pairwise:transient)
coex.rank.long <- subset(coex.rank.long,value == TRUE)
coex.rank.long$coexistence_type <- factor(coex.rank.long$coexistence_type,levels = c("pairwise","indirect","transient","dominant"))

# plot --------------------------------------------------------------------

my.palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
palette.values <- my.palette[c(2,6,4,3,7)]

coex.rank.boxplot <- ggplot(coex.rank.long,aes(x = coexistence_type, y = sp.rank))+
  geom_boxplot() + 
  NULL
coex.rank.boxplot

# average coex-rank plot --------------------------------------------------

avg.coex.rank <- coex.rank %>% group_by(sp) %>% summarise(avg.rank = mean(sp.rank),prob.coex = sum(pairwise | indirect)/n())

avg.coex.rank.plot <- ggplot(avg.coex.rank,aes(x = avg.rank,y = prob.coex, label = sp))+
  geom_point(aes(color = sp))+
  geom_smooth(method = "lm")+
  geom_text(aes(color = sp),vjust = 0, nudge_y = 0.01)+
  NULL
avg.coex.rank.plot


