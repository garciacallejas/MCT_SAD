# returns a list of matrices with decreasing kurtosis in their elements
# or increasing kurtosis if decreasing = FALSE

# NOTE: lower tau values indicate HIGHER kurtosis, and viceversa
# and the "decreasing" argument refers to the kurtosis, not to the tau values

gradient_strength_dist <- function(A,steps,
                                   # min.value = 0.01,
                                   # max.value = 0.99,
                                   initial.tau = NULL,
                                   min.tau = 0.5,
                                   max.tau = 1.5,
                                   decreasing = TRUE){
  
  if(decreasing){
    # decreasing kurtosis -> sampling tau from initial/min to max
    if(!is.null(initial.tau)){
      my.max.tau = max.tau
      my.min.tau = initial.tau
    }else{
      my.max.tau = max.tau
      my.min.tau = min.tau
    }
  }else{
    # increasing kurtosis -> sampling tau from initial/max to min
    if(!is.null(initial.tau)){
      my.max.tau = min.tau
      my.min.tau = initial.tau
    }else{
      my.max.tau = min.tau
      my.min.tau = max.tau
    }
  }
  
  # all A values are changed here
  # this is a test
  # a.rows <- 50
  # a.cols <- 50
  # c <- 0.2
  # l <- c * (a.rows*a.cols)
  
  tau.steps <- seq(from = min.tau,
                   to = max.tau,
                   length.out = steps)
  res <- list()
  for(i.tau in 1:tau.steps){
    A <- matrix(0,nrow = a.rows,ncol = a.cols)
    ints <- abs(gamlss.dist::rSHASHo(l, mu = 0, 
                                     sigma = 1, nu = 0, tau = tau.steps[i.tau]))
    # hist(ints,breaks = 20)
    for(i in 1:l){
      my.sample.row <- sample(1:a.rows,1,replace = T)
      my.sample.col <- sample(1:a.cols,1,replace = T)
      
      while(A[my.sample.row,my.sample.col] != 0 & 
            my.sample.row == my.sample.col){
        my.sample.row <- sample(1:a.rows,1,replace = T)
        my.sample.col <- sample(1:a.cols,1,replace = T)
      }
      A[my.sample.row,my.sample.col] <- ints[i]
    }

    A.list[[i.tau]] <- A
    
  }# for i.tau
  
  return(res)
  
  ## previous version
  # for(i.step in 1:steps){
  #   step.matrix <- A
  #   
  #   for(i in 1:nrow(A)){
  #     for(j in 1:ncol(A)){
  #       my.gradient <- seq(from = A[i,j],to = mean.interaction,
  #                          length.out = steps)
  #       step.matrix[i,j] <- my.gradient[i.step]
  #     }
  #   }
  #   
  #   res[[i.step]] <- step.matrix
  #   
  # }# for i.step
  
}

library(tidyverse)
library(moments)
x <- 1e5
tau_list <- c(0.25, 0.5, 0.75, 1, 1.5)
df <- data.frame()
for (tau in tau_list) {
  temp_df <- data.frame(x = x,
                        # y = gamlss.dist::dSHASHo(x, mu = 0, sigma = 1, nu = 0, tau = tau))
                        y = abs(gamlss.dist::rSHASHo(x, mu = 0, sigma = 1, nu = 0, tau = tau)))

  temp_df$tau <- tau
  df <- rbind(df, temp_df)
}
ggplot(data = df) +
  geom_histogram(aes(color = factor(tau), y = y)) +
  facet_grid(factor(tau)~.,scales = "free_y") +
  theme_bw()


df %>%
  group_by(tau) %>%
  summarise(kurtosis = kurtosis(y),max = max(y),min = min(y))




