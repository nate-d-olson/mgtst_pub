## MC Mix Specific counts
library(tidyverse)
pre_total <- 10000
post_total <- 10000
max_prop <- 1/pre_total

x <- runif(10000,min = 0, max = max_prop)
ggplot(data_frame(x)) + geom_density(aes(x = x))


mc <- x %>% map(~rbinom(n = 1000, size = 10000, prob = .))

mc_sim <- data_frame(prop = x, sim_count = mc) %>% unnest()

ggplot(mc_sim) + geom_density(aes(x = sim_count))

mc_sim$sim_count %>% summary()

quantile(mc_sim$sim_count,probs = c(0.95, 0.99))

sum(mc_sim$sim_count > 4.5)/10000000
sum(mc_sim$sim_count > 6.5)/10000000
length(mc_sim$sim_count)


sim_count <- function(min_prop, max_prop, total_abundance, obs_count){
      sim_prob <- runif(1, min = 0, max = max_prop) 
      sim_count <- rbinom(n = 1, size = total_abundance, prob = sim_prob) 
      list(sim_count = sim_count)
}

install.packages("MonteCarlo")
library(MonteCarlo)

# collect parameter grids in list:
param_list=list("min_prop" = 0, "max_prop" = max_prop, "total_abundance" = 50000, "obs_count" = 10)  

MC_result<-MonteCarlo(func=sim_count, nrep=10000, param_list=param_list)  

summary(MC_result)


feature_exp_df <- readRDS("../data/nb_expected_eo_metric_feature_df.rds") 
min_prop <- feature_exp_df %>% filter(pre != 0, post != 0) %>% 
      group_by(pipe, biosample_id, t_fctr) %>% 
      summarise(min_pre = min(pre),
                min_post = min(post)) %>% 
      mutate(theta = 2^-as.numeric(t_fctr),
             min_exp_prop = min_post*theta + min_pre * (1-0))  

eo_pos1 <- feature_exp_df  %>% filter(eo_metric == 1) %>% left_join(min_prop)

max_prop = eo_pos1[1,]$min_exp_prop
total_abu <- eo_pos1[1,]$total_abu 
obs_count <- eo_pos1[1,]$count

sim_count <- function(min_prop, max_prop, total_abundance, obs_count){
      sim_prob <- runif(1, min = 0, max = max_prop) 
      sim_count <- rbinom(n = 1, size = total_abundance, prob = sim_prob) 
      list(sim_count = sim_count)
}

library(MonteCarlo)

# collect parameter grids in list:
param_list=list("min_prop" = 0, "max_prop" = max_prop, "total_abundance" = total_abu, "obs_count" = obs_count)  

MC_result<-MonteCarlo(func=sim_count, nrep=100000, param_list=param_list)  

summary(MC_result)


## Wald Test 
mc_count <- MC_result$results$sim_count
mean_mc <- mean(mc_count)
var_mc <- var(mc_count)
wald_stat <- abs(mean_mc - obs_count)/sqrt(var_mc)
1-pchisq(wald_stat, 1)
qqnorm(mc_count)
### Using the likelihood ratio test
phat <- mean_mc/total_abu
loglik <- function(p, obs_count = 25, total_abu = 208790) dbinom(obs_count, total_abu, p, log = TRUE)

LRstats=2*(loglik(p = phat)-loglik(0))
LRstats


l_hat <- dbinom(obs_count, size = 208790,prob = mean_mc/total_abu, log = TRUE)
l_hat
l_obs <- dbinom(obs_count, size = 208790, prob = obs_count/total_abu, log = TRUE)
l_obs

LRstats=2*(l_hat - l_obs)

pchisq(LRstats, 1)
pvalue


LRstats=2*(loglik(0.48)-loglik(0.5))
LRstats

##### p-value from chisq distribution with degrees of freedom=1 

pvalue<-1-pchisq(LRstats, 1)
pvalue

##### p-value from chisq distribution with degrees of freedom=1 

pvalue<-1-pchisq(LRstats, 1)
pvalue

###### Likelihood ratio based 95% CI #####

obs_dist <- rbinom(n = 10000, size = total_abundance, prob = obs_count/total_abundance)
hist(obs_dist)

obs_count
#################################### Bayesian Approach #########################
# This is the best approach yet for obtaining a p-value but not sure it is better than just using the proportion of simulated counts >= to the observed counts 
## 
sim_count <- function(min_prop, max_prop, total_abundance){
      sim_prob <- runif(1, min = 0, max = max_prop) 
      sim_count <- rbinom(n = 1, size = total_abundance, prob = sim_prob) 
      list(sim_count = sim_count)
}


## sim distribution for P(Pi > max_prop)
param_list=list("min_prop" = max_prop, "max_prop" = 1, "total_abundance" = total_abu)  
mc_gt_max <-MonteCarlo(func=sim_count, nrep=10000, param_list=param_list)  
mc_gt_count <- mc_gt_max$results$sim_count
p_gt_max <- sum(mc_gt_count >= obs_count)/length(mc_gt_count)

## sim distribution for P(Pi <= max_prop)
param_list=list("min_prop" = 0, "max_prop" = max_prop, "total_abundance" = total_abu)  
mc_lt_max <-MonteCarlo(func=sim_count, nrep=10000, param_list=param_list)  
mc_lt_count <- mc_lt_max$results$sim_count
p_lt_max <- sum(mc_lt_count >= obs_count)/length(mc_lt_count)

## Bayesian Hypothesis test 

(p_lt_max * 0.5)/sum((p_lt_max * 0.5), (p_gt_max * 0.5))

################################################################################

#### Simulating pi
sim_prop <- function(obs_count, total_abundance){
      sim_count <- rbinom(n = 10000, 
                          size = total_abundance, 
                          prob = obs_count/total_abundance)
      sim_count/total_abundance
}

mc_prop <- sim_prop(obs_count, total_abu)

mc_prop %>% hist() 
## Wald Test
abs(mean(mc_prop) - max_prop)/sqrt(var(mc_prop))

## Likelihood Ratio  
# https://onlinecourses.science.psu.edu/stat504/node/39 