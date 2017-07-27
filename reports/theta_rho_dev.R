library(tidyverse)
x <- 1:5
theta <- 2^-x
theta
rho_pre <- rep(sample(seq(0,0.75,0.01), 10), each = 10)
rho_post <- rep(sample(seq(0,0.75,0.01), 10), 10)
dat <- data_frame(x, theta, rho_pre = list(rho_pre), rho_post = list(rho_post)) %>% 
    add_column(sim = 1:nrow(.)) %>% group_by(sim) %>% nest()

n_features <- 10
pre <- sample(seq(to = 0.1, from = 0.001, by = 0.001), size = n_features)
post <- sample(seq(to = 0.1, from = 0.001, by = 0.001), size = n_features)

sim_count <- data_frame(sim = 1:nrow(dat), pre_f_prob = list(pre), post_f_prob = list(post)) %>%
    unnest() %>% rowwise() %>% 
    mutate(pre_f_count = rbinom(n = 1, size = 10000, prob = pre_f_prob),
           post_f_count = rbinom(n = 1, size = 10000, prob = post_f_prob))

sim_dat <-sim_count %>% left_join(dat) %>% unnest() %>% unnest()

## calculate ng/ul DNA for pre and post
sim_dat <- sim_dat %>% mutate(pre_ng = 12.5*rho_pre*(1-theta), 
                              post_ng = 12.5*rho_post*theta,
                              rho_theta = post_ng/(pre_ng + post_ng))

sim_dat %>% ggplot() + 
    geom_point(aes(x = factor(rho_post), 
                   y = rho_theta, color = factor(rho_pre))) + 
    geom_hline(aes(yintercept = theta)) + 
    facet_wrap(~theta) + theme(axis.text.x = element_text(angle = 90))

sim_dat %>% group_by(sim, theta, rho_pre, rho_post) %>% nest() #%>% 
    # mutate(fit = map(data, ~lm()))
    
## Simulate counts using rho 
## Fit theta estimate model to see how well it agrees with rho_theta  


