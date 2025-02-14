#wrangle data to appropriate format
library(tidyverse)
library(brms)
data <- readRDS("../data/c_data.rds")
library(boot)

#TTB in minority list
minor_list <- seq(5,74,by = 3) - 2


data_glm <- data %>% 
  #adjusting so it doesn't start with turn 3
  mutate(turn = turn - 2,
         #adding the intervention 0s and 1s
         inter = ifelse((induction_time == 1 & turn > 23) | (induction_time == 2 & turn > 47),
                        1,
                        0),
         minor = ifelse(turn %in% minor_list,1,0))

#get important cols

data_glm <- data_glm %>% 
  select(ID,induction_time,turn,inter,minor,compare) %>% 
  rename( y = compare ) %>% 
  mutate( induction_time = as.factor(induction_time),
          inter= as.factor(inter),
          minor = as.factor(minor))



f_decision <- brmsformula( y ~ 0 + induction_time*inter + induction_time*minor + ( 1 + inter + minor | gr(ID, by = induction_time)),
                           family = bernoulli())

#get_prior(f_decision,data_glm)

p_decision <- c(
  prior(normal(0,1),class = b)
)

#get_prior(f_decision,data_glm)

#pp_1 <- brm(formula = f_decision,
#            data = data_glm,
#           family = bernoulli(),
#           prior = p_decision,
#           sample_prior = "only",
#           cores = 4,
#           chains = 4,
#            iter = 4000,
#            warmup = 2000)

#pp_check(pp_1, ndraws = 100)


mod <- brm(formula = f_decision,
           data = data_glm,
           family = bernoulli(),
           prior = p_decision,
           sample_prior = T,
           cores = 4,
           chains = 4,
           iter = 4000,
           warmup = 2000)

pp_check(mod, ndraws = 100)


mod
########################## mod check
#baseline
hypothesis(mod,"induction_time1 < induction_time2")
#some? definitive diff between the two baselines
#could indicate learning !

#intervention
hypothesis(mod,"induction_time1 < inter1 + induction_time1")
hypothesis(mod,"induction_time2 < inter1 + induction_time2 + induction_time2:inter1")
#intervention works for both
hypothesis(mod, "inter1 < inter1 + induction_time2:inter1")
#later intervention might have a higher impact, not definitive


#minority
hypothesis(mod,"induction_time1 > minor1 + induction_time1")
hypothesis(mod,"induction_time2 > minor1 + induction_time2 +induction_time2:minor1")
#minority works for both
hypothesis(mod, "minor1 > minor1 + induction_time2:minor1")
#later intervention lessen the impact of TTB in minority, not definitive


###############################
ests <- as_draws_df(mod)

i_t1 <- ests$b_induction_time1
i_t2 <- ests$b_induction_time2
i_t1_sd <- ests$`sd_ID__Intercept:induction_time1`
i_t2_sd <- ests$`sd_ID__Intercept:induction_time2`

b_inter_mean_t1 <- ests$b_inter1
b_inter_mean_t2 <- ests$b_inter1 + ests$`b_induction_time2:inter1`
b_inter_sd_t1 <- ests$`sd_ID__inter1:induction_time1`
b_inter_sd_t2 <- ests$`sd_ID__inter1:induction_time2`

b_minor_mean_t1 <- ests$b_minor1
b_minor_mean_t2 <- ests$b_minor1 + ests$`b_induction_time2:minor1`
b_minor_sd_t1 <- ests$`sd_ID__minor1:induction_time1`
b_minor_sd_t2 <- ests$`sd_ID__minor1:induction_time2`

### priors

b_priors <-ests$prior_b
sd_priors <- ests$prior_sd_ID

prior_df <- tibble(
  type = rep("prior",8000*6),
  param = rep(c("Intercept","inter_mean","minor_mean",
                "Intercept_sd","inter_sd","minor_sd"),each = 8000),
  vals = rep(c(b_priors,sd_priors),each = 3)
)

vis_ests <- tibble(
  
  type = rep(c("X1","X2"),each = 8000*6),
  param = rep(rep(c("Intercept","Intercept_sd",
                    "inter_mean","inter_sd",
                    "minor_mean","minor_sd"),each = 8000),2),
  vals = c(i_t1,i_t1_sd,
           b_inter_mean_t1,b_inter_sd_t1,
           b_minor_mean_t1,b_minor_sd_t1,
           #t2 vals
           i_t2,i_t2_sd,
           b_inter_mean_t2,b_inter_sd_t2,
           b_minor_mean_t2,b_minor_sd_t2)
) #%>% 
#mutate(vals = inv.logit(vals))


pp_u_df <- rbind(prior_df,vis_ests)
#pp_u
pp_u_df %>% 
  ggplot(aes(x= vals, fill = type)) +
  geom_density(alpha = .4) +
  ylab("")+
  xlab("")+
  scale_fill_manual(values = c("lightgreen","#cd7e3c","#3C8BCD"),
                    name = "Group") +
  theme_classic()+
  facet_wrap(~param, scales = "free", nrow = 3)


#visualize estimates
vis_ests %>% 
  ggplot(aes(x= vals, fill = type)) +
  geom_density(alpha = .4) +
  ylab("")+
  xlab("")+
  scale_fill_manual(values = c("#cd7e3c","#3C8BCD"),
                    name = "Group") +
  theme_classic()+
  facet_wrap(~param, scales = "free", nrow = 3)

#make results more interpretable
vis_ests_p <- tibble(
  
  type = rep(c("X1","X2"),each = 8000*4),
  param = rep(rep(c("Intercept",
                    "Intercept+m",
                    "Intercept+i",
                    "Intercept+i+m"),each = 8000),2),
  vals = c(i_t1,
           i_t1 + b_minor_mean_t1,
           i_t1 + b_inter_mean_t1,
           i_t1 + b_minor_mean_t1 + b_inter_mean_t1,
           #t2 vals
           i_t2,
           i_t2 + b_minor_mean_t2,
           i_t2 + b_inter_mean_t2,
           i_t2 + b_minor_mean_t2 + b_inter_mean_t2)
) %>% 
  mutate(vals = inv.logit(vals))

vis_ests_p %>% 
  ggplot(aes(x= vals, fill = type)) +
  geom_density(alpha = .4) +
  ggtitle("Chances of choosing TTB by predictor beta estimates") +
  ylab("")+
  xlab("")+
  scale_fill_manual(values = c("#cd7e3c","#3C8BCD"),
                    name = "Group") +
  theme_classic()+
  facet_wrap(~param, nrow = 2)

#get estimates in df
source("../u_functions.R")
glm_ests <- vis_ests %>% 
  group_by(type, param) %>% 
  reframe( est = MPD(vals),
           CI_10 = quantile(vals,c(0.1,0.9))[1],
           CI_90 = quantile(vals,c(0.1,0.9))[2])

write_rds(glm_ests,"glm_estimates.rds")


###################

pp_check(mod,ndraws = 100,)+
  ggtitle("Posterior Predictions - GLM")

##################

average_acc_glm <-function(model, data){
  seed_v <- 1:100
  acc_list <- array(NA,100)
  for (s in seed_v){
    set.seed(s)
    data_glm <- data
    mod <- model
    
    test <- predict(mod)
    ests <- ifelse(test[,1] > 0.5,1,0)
    data_glm$pred_y <- ests
    data_glm <- data_glm %>% 
      mutate(acc = ifelse(y == pred_y,1,0))
    
    acc <- sum(data_glm$acc)/nrow(data_glm)
    acc_list[s] <- acc
  }
  return(acc_list)
}

accuracy <- average_acc_glm(mod,data_glm)
write_rds(accuracy,"model_glm_acc.rds")