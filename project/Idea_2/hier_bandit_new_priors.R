#take argument from bash
arg <- commandArgs(TRUE)
#make it an integer
seed <- as.integer(arg)
# setting up a bandit-model
#install.packages("pacman")
pacman::p_load(tidyverse,hesim, extraDistr, R2jags, parallel,truncnorm)
# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#import payoff structure
payoff <- readRDS("./data/experiment_payoff.rds")
minority <- readRDS("./data/TTB_minority.rds")
intervention_1 <- readRDS("./data/inter_1.rds")
intervention_2 <- readRDS("./data/inter_2.rds")

inter_scheme <- list(intervention_1,intervention_2)
#get sim functions
source("Simulate_bandit.R")


set.seed(seed)
n_runs <- 1:24  #12
ntrials <- 72
n_subj <-  50 #50

for (i in n_runs){
  for (i_time in 1:length(inter_scheme)){
# values for simulations
a_mean <- runif(1,0,0.2)#.1
a_sd <- rgamma(1,4,1) #0.01
beta_mean <- rgamma(1,3,3)#.05 # higher number means more consistent choice behavior (aka. less exploration)
beta_sd <- rgamma(1,1.5,15) #0.1
  
inter_mean <- rnorm(1,0.6,0.2) #0.2
inter_sd <- rgamma(1,4,1) #0.05

minor_mean <- rgamma(1,0.6,0.2) #0.3
minor_sd <- rgamma(1,4,1) #0.1

RW_sims <- RW_hier_inter_minor(payoff = payoff,
                               intervention = inter_scheme[[1]],
                               minority = minority,
                               ntrials = ntrials,
                               n_subj = n_subj,
                               inter_mean = inter_mean,
                               inter_sd = inter_sd,
                               minor_mean = minor_mean,
                               minor_sd = minor_sd,
                               a_mean = a_mean,
                               a_sd = a_sd,
                               beta_mean = beta_mean,
                               beta_sd = beta_sd
                         
)


test <- tibble( arm = rep(c("TTB","non_TTB"), each = 72),
                turn = c(1:72,1:72),
                Q = c(RW_sims$Q[1,,1],
                          RW_sims$Q[1,,2]))
#agent assessment 
test %>% 
  ggplot(aes(y = Q, x = turn, colour = arm))+
  geom_point(alpha = .7)+
  #geom_smooth(method = "lm") +
  geom_vline(xintercept = 25, col = "red") +
  scale_color_manual(values = c("#cd7e3c","#3C8BCD")) +
  #facet_wrap(~arm) +
  theme_classic()
########################## plots for priors
# curve
#curve(dgamma(x,1,15), from = 0, to = 0.5)
#curve(dgamma(x,2,15), from = 0, to = 0.5, add = T)
#curve(dgamma(x,3,15), from = 0, to = 0.5, add = T)
############################################# 

x <- RW_sims$x
r <- RW_sims$r
intervention <- inter_scheme[[i_time]]
data <- list("x","r","intervention","minority","ntrials","n_subj") 
params<-c( "a_mean", "a_sd", 
          "beta_mean", "beta_sd",
          "inter_mean","inter_sd",
          "minor_mean","minor_sd")


start_t <- Sys.time() 
samples_2 <- jags.parallel(data, inits=NULL, params,
                           model.file ="bandit_model.txt",
                           n.chains=3,
                           n.iter=4000,
                           n.burnin=2000,
                           n.thin=1,
                           n.cluster = 3)


X <- samples_2$BUGSoutput$sims.list

if(i==1 && i_time == 1){
  
  output_df <- make_df_from_model(X,RW_sims,i,i_time,seed)
  
}else{
  
  temp_df <- make_df_from_model(X,RW_sims,i,i_time,seed)
  output_df <- rbind(output_df,temp_df)
  
}

end_t <- Sys.time()

print(end_t - start_t)


  }
  
  print("number:")
  print(i)
  
}

#output_df %>% 
#  filter(i_time == 2) %>% 
#  ggplot(aes(x=t_val, y =i_val))+
#  geom_point()+
#  geom_abline(intercept = 0, slope = 1, color = "red")+
#  geom_smooth(method = "lm")+
#  ylab("Inferred")+
#  xlab("True Value") +
#  facet_wrap(~param, scales = "free")

#write output to folder
name <- paste("seed",seed, sep = "_")
out_f <- paste("./results_np/",name,".rds", sep = "")
write_rds(output_df,out_f)
