library(tidyverse)
library(R2jags)

########### read in data and else
df_real <- read_rds("../c_data.rds")
minority <- readRDS("../data/TTB_minority.rds")
intervention_1 <- readRDS("../data/inter_1.rds")
intervention_2 <- readRDS("../data/inter_2.rds")
payoff <- readRDS("../data/experiment_payoff.rds")
#separate data by induction/intervention scheme
df_1 <- df_real %>% 
  filter(induction_time == 1)

df_2 <- df_real %>% 
  filter(induction_time == 2)

#find x for both df_1 and df_2
########## wrangle out rewards based on payoff and x
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}
make_x <- function(df){
  n_subj <- length(unique(df$ID))
  ntrials <- length(df$ID)/n_subj
  x_m <- vector("list", n_subj)
  
  for(i in 1:n_subj){
    #separate 1 subject at a time
    x_1_subj <- df$x[(((i-1)*ntrials)+1):(i*ntrials)]
    #make it a list of vectors
    x_m[[i]] <- x_1_subj
    
  }
  x_matrix <- do.call(rbind, x_m)
  
  return(x_matrix)
}
####################################################
df_1_x <- make_x(df_1)
df_2_x <- make_x(df_2)

df_1_subj <- length(unique(df_1$ID))
df_2_subj <- length(unique(df_2$ID))
#find r
df_1_r <- array(NA, c(df_1_subj,72))
df_2_r <- array(NA, c(df_2_subj,72))

for (s in 1:df_1_subj){
  for (t in 1:72){
    
    df_1_r[s,t] <- payoff[t, df_1_x[s,t]]
    
  }
}

for (s in 1:df_2_subj){
  for (t in 1:72){
    
    df_2_r[s,t] <- payoff[t, df_2_x[s,t]]
    
  }
}
##########################################
#model 1/3 intervention 
x <- df_1_x
r <- df_1_r
intervention <- intervention_1
ntrials <-72
n_subj <- df_1_subj
data <- list("x","r","intervention","minority","ntrials","n_subj") 
params<-c( "a_mean", "a_sd", 
           "beta_mean", "beta_sd",
           "inter_mean","inter_sd",
           "minor_mean","minor_sd",
           #get individual estimates for inspection
           "inter","minor","a","b")


start_t <- Sys.time() 
samples <- jags.parallel(data, inits=NULL, params,
                           model.file ="../bandit_model.txt",
                           n.chains=3,
                           n.iter=3000,
                           n.burnin=1000,
                           n.thin=1,
                           n.cluster = 3)

end_t <- Sys.time()
end_t-start_t

X_1 <- samples$BUGSoutput$sims.list

##########################################
#model 2/3 intervention 
x <- df_2_x
r <- df_2_r
intervention <- intervention_2
ntrials <-72
n_subj <- df_2_subj
data <- list("x","r","intervention","minority","ntrials","n_subj") 
params<-c( "a_mean", "a_sd", 
           "beta_mean", "beta_sd",
           "inter_mean","inter_sd",
           "minor_mean","minor_sd",
           #get individual estimates for inspection
           "inter","minor","a","b")


start_t <- Sys.time() 
samples_2 <- jags.parallel(data, inits=NULL, params,
                         model.file ="../bandit_model.txt",
                         n.chains=3,
                         n.iter=3000,
                         n.burnin=1000,
                         n.thin=1,
                         n.cluster = 3)

end_t <- Sys.time()
end_t-start_t

X_2 <- samples_2$BUGSoutput$sims.list


#prior-posterior update checks
#priors from the model spec
prior_a_mean <- runif(6000,0,1)
prior_a_sd <- rgamma(6000,3,3)
prior_beta_mean <- rgamma(6000,3,3)
prior_beta_sd <- rgamma(6000,3,15)
prior_inter_mean <- rnorm(6000,0,1)
prior_inter_sd <- rgamma(6000,4,1)
prior_minor_mean <- rnorm(6000,0,1)
prior_minor_sd <- rgamma(6000,4,1)

params <-c("a_mean","a_sd","beta_mean","beta_sd",
           "inter_mean","inter_sd","minor_mean","minor_sd")
types <- c("prior","X1","X2")

priors <- tibble(param = rep(params, each = 6000),
                 type = rep("prior", length(params)*6000),
                 values = c(prior_a_mean,
                            prior_a_sd,
                            prior_beta_mean,
                            prior_beta_sd,
                            prior_inter_mean,
                            prior_inter_sd,
                            prior_minor_mean,
                            prior_minor_sd))

X_1_estimates <- tibble(param = rep(params, each = 6000),
                        type = rep("X1", length(params)*6000),
                        values = c(X_1$a_mean,
                                   X_1$a_sd,
                                   X_1$beta_mean,
                                   X_1$beta_sd,
                                   X_1$inter_mean,
                                   X_1$inter_sd,
                                   X_1$minor_mean,
                                   X_1$minor_sd))

X_2_estimates <- tibble(param = rep(params, each = 6000),
                        type = rep("X2", length(params)*6000),
                        values = c(X_2$a_mean,
                                   X_2$a_sd,
                                   X_2$beta_mean,
                                   X_2$beta_sd,
                                   X_2$inter_mean,
                                   X_2$inter_sd,
                                   X_2$minor_mean,
                                   X_2$minor_sd))

vis_estimates <- rbind(priors,X_1_estimates,X_2_estimates)

vis_estimates %>% 
  ggplot(aes(x = values, fill = type))+
  geom_density(alpha=0.4) +
  facet_wrap(~param,scales = "free")+
  ylab("")+
  xlab("")+
  scale_fill_manual(values = c("lightgreen","#cd7e3c","#3C8BCD")) +
  theme_classic()

###########inspect inter and minor vals group level
inter_vect <- c()
minor_vect <- c()
inter_vect_2 <- c()
minor_vect_2 <- c()

for (i in 1:51){
subject_inter <- MPD(X_1$inter[,i])
subject_minor <- MPD(X_1$minor[,i])
minor_vect <- c(minor_vect,subject_minor)
inter_vect <- c(inter_vect, subject_inter)
}

for (i in 1:48){
  subject_inter <- MPD(X_2$inter[,i])
  subject_minor <- MPD(X_2$minor[,i])
  minor_vect_2 <- c(minor_vect_2,subject_minor)
  inter_vect_2 <- c(inter_vect_2, subject_inter)
}

ind_ests_1 <- tibble(
  param = rep(c("inter","minor"),each = 51),
  model = rep(c("X1"),  51*2),
  values = c(inter_vect,minor_vect)

)
ind_ests_2 <- tibble(
  param = rep(c("inter","minor"),each = 48),
  model = rep(c("X2"),  48*2),
  values = c(inter_vect_2,minor_vect_2)
  
)
ind_ests <- rbind(ind_ests_1,ind_ests_2)


ind_ests %>% 
ggplot(aes(x = values, fill = model))+
  geom_density(alpha=0.4)+
  ylab("")+
  xlab("")+
  scale_fill_manual(values = c("#cd7e3c","#3C8BCD")) +
  facet_wrap(~param,scales = "free")+
  theme_classic()

#################################
# compare population means

vis_estimates %>% 
  filter(str_detect(type,"X")) %>% 
  filter(str_detect(param,"mean|sd")) %>% 
  ggplot(aes(x = values, fill = type))+
  geom_density(alpha=0.4) +
  scale_fill_manual(values = c("#cd7e3c","#3C8BCD"),
                    name = "Group") +
  ylab("")+
  xlab("")+
  facet_wrap(~param,scales = "free")+
  theme_classic()

#get values estimates

ests <- vis_estimates %>% 
  filter(str_detect(type,"X")) %>% 
  group_by(param,type) %>% 
  reframe(Max_dens = MPD(values),
          CI_10 = quantile(values,c(0.1,0.9))[1],
          CI_90 = quantile(values,c(0.1,0.9))[2])
          


  
write_rds(ests,"estimates.rds")
  
  
  
  
  