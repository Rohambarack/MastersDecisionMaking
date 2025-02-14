#take argument from bash
arg <- commandArgs(TRUE)
#make it an integer
seed <- as.integer(arg)
#simulate data based on a glm architecture
library(tidyverse)
library(brms)
pacman::p_load(LaplacesDemon)

source("u_functions.R")

n_runs <- 13
set.seed(seed)

#for loop 
for (i in 1:n_runs){
  start_t <- Sys.time()

Intercept_mu <- rnorm(1,0.645,0.3)
Intercept_sigma <- rgamma(1,4.5,35)# ~0.15
inter_mu <- rnorm(1,1,0.3)
inter_sigma <-rgamma(1,10,35) #~0.33
minor_mu <- rnorm(1,-1,0.3)
minor_sigma <- rgamma(1,10,35)

sim_params<- c(Intercept_mu,Intercept_sigma,
               inter_mu,inter_sigma,
               minor_mu,minor_sigma)
  
  
data_glm <- simulate_glm_data_non_diff(intercept = sim_params[1],
                                       sigma_ = sim_params[2],
                                       b_inter = sim_params[3],
                                       sigma_inter = sim_params[4],
                                       b_minor = sim_params[5],
                                       sigma_minor = sim_params[6]) 

#glm formula
f_decision <- brmsformula( y ~ 0 + induction_time*inter + induction_time*minor + ( 1 + inter + minor | gr(ID, by = induction_time)),
                           family = bernoulli())

#get_prior(f_decision,data_glm)
#priors -3 and 3 should cover most of the logit scale
p_decision <- c(
  prior(normal(0,1),class = b)
)

mod <- brm(formula = f_decision,
           data = data_glm,
           family = bernoulli(),
           prior = p_decision,
           sample_prior = T,
           cores = 4,
           chains = 4,
           iter = 4000,
           warmup = 2000)


##### make df

if (i ==1){
  
  output_df <- make_df_from_model_glm(mod,sim_params,i,seed)
  
}else{
  
  temp_df <- make_df_from_model_glm(mod,sim_params,i,seed)
  output_df <- rbind(output_df,temp_df)
  
}

end_t <- Sys.time()

print(end_t - start_t)
print("number:")
print(i)

}

####write output
name <- paste("seed",seed, sep = "_")
out_f <- paste("./results_glm/",name,".rds", sep = "")
write_rds(output_df,out_f)





