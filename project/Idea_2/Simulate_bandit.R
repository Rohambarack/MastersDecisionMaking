#modified rescorla Wagner with group level effects for 2 bandits
RW_hier <- function(payoff, ntrials, n_subj,
                          a_mean, a_sd, 
                          beta_mean, beta_sd) {
  
  a <- array(NA, c(n_subj))
  beta <- array(NA, c(n_subj))
  x <- array(NA, c(n_subj,ntrials))
  r <- array(NA, c(n_subj,ntrials))
  Q <- array(NA, c(n_subj,ntrials,2))
  # bias
  bias <- array(NA, c(n_subj,ntrials,2))
  Qupdate <- array(NA, c(n_subj,ntrials, 2))
  exp_p <- array(NA, c(n_subj,ntrials, 2))
  p <- array(NA, c(n_subj,ntrials, 2))
  
  for (s in 1:n_subj){
    #define participant a and beta
    a[s] <- rtruncnorm(1,0.001,1,a_mean,a_sd)
    beta[s] <- rtruncnorm(1,0.001,,beta_mean,beta_sd)
    
    # set value for trial 1 roughly based on exp
    # TTB choice is 5 out of 6
    # not taking is adding up the remaining two, so slightly worse
    Q[s,1,1] <- rbeta(1,5,1)
    Q[s,1,2] <- rbeta(1,7,5)
    ################manually set x on first trial
    for (k in 1:2){
      exp_p[s,1,k] <- exp(beta[s]*Q[s,1,k])
    }
    for (k in 1:2) {
      
      p[s,1,k] <- exp_p[s,1,k]/sum(exp_p[s,1,])
      
    }
    
    x[s,1] <- rcat(1, p[s,1,])
    
    r[s,1] <- payoff[1, x[s,1]]
    ##########################################
    for (t in 2:ntrials) {
      
      for (k in 1:2) {
        
        Qupdate[s,t,k] <- Q[s,t-1, k] + (a[s]*(r[s,t-1]-Q[s,t-1,k]))
        Q[s,t,k] <- ifelse(k==x[s,t-1],
                           Qupdate[s,t,k],
                           Q[s,t-1,k])
        
        exp_p[s,t,k] <- exp(beta[s]*Q[s,t,k])
        
      }
      
      for (k in 1:2) {
        
        p[s,t,k] <- exp_p[s,t,k]/sum(exp_p[s,t,])
        
      }
      
      x[s,t] <- rcat(1, p[s,t,])
      
      r[s,t] <- payoff[t, x[s,t]]
      
    }
    
  }
  
  result <- list(x=x,
                 r=r,
                 Q=Q,
                 p=p,
                 a=a,
                 beta=beta)
  return(result)
  
}

#modified rescorla Wagner with group level effects for 2 bandits
# and group level intervention effect, and truth in minority
RW_hier_inter_minor <- function(payoff, intervention, minority,
                    ntrials, n_subj,
                    inter_mean, inter_sd,
                    minor_mean, minor_sd,
                    a_mean, a_sd, 
                    beta_mean, beta_sd) {
  
  a <- array(NA, c(n_subj))
  beta <- array(NA, c(n_subj))
  inter <- array(NA, c(n_subj))
  minor <- array(NA, c(n_subj))
  x <- array(NA, c(n_subj,ntrials))
  r <- array(NA, c(n_subj,ntrials))
  Q <- array(NA, c(n_subj,ntrials,2))
  Qupdate <- array(NA, c(n_subj,ntrials, 2))
  exp_p <- array(NA, c(n_subj,ntrials, 2))
  p <- array(NA, c(n_subj,ntrials, 2))
  
  for (s in 1:n_subj){
    #define participant a and beta, inter and min effect
    a[s] <- rtruncnorm(1,0.001,1,a_mean,a_sd)
    beta[s] <- rtruncnorm(1,0.001,,beta_mean,beta_sd)
    
    inter[s] <- rnorm(1,inter_mean,inter_sd) 
    minor[s] <- rnorm(1,minor_mean,minor_sd)
    
    # set value for trial 1 roughly based on exp
    # TTB choice is 5 out of 6
    # not taking is adding up the remaining two, so slightly worse
    # not adding the effect directly to Q because it shouldn't influence the
    # next turn
    Q[s,1,1] <- rbeta(1,5,1) 
    Q[s,1,2] <- rbeta(1,7,5)
    ################manually set x on first trial
    for (k in 1:2){
      
      # add intervention effect to TTB arm, minority to non_TTB if the
      # conditions apply
      exp_p[s,1,k] <- exp(beta[s]*
                            ifelse(k==1,
                                   Q[s,1,k] * (intervention[1]*inter[s]),
                                   Q[s,1,k] * (minority[1]*minor[s]))
                          )
    }
    for (k in 1:2) {
      
      p[s,1,k] <- exp_p[s,1,k]/sum(exp_p[s,1,])
      
    }
    
    x[s,1] <- rcat(1, p[s,1,])
    
    r[s,1] <- payoff[1, x[s,1]]
    ##########################################
    for (t in 2:ntrials) {
      
      for (k in 1:2) {
        #learning component
        Qupdate[s,t,k] <- Q[s,t-1, k] + (a[s]*(r[s,t-1]-Q[s,t-1,k]))
        Q[s,t,k] <- ifelse(k==x[s,t-1],
                           Qupdate[s,t,k],
                           Q[s,t-1,k])
        
        exp_p[s,t,k] <- exp(beta[s]*
                              ifelse(k==1,
                                     Q[s,t,k] * (intervention[t]*inter[s]),
                                     Q[s,t,k] * (minority[t]*minor[s]))
        )
        
      }
      
      for (k in 1:2) {
        
        p[s,t,k] <- exp_p[s,t,k]/sum(exp_p[s,t,])
        
      }
      
      x[s,t] <- rcat(1, p[s,t,])
      
      r[s,t] <- payoff[t, x[s,t]]
      
    }
    
  }
  
  result <- list(x=x,
                 r=r,
                 Q=Q,
                 p=p,
                 a=a,
                 a_mean=a_mean,
                 a_sd=a_sd,
                 beta=beta,
                 beta_mean=beta_mean,
                 beta_sd=beta_sd,
                 minor=minor,
                 minor_mean=minor_mean,
                 minor_sd=minor_sd,
                 inter=inter,
                 inter_mean=inter_mean,
                 inter_sd=inter_sd,
                 exp_p=exp_p,
                 p=p,
                 n_subj=n_subj,
                 ntrials=ntrials)
  return(result)
  
}

make_df_from_model <- function(X,RW_sims,i,i_time,seed){
  
  n_subj = RW_sims$n_subj
  ntrials = RW_sims$ntrials
  
  new_params <- c("inter_mean",
                  "inter_sd","minor_mean","minor_sd",
                  "a_mean", "a_sd", "beta_mean", "beta_sd")
  #make template for
  
  output_df <- tibble(
    seed = rep(seed,length(new_params)),
    run = rep(i,length(new_params)),
    i_time =rep(i_time,length(new_params)),
    param = new_params
  ) %>% 
    mutate(
      
      t_val = case_when(
                        param == "inter_mean" ~ RW_sims$inter_mean,
                        param == "minor_mean" ~ RW_sims$minor_mean,
                        param == "inter_sd" ~ RW_sims$inter_sd,
                        param == "minor_sd" ~ RW_sims$minor_sd,
                        
                       
                        param == "a_mean" ~ RW_sims$a_mean,
                        param == "a_sd" ~ RW_sims$a_sd,
                        param == "beta_mean" ~ RW_sims$beta_mean,
                        param == "beta_sd" ~ RW_sims$beta_mean
      ),
      # cheeat sheet
      # estimate distribution,index_participant,index_trial,index_arm
      # MPD(X$Q[,5,72,3])
      i_val = case_when(
                        param == "inter_mean" ~ MPD(X$inter_mean),
                        param == "minor_mean" ~ MPD(X$minor_mean),
                        param == "inter_sd" ~ MPD(X$inter_sd),
                        param == "minor_sd" ~ MPD(X$minor_sd),
                        
                        param == "a_mean" ~ MPD(X$a_mean),
                        param == "a_sd" ~ MPD(X$a_sd),
                        param == "beta_mean" ~ MPD(X$beta_mean),
                        param == "beta_sd" ~ MPD(X$beta_sd)
      )
      
    )
  
  return(output_df)
  
}

