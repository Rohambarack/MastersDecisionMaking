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

#modified rescorla Wagner with group level effects for 2 bandits
# and group level intervention effect, and truth in minority
# only on Q for th TTB arm
RW_hier_inter_minor_2 <- function(payoff, intervention, minority,
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
                                   Q[s,1,k] + (intervention[1]*inter[s]) + (minor[1]*minor[s]),
                                   Q[s,1,k])
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
                                     Q[s,t,k] + (intervention[t]*inter[s]) + (minority[t]*minor[s]),
                                     Q[s,t,k])
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

retrieve_value <- function(parameter,n_subj,n_trial,n_arm){
  #all participant, works for 2 arms
  test_4 <- parameter
  Q_list <- array(NA,c(n_subj,n_trial,n_arm))
  #loop trhrough the estimates by participant, arm, and turn
  #make a list of mean estimates my finding the maximum probability density
  for (p in 1:n_subj){
    for (i in 1:n_arm){
      for (k in 1:n_trial){
        i_arm_k_turn <- MPD(test_4[,p,k,i])
        Q_list[p,k,i] <- i_arm_k_turn
        
      }
    }
  }
  
  #wrangle estimates into a ggplot acceptable format
  for (p in 1:n_subj){
    if(p ==1){
      all_Q <- tibble(mu = rep("individual",n_trial*n_arm),
                      participant = rep(p,,n_trial*n_arm),
                      arm = rep(c("TTB","non_TTB"), each = n_trial),
                      turn = c(1:n_trial,1:n_trial),
                      p_val = c(Q_list[p,,1],
                                Q_list[p,,2]))
    } else {
      temp_Q <- tibble(
        mu = rep("individual",n_trial*n_arm),
        participant = rep(p,n_trial*n_arm),
        arm = rep(c("TTB","non_TTB"), each = n_trial),
        turn = c(1:n_trial,1:n_trial),
        p_val = c(Q_list[p,,1],
                  Q_list[p,,2]))
      all_Q <- rbind(all_Q, temp_Q)
    }
  }
  
  #also find average for all turns
  mean_Q_list <- array(NA,c(n_trial,n_arm))
  for (t in 1:n_trial){
    for (a in 1:n_arm){
      mean_Q_list[t,a] <- mean(Q_list[,t,a])
    }
  }
  
  mean_Q_df <- tibble(mu = rep("mean",n_trial*n_arm),
                      participant = rep("mean",n_trial*n_arm),
                      arm = rep(c("TTB","non_TTB"), each = n_trial),
                      turn = c(1:n_trial,1:n_trial),
                      p_val = c(mean_Q_list[,1],
                                mean_Q_list[,2]))
  
  all_Q <- rbind(all_Q,mean_Q_df)
  return(all_Q)
}

retrieve_value_sim <- function(parameter,n_subj,n_trial,n_arm){
  #all participant, works for 2 arms
  Q_list <- parameter

  #wrangle estimates into a ggplot acceptable format
  for (p in 1:n_subj){
    if(p ==1){
      all_Q <- tibble(mu = rep("individual",n_trial*n_arm),
                      participant = rep(p,,n_trial*n_arm),
                      arm = rep(c("TTB","non_TTB"), each = n_trial),
                      turn = c(1:n_trial,1:n_trial),
                      p_val = c(Q_list[p,,1],
                                Q_list[p,,2]))
    } else {
      temp_Q <- tibble(
        mu = rep("individual",n_trial*n_arm),
        participant = rep(p,n_trial*n_arm),
        arm = rep(c("TTB","non_TTB"), each = n_trial),
        turn = c(1:n_trial,1:n_trial),
        p_val = c(Q_list[p,,1],
                  Q_list[p,,2]))
      all_Q <- rbind(all_Q, temp_Q)
    }
  }
  
  #also find average for all turns
  mean_Q_list <- array(NA,c(n_trial,n_arm))
  for (t in 1:n_trial){
    for (a in 1:n_arm){
      mean_Q_list[t,a] <- mean(Q_list[,t,a])
    }
  }
  
  mean_Q_df <- tibble(mu = rep("mean",n_trial*n_arm),
                      participant = rep("mean",n_trial*n_arm),
                      arm = rep(c("TTB","non_TTB"), each = n_trial),
                      turn = c(1:n_trial,1:n_trial),
                      p_val = c(mean_Q_list[,1],
                                mean_Q_list[,2]))
  
  all_Q <- rbind(all_Q,mean_Q_df)
  return(all_Q)
}

homebrew_pp_check <- function(recovered_preds,data){
  probs <- recovered_preds
  df_data <- data
  #remove mean, pivot wider
  probs_preds <- probs %>% 
    filter(mu == "individual") %>% 
    pivot_wider(names_from = arm, values_from = p_val)
  
  #remove 0s
  probs_preds <- probs_preds %>% 
    mutate( TTB = ifelse(TTB < 0,0,TTB),
            non_TTB = ifelse(non_TTB < 0,0,non_TTB))
  
  seed_v <- 1:100
  for (s in seed_v){
    
    set.seed(s)
    
    probs_preds <- probs_preds %>% 
      rowwise() %>% 
      mutate( pred_x = rcat(1,c(TTB,non_TTB)))
    
    
    if (s == 1){
      
      pred_x_df <- tibble(
        pred_x = probs_preds$pred_x,
        seed = s,
        type = "draw",
      )
      
    }else{
      
      temp_x_df <- tibble(
        pred_x = probs_preds$pred_x,
        seed = s,
        type = "draw",
      )
      
      pred_x_df <- rbind(pred_x_df,temp_x_df)
      
    }
    
  }
  
  #add true vals
  true_x_df <- tibble(
    pred_x = df_data$x,
    seed = 101,
    type = "real",
  )
  
  pred_x_df <- rbind(pred_x_df,true_x_df)
  
  plot_x <- pred_x_df %>% 
    mutate(pred_x = ifelse(pred_x == 1,1,0)) %>% 
    ggplot(aes(x = pred_x, group = seed, colour = type))+
    scale_color_manual(values = c("lightblue","black"),
                       name = "",
                       labels = c("Y_rep","Y")) +
    scale_alpha_manual(values = c(0.05, 1)) +
    geom_density() +
    ylab("")+
    xlab("") +
    theme_classic() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  return(plot_x)
  
}

simulate_glm_data_non_diff <- function(intercept = logit(0.6),
                                       sigma_ = 0.15,
                                       b_inter = 1,
                                       sigma_inter = 0.33,
                                       b_minor = -1,
                                       sigma_minor= 0.33){
  #mu and sigma values chosen with logit scale in mind!!!
  b_intercept_X1 <- intercept
  #above chance TTB by default 60%
  sigma_X1 <- sigma_ 
  # almost all probabilty mass between 3 sd, 3*0.15=0.45
  # invlogit(logit(0.6)+0.45) = 0.48
  # invlogit(logit(0.6)-0.45) = 0.70
  # values between 48% and 70% for choosing TTB by default seem reasonable
  # same logic applies for the rest..
  # b_inter and b_minor mu
  # invlogit(logit(0.6)+1) = 0.8
  # invlogit(logit(0.6)-1) = 0.35
  # sigma 
  # invlogit(logit(0.6)+2) = 91, so 3*sd should reach 1
  b_inter_X1 <- b_inter
  sigma_inter_X1 <- sigma_inter
  b_minor_X1 <- b_minor
  sigma_minor_X1 <- sigma_minor
  # same for X_2, no difference assumed
  b_intercept_X2 <-  b_intercept_X1
  sigma_X2 <- sigma_X1
  b_inter_X2 <- b_inter_X1
  sigma_inter_X2 <- sigma_inter_X1
  b_minor_X2 <- b_minor_X1
  sigma_minor_X2 <- sigma_minor_X1 
  
  #make intervention templates, minority template
  inter_1_temp <- c(rep(0,24),rep(1,48))
  inter_2_temp <- c(rep(0,48),rep(1,24))
  #every 3rd is a minority case, 2 practice trials
  minority <- seq(5,74,by = 3) - 2
  minor_temp <- rep(0,72)
  minor_temp[minority] <- 1
  
  #make base df, ID, induction_time
  data <- tibble(ID = 1:100,
                 induction_time = rep(c(1,2),each = 50))
  #add betas and sigmas by group
  data <- data %>% 
    mutate(intercept_mu = ifelse(induction_time == 1,
                                 b_intercept_X1,
                                 b_intercept_X2),
           intercept_sigma = ifelse(induction_time == 1,
                                    sigma_X1,
                                    sigma_X2),
           b_inter = ifelse(induction_time == 1,
                            b_inter_X1,
                            b_inter_X2),
           b_minor = ifelse(induction_time == 1,
                            b_minor_X1,
                            b_minor_X2),
           sigma_inter = ifelse(induction_time == 1,
                                sigma_inter_X1,
                                sigma_inter_X2),
           sigma_minor = ifelse(induction_time == 1,
                                sigma_minor_X1,
                                sigma_minor_X2),
           
    )
  
  #calculate individual values by ID
  
  data <- data %>% 
    rowwise() %>% 
    mutate(
      ind_intercept = rnorm(1,intercept_mu,intercept_sigma),
      ind_b_inter = rnorm(1,b_inter,sigma_inter),
      ind_b_minor = rnorm(1,b_minor,sigma_minor),
    ) %>% 
    ungroup()
  
  #expand on individual params, drop global (not necessary to drop them)
  
  data <- data %>% 
    expand(nesting(ID,
                   induction_time,
                   ind_intercept,
                   ind_b_inter,
                   ind_b_minor),turn =1:72)
  
  #add predictors intervention and minority
  data <-data %>% 
    mutate(minor = rep(minor_temp,100),
           inter = c(rep(inter_1_temp,50),rep(inter_2_temp,50))
    )
  
  # calculate y
  data <- data %>% 
    rowwise() %>% 
    mutate(
      logit_y = ind_intercept + ind_b_inter*inter + ind_b_minor*minor,
      prob_y = invlogit(logit_y),
      y =rbern(1,prob_y)
      
      
    ) %>% 
    ungroup()
  #factor factors
  data <- data %>% 
    mutate(
      induction_time = as.factor(induction_time),
      inter= as.factor(inter),
      minor= as.factor(minor)
    )
  return(data)
}

make_df_from_model_glm <- function(mod,sim_params,i,seed){
  
  mod_df <- as_draws_df(mod)
  
  new_params <- c("intercept_mean", "intercept_sd",
                  "inter_sd","inter_mean",
                  "minor_mean","minor_sd")
  seed = 1
  i = 1
  
  output_df <- tibble(
    seed = rep(seed,length(new_params)*2),
    run = rep(i,length(new_params)*2),
    i_time =rep(c(1,2), each = length(new_params)),
    param = rep(new_params,2)
  ) %>% 
    mutate(
      
      t_val = case_when(
        param == "intercept_mean" ~ sim_params[1],
        param == "intercept_sd" ~ sim_params[2],
        param == "inter_mean" ~ sim_params[3],
        param == "inter_sd" ~ sim_params[4],
        param == "minor_mean" ~ sim_params[5],
        param == "minor_sd" ~ sim_params[6]
      ),
      
      i_val = case_when(
        param == "intercept_mean" & i_time == 1 ~ MPD(mod_df$b_induction_time1),
        param == "intercept_sd" & i_time == 1 ~ MPD(mod_df$`sd_ID__Intercept:induction_time1`),
        param == "inter_mean" & i_time == 1 ~ MPD(mod_df$b_inter1),
        param == "inter_sd" & i_time == 1 ~ MPD(mod_df$`sd_ID__inter1:induction_time1`),
        param == "minor_mean" & i_time == 1 ~ MPD(mod_df$b_minor1),
        param == "minor_sd" & i_time == 1 ~ MPD(mod_df$`sd_ID__minor1:induction_time1`),
        
        param == "intercept_mean" & i_time == 2 ~ MPD(mod_df$b_induction_time2),
        param == "intercept_sd" & i_time == 2 ~ MPD(mod_df$`sd_ID__Intercept:induction_time2`),
        param == "inter_mean" & i_time == 2 ~ MPD(mod_df$b_inter1) + MPD(mod_df$`b_induction_time2:inter1`),
        param == "inter_sd" & i_time == 2 ~ MPD(mod_df$`sd_ID__inter1:induction_time2`),
        param == "minor_mean" & i_time == 2 ~ MPD(mod_df$b_minor1) + MPD(mod_df$`b_induction_time2:minor1`),
        param == "minor_sd" & i_time == 2 ~ MPD(mod_df$`sd_ID__minor1:induction_time2`)
        
      )
      
    )
  
  return(output_df)
  
}


RW_prediction <- function(recovered_preds,data){
  probs <- recovered_preds
  df_data <- data
  #remove mean, pivot wider
  probs_preds <- probs %>% 
    filter(mu == "individual") %>% 
    pivot_wider(names_from = arm, values_from = p_val)
  
  #remove 0s
  probs_preds <- probs_preds %>% 
    mutate( TTB = ifelse(TTB < 0,0,TTB),
            non_TTB = ifelse(non_TTB < 0,0,non_TTB))
  
  seed_v <- 1:100
  for (s in seed_v){
    
    set.seed(s)
    
    probs_preds <- probs_preds %>% 
      rowwise() %>% 
      mutate( pred_x = rcat(1,c(TTB,non_TTB)))
    
    
    if (s == 1){
      
      pred_x_df <- tibble(
        pred_x = probs_preds$pred_x,
        seed = s,
        type = "draw",
      )
      
    }else{
      
      temp_x_df <- tibble(
        pred_x = probs_preds$pred_x,
        seed = s,
        type = "draw",
      )
      
      pred_x_df <- rbind(pred_x_df,temp_x_df)
      
    }
    
  }
  
  pred_x_df$true_val <- rep(df_data$x,100)
  
  pred_x_df <- pred_x_df %>% 
    mutate(hits = ifelse(pred_x == true_val,1,0)) %>% 
    group_by(seed) %>% 
    reframe( acc = sum(hits)/nrow(df_data))
  
  
  
  
  return(pred_x_df)
}