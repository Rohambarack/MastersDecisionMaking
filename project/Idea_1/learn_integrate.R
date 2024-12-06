
### to do
#biases for animals
#update other slots?
# Model 2, 
# intervention changes the probabilities of choosing, but choice also ends
# up modifying Q_exp
# 


RW <- function(payoff,
               ntrials,
               alpha,
               beta,
               w_s,
               interv,
               interv_values,
               w_interv
) {  
  
  
  #alpha <- 0.001
  #beta <- 0.5
  
  x <- array(0, c(ntrials))
  r <- array(0, c(ntrials))
  Q <- array(0,c(ntrials,3))
  Qupdate <- array(0, c(ntrials, 3))
  exp_p <- array(0,c(ntrials,3))
  exp_i <- array(0,c(ntrials,3))
  i <- array(0, c(ntrials,3))
  p <- array(0, c(ntrials,3))
  
  #set trial to 1 according to setup, estimate from known hits and misses
  #estimating how good they are at guessing true Q
  Q[1,1] <- rbeta(1,3*w_s,3*w_s)
  Q[1,2] <- rbeta(1,4*w_s,2*w_s)
  Q[1,3] <- rbeta(1,5*w_s,1*w_s)
  
  for (t in 2:ntrials) {
    
    for(k in 1:3) {
      
      Qupdate[t,k] <- Q[t-1,k] + (alpha*(r[t-1]- Q[t-1,k]))
      #only update the one chosen
      Q[t,k] <- ifelse(k==x[t-1],
                       Qupdate[t,k],
                       Q[t-1,k])
      #softmax
      exp_p[t,k] <- exp(beta*Q[t,k])
      
    }
    
    for (k in 1:3) {
      
      p[t,k] <- exp_p[t,k]/sum(exp_p[t,])
      
    }
    
    #before intervention
    if(interv[t,1] == 0){
      x[t] <- rcat(1,p[t,])
    } else {
      #prob? between own and intervened
      inter_choic <- rbinom(1,1,w_interv)
      
      if(inter_choic == 1){
        x[t] <- rcat(1,interv_values[t,])
      } else{
        x[t] <- rcat(1,p[t,])
      }
    }
    
    r[t] <- payoff[t, x[t]]
    
  }
  
  result <- list(x=x, r=r, Q=Q)
  
  return(result)
}
