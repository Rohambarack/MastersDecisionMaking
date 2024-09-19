#slightly offset step function
offset_step <- function(x,step_change){
  x_2 <- x-step_change
  val <- ((x_2+abs(x_2))/(2*x_2))
  
  return(val)
}

curve(offset_step(x,0.7),from = -1, to = 3)

# model 4 --"alpha decay, transcendental intervention paradigm-shift"
ntrials <- 120
Glearn <- array(NA, c(ntrials))
fails <- array(NA, c(ntrials))
#learning in time
theta_learn <- array(NA, c(ntrials))
theta_1 <- 0.5
theta_learn[1] <- theta_1
#alpha decays with each guess
alpha_decay <- array(NA, c(ntrials))
alpha_1 <- 0.05
alpha_decay[1] <- alpha_1
d_rate <- 0.94
#newcastle chicken sexer apparition chance
newcastle_chance <- array(NA, c(ntrials))
nci <- array(NA, c(ntrials))
n_limit <- 0.7
n_rate <- 10
newcastle_transcend <- array(NA, c(ntrials))
newcastle_count <- array(NA, c(ntrials))
newcastle_lock <- array(NA, c(ntrials))

# set first guess
Glearn[1] <- rbinom(1,1,theta_1)
fails[1] <- abs(Glearn[1]-1)
# determine first nci, newcastle_chance, newcastle_transcend
nci[1] <- offset_step(theta_learn[1],n_limit) * fails[1]
newcastle_chance[1] <- nci[1]/(nci[1] + n_rate)
newcastle_transcend[1] <- rbinom(1,1,newcastle_chance[1])
newcastle_count[1] <- cumsum(newcastle_transcend)[1]
newcastle_lock[1] <- offset_step(newcastle_count[1],0.9)


# 2:ntrials for a lagged time
for (t in 2:ntrials){
  #learning below 0.7 theta
  #alpha decays, but if newcastle_lock is on, no decay, 
  #      + alpha becomes 0.5 again
  alpha_decay[t] <- alpha_decay[t-1]*((1-newcastle_lock[t-1])*d_rate) + (abs(0-newcastle_lock[t-1])*alpha_decay[1])
  
  #diff(alpha_decay is only + when it resets to 0.5)
  #so offset_step(diff(alpha_decay)[t-1],0.01) is only 1 when that happens so 
  #(1-offset_step(diff(alpha_decay)[t-1],0.01)) multiplication is 1, except at change
  #(abs(0-offset_step(diff(alpha_decay)[t-1],0.01))*0.5) that part of the equation is always 0, except at change whtn it is 0.5
  theta_learn[t] <- (theta_learn[t-1]^(1/(1+alpha_decay[t]))) * (1-offset_step(diff(alpha_decay)[t-1],0.01)) + (abs(0-offset_step(diff(alpha_decay)[t-1],0.01))*0.5)
  Glearn[t] <- rbinom(1,1,theta_learn[t])
  fails[t] <- abs(Glearn[t]-1)
  #newcastle chicke sexer prob, only above 0.7 theta
  #offset_step function is 0 below 0.7, 1 above 0.7
  #multiplied by Guess to find successes, inverted to see fails
  #so nci is 1 when there is a fail above 0.7 theta
  #nci[t-1] is added to make it cumulative
  nci[t] <- (offset_step(theta_learn[t],n_limit) * fails[t]) + nci[t-1]
  #newcastle_chance is is always between 0 and 1 based on (x/x+a)
  newcastle_chance[t] <- nci[t]/(nci[t] + n_rate)
  newcastle_transcend[t] <- rbinom(1,1,newcastle_chance[t])
  #lock keeps it 1 after 1st occurence
  #counts occurences, step fuction makes it 1 if bigger than 0.9, so
  #everytime above 0 since it is only integers
  newcastle_count[t] <- cumsum(newcastle_transcend)[t]
  newcastle_lock[t] <- offset_step(newcastle_count[t],0.9)
  
}


plot(theta_learn, ylim = c(0,1))
points(Glearn, col = 2)
#plots a line where alpha decay becomes 0.05 again
abline(v = which(alpha_decay == 0.05)[2], lty = 3)
title("T.I.P.S")
