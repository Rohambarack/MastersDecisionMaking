#setup
library(R2jags)
# building models

set.seed(123)

# model 1 --"fixed skill level"
# simulation, G is "Guess", t is "trial"
ntrials <- 100
Gfixed <- array(NA, c(ntrials))
theta <- 0.7

for (t in 1:ntrials){
  Gfixed[t] <- rbinom(1,1,theta)
}

plot(Gfixed)

# model 2 --"learning model"
ntrials <- 100
Glearn <- array(NA, c(ntrials))
#learning in time
theta_learn <- array(NA, c(ntrials))
theta_1 <- 0.5
theta_learn[1] <- theta_1

alpha <- 0.05
# set first guess
Glearn[1] <- rbinom(1,1,theta_1)

# 2:ntrials for a lagged time
for (t in 2:ntrials){
  #learning
  theta_learn[t] <- theta_learn[t-1]^(1/(1+alpha))
  #guessing
  Glearn[t] <- rbinom(1,1,theta_learn[t])
}


plot(theta_learn, ylim = c(0,1))
points(Glearn, col = 2)
title("M2")

# model 3 --"learning model_2, decaying alpha?"
ntrials <- 100
Glearn <- array(NA, c(ntrials))
#learning in time
theta_learn <- array(NA, c(ntrials))
theta_1 <- 0.5
theta_learn[1] <- theta_1

alpha_decay <- array(NA, c(ntrials))
alpha_1 <- 0.05
alpha_decay[1] <- alpha_1

d_rate <- 0.98
# set first guess
Glearn[1] <- rbinom(1,1,theta_1)

# 2:ntrials for a lagged time
for (t in 2:ntrials){
  #learning
  alpha_decay[t] <- alpha_decay[t-1]*d_rate
  theta_learn[t] <- theta_learn[t-1]^(1/(1+alpha_decay[t]))
  #guessing
  Glearn[t] <- rbinom(1,1,theta_learn[t])
}

plot(theta_learn, ylim = c(0,1))
points(Glearn, col = 2)
title("M3, alpha[t] <- alpha[t-1]*0.98")