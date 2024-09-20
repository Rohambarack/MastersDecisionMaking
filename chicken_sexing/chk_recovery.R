# Building models
#install.packages("R2jags")
library(R2jags)

set.seed(1983)

#function to iterate over data
MPD <- function(x){
  
  density(x)$x[which(density(x)$y == max(density(x)$y))]
  
  
}

#-------- Model 1 - fixed skill level (theta) ------
#simulation - running an experiment as if we knew the process and the parameter
# "G" is for "Guess" 

ntrials <- 100
Gfixed <- array(NA, c(ntrials))
theta <- 0.7

for (t in 1:ntrials) {
  Gfixed[t] <- rbinom(1,1,theta)
}

plot(Gfixed)


#----- Condition model on data (JAGS) -----


data <- list("Gfixed", "ntrials")
params <- c("theta")

fixed_samples <- jags(data,
                      inits = NULL,
                      params,
                      model.file = "chick_jags.txt",
                      n.chains = 3,
                      n.iter = 5000,
                      n.burnin = 1000,
                      n.thin = 1)




#-------- Model 2 - learning model ------
#simulation - running an experiment as if we knew the process and the parameter

ntrials <- 100
Glearn <- array(NA, c(ntrials))

theta_learn <- array(NA, c(ntrials))
alpha <- 0.05
theta1 <- 0.5

theta_learn[1] <- theta1
Glearn[1] <- rbinom(1,1,theta1)

for (t in 2:ntrials) {
  theta_learn[t] <- theta_learn[t-1]^(1/(1+alpha))
  Glearn[t] <- rbinom(1,1,theta_learn[t])
}

plot(theta_learn)
plot(Glearn, col=2) # plotting in red

# mimic "hold on" - i.e. overlay two plots
plot(Gfixed)
points(Glearn, col=2) # overlaying in red


#------------model it

data <- list("Glearn", "ntrials")
params <- c("theta1","theta","alpha")

learn_samples <- jags(data,
                      inits = NULL,
                      params,
                      model.file = "chick_l_jags.txt",
                      n.chains = 3,
                      n.iter = 5000,
                      n.burnin = 1000,
                      n.thin = 1)

plot(density(learn_samples$BUGSoutput$sims.list$alpha))
plot(density(learn_samples$BUGSoutput$sims.list$theta1))


X <- learn_samples$BUGSoutput$sims.list
MPD(X$alpha)


