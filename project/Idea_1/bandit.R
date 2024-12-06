# setting up a bandit-model
install.packages("pacman")
pacman::p_load(hesim, extraDistr, R2jags)
#pacman::p_load(hesim, extraDistr, R2jags, parallel, ggpubr)

set.seed(1983)

# payoff structure
# generate a payoff matrix for the bandit task
# Choice of bandit A = 30% chance of 2. Choice of bandit B = 70% chance of 1. Otherwise nothing.
ntrials <- 74
C1prob <- 0.5
C2prob <- 0.63
C3prob <- 0.83
rew <- 1


# there's all sorts of ways to generate a payoff matrix. We can use conditions or whatever.
# but let's use binomials because it's really efficient
# why do we mulitply the output of the binomial draw?
payoff <- cbind(rbinom(ntrials,1,C1prob)*rew,
                rbinom(ntrials,1,C2prob)*rew,
                rbinom(ntrials,1,C3prob)*rew)

inter_matrix <- cbind(rep(c(0,1), each = ntrials/2),
                      rep(c(0,1), each = ntrials/2),
                      rep(c(0,1), each = ntrials/2))

inter_matrix_values <- cbind(rep(0, ntrials),
                            rep(0, ntrials),
                            rep(1, ntrials))


# let's look at which option is best in the long run
colSums(payoff)

#theta <- .7 # bias for one bandit over the other
#b <- c(theta, 1-theta) # stated in terms of our categorical choice probabilities

# Don't forget to set your working directory
#setwd("Where/Are/My/Files")

source("learn_sep.R") # NB! working directory
random_sims <- RW(payoff = payoff,
                  ntrials = 74,
                  alpha = 0.05,
                  beta = 0.5,
                  w_s = 10,
                  interv = inter_matrix,
                  interv_values = inter_matrix_values,
                  w_interv = 0.7,
                  bias_vect = c(0.5,0.1,-0.2))

#vistest
curve(dbeta(x,3,3),from = 0, to = 1.5)
curve(dbeta(x,4,2),from = 0, to = 1, col = "red", add = T)
curve(dbeta(x,5,1),from = 0, to = 1, col = "blue", add = T)


curve(dbeta(x,3000,3000), from = 0, to = 1)
curve(dbeta(x,3,3),from = 0, to = 1, col = "red", add = T)

### test softmax
temp <- .1
values <- c(0,0,1)
exp_values <- exp(temp*values)
sum_v <- sum(exp_values)
soft_values <- exp_values/sum_v

curve(dbinom(x,1,0.7), from =0, to = 1)
