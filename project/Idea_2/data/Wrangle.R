
library(stringr)
library(tidyverse)


################# Extract values from data                     
ntrials <- 72
practice <- 2
#count when TTB fails and adjust for 2 practice rounds
TTB_fail <- c(10,11,17,25,37,38,40,44,
              59,67,68,73) - practice
TTB <- rep(1,ntrials)
TTB[TTB_fail] <- 0
#count when not listening to TTB causes Treasure TTB fails + dublicate rewards
duplicate_reward <- c(3,8,14,18,
                      19,20,22,26,
                      28,29,31,33,
                      35,39,41,47,
                      53,54,55,56,
                      61,65,71,72) - practice
non_TTB_win <- c(TTB_fail,duplicate_reward)
non_TTB <- rep(0,ntrials)
non_TTB[non_TTB_win] <- 1
##### make payoff
payoff <- cbind(TTB,non_TTB)
##### find Cue3 in minority
minority <- seq(5,74,by = 3) - practice
m_indicator <- rep(0,ntrials)
m_indicator[minority] <- 1

#### get TTB actual choices
TTB_Choices <- c(1,1,1,
                 2,2,2,
                 1,2,2,
                 2,1,1,
                 1,2,1,
                 2,1,2,
                 1,2,2,
                 2,1,1,
                 #b1
                 2,1,2,
                 2,2,2,
                 2,1,1,
                 1,1,1,
                 1,2,1,
                 1,2,2,
                 2,2,2,
                 1,1,1,
                 #b2
                 2,1,2,
                 1,1,1,
                 1,1,2,
                 2,2,1,
                 1,2,1,
                 1,2,1,
                 2,2,2,
                 2,1,2)

#######################
#write_rds(payoff,"experiment_payoff.rds")
#write_rds(m_indicator,"TTB_minority.rds")
intervention_1 <- c(rep(0,24),rep(1,48))
intervention_2 <- c(rep(0,48),rep(1,24))
#write_rds(intervention_2,"inter_2.rds")

#wrangle the data
data <- read.csv("TeachingTTB_anonym.csv", sep = ";")
#clean out NAs
data <- drop_na(data)
#select relevant columns
c_data <- data %>% 
  select(ID,induction_time,MC1_smartestanimal,ends_with("Choice"))

c_data %>% 
  group_by(MC1_smartestanimal) %>% 
  reframe(n())
#99 did infer correctly, use them !
c_data <- c_data %>% 
  filter(MC1_smartestanimal == 3)

c_data %>% 
  group_by(induction_time) %>% 
  reframe(n())

c_data <- c_data %>% 
  select(!MC1_smartestanimal) %>% 
  select(!D1_Choice) %>% 
  select(!D2_Choice)

ID_list <- unique(c_data$ID)

test_1 <- c_data %>% 
  filter(induction_time == 1)

test_1 <- pivot_longer(c_data, cols = starts_with("D"), names_to = "turn")
test_1$turn <- str_replace_all(test_1$turn, "_Choice","")
test_1$turn <- str_replace_all(test_1$turn, "D","")
test_1$turn <- as.integer(test_1$turn)


#separate by ID
for (i in ID_list){
  
  test <- c_data %>% 
    filter(ID == i)
  #pivot long
  test <- pivot_longer(test, cols = starts_with("D"), names_to = "turn")
  test$turn <- str_replace_all(test$turn, "_Choice","")
  test$turn <- str_replace_all(test$turn, "D","")
  test$turn <- as.integer(test$turn)
  #compare to TTB
  test$compare <- as.integer(test$value == TTB_Choices)
  #rebrand as if 1 is TTB 2 is non_TTB 
  test$x <- ifelse(test$compare == 1,1,2)
  
  
  #recompile df
  
  if(i == ID_list[1]){ 
    
    data_df <- test
    
  }else{
    temp_df <- test
    data_df <- rbind(data_df,test)
    
  }
  
}

write_rds(data_df,"c_data.rds")


# 2414 ID <- 2/3 intervention
# 1710 ID <- TTB by default
# 2515 ID <- 1/3 intervention
# 901 ID <- Not really bothered

test_2 <- test_1 %>% 
  filter(ID == 901)
unique(test_1$ID)
test_2 %>% 
  ggplot(aes(y=value, x = turn)) +
  geom_point()
#compare to TTB

test_2 %>% 
  ggplot(aes(y=compare, x = turn)) +
  geom_point()

#payoff <- readRDS("data/experiment_payoff.rds")
payoff

########## wrangle out rewards based on payoff and x
n_subj <- 99
ntrials <- 72
#cut up the df to make a (nsubj*ntrials) matrix of x
x_m <- vector("list", n_subj)
for(i in 1:n_subj){
  #separate 1 subject at a time
  x_1_subj <- data_df$x[(((i-1)*ntrials)+1):(i*ntrials)]
  #make it a list of vectors
  x_m[[i]] <- x_1_subj
  
}
x_matrix <- do.call(rbind, x_m)

x_matrix[1,] == data_df %>% filter(ID==1603) %>% select(x)

#find r
for (s in 1:n_subj){
  for (t in 1:ntrials){

r[s,t] <- payoff[t, x[s,t]]

}
}