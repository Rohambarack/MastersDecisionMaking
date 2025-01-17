##############read in data
library(tidyverse)
library(stringr)
list_of_files <- list.files("../results")
path <- "../results"

counting <- 1
for (i in list_of_files){
  
  if (counting == 1){
    temporary_path <- paste(path,i,sep="/")
    df <- readRDS(temporary_path)
  } else {
    temporary_path <- paste(path,i,sep="/")
    temporary_df <- readRDS(temporary_path)
    df <- rbind(df,temporary_df)
  }
  
  counting <- counting +1
  
}

#############
df %>% 
  #filter(str_detect(param,"inter|minor")) %>% 
  ggplot(aes(x = t_val, y =i_val, colour = as.factor(i_time))) +
  geom_point(alpha = .4) +
  #lm line
  geom_smooth(method = "lm")+
  #ref line
  geom_abline(color="red",intercept = 0, slope = 1 ) +
  scale_colour_manual(values = c("#cd7e3c","#3C8BCD"),
                      name = "Group") +
  xlab("True Value") +
  ylab("Inferred Value")+
  theme_classic() +
  facet_wrap(~param, scales = "free")
