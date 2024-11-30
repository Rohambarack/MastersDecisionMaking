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
###############wrangle it to format
w_df <- pivot_longer(df,cols = colnames(df)[3:length(colnames(df))],
                    names_to = "label",values_to = "values")

w_df_t <- w_df %>% 
  filter(str_detect(label, 'true')) %>% 
  rename(t_values = values) %>% 
  rowwise() %>% 
  mutate_at("label", str_replace,"true_", "")
w_df_i <- w_df %>% 
  filter(str_detect(label, 'infer')) %>% 
  rename(i_values = values) %>% 
  rowwise() %>% 
  mutate_at("label", str_replace,"infer_", "")

w_df_t$i_values <- w_df_i$i_values


############## do plots

w_df_t %>% 
  ggplot(aes(x = t_values, y =i_values)) +
  geom_point() +
  #lm line
  geom_smooth(method = "lm")+
  #ref line
  geom_abline(color="red",intercept = 0, slope = 1 ) +
  xlab("True Value") +
  ylab("Inferred Value")+
  theme_classic() +
  facet_wrap(~label,scales = "free")





