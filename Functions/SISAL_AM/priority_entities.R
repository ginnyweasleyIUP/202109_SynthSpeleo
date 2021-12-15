library(rJava)
library(xlsx)
#library(dplyr)
library(tidyverse)
library(janitor)

wd <- '/home/ariana/SISAL Data/sisalv1b_csv'
setwd(wd)
setwd("~/SISAL Data/interim_v2_csv_20190812")

dating <- read.csv('dating.csv', header = T, stringsAsFactors = F)
original_chronology <- read.csv('original_chronology.csv', header = T, stringsAsFactors = F)
hiatus <- read.csv('hiatus.csv', header = T, stringsAsFactors = F)
sample <- read.csv('sample.csv', header = T, stringsAsFactors = F)
entity <- read.csv('entity.csv', header = T, stringsAsFactors = F)

dating_tb <- dating %>% left_join(.,entity, by = "entity_id") %>% filter(entity_status == 'current') %>%
  mutate_at(vars(dating_id, depth_dating, dating_thickness, X14C_correction, corr_age, corr_age_uncert_pos, corr_age_uncert_neg), as.numeric)
original_tb <- sample %>% left_join(., original_chronology, by = 'sample_id') %>% left_join(., entity, by = "entity_id") %>% filter(entity_status == 'current')
sample_tb <- sample %>% left_join(., hiatus, by = "sample_id") %>% left_join(., entity, by = "entity_id") %>% filter(entity_status == 'current') 
hiatus <- sample_tb %>% filter(entity_id %in% dating_tb$entity_id & hiatus == 'H' )

mis1 <- c(-70,10000) 
mis2 <- c(19000,29000)
mis3 <- c(29000,57000)
mis4 <- c(57000,71000)
mis5 <- c(71000,130000)



# quantify MIS stages
mis_periods <- dating_tb %>% select(entity_id, depth_dating, corr_age, corr_age_uncert_pos, corr_age_uncert_neg, speleothem_type) %>% group_by(entity_id) %>%
  summarise(min_depth_dating = round(min(depth_dating, na.rm = T),digits = 2),
            max_depth_dating = round(max(depth_dating, na.rm = T),digits= 2),
            depth = max_depth_dating - min_depth_dating,
            min_corr_age = round(min(corr_age, na.rm = T),digits = 2),
            max_corr_age = round(max(corr_age, na.rm = T),digits = 2),
            age = max_corr_age-min_corr_age,
            corr_age_uncert = if_else(any(!is.na(corr_age_uncert_neg) | !is.na(corr_age_uncert_pos)), TRUE, FALSE),
            number_dates = length(corr_age),
            #number_hiatus = if_else(entity_id %in% hiatus$entity_id, sum(entity_id %in% hiatus$entity_id), integer(NA)),
            MIS1 = if_else(any(corr_age > mis1[1] & corr_age < mis1[2], na.rm = T),TRUE, FALSE),
            MIS2 = if_else(any(corr_age > mis2[1] & corr_age < mis2[2], na.rm = T),TRUE, FALSE),
            MIS3 = if_else(any(corr_age > mis3[1] & corr_age < mis3[2], na.rm = T),TRUE, FALSE),
            MIS4 = if_else(any(corr_age > mis4[1] & corr_age < mis4[2], na.rm = T),TRUE, FALSE),
            MIS5 = if_else(any(corr_age > mis5[1] & corr_age < mis5[2], na.rm = T),TRUE, FALSE),
            MIS_plus = if_else(any(corr_age > mis5[2], na.rm = T),TRUE, FALSE),
            speleothem_type = speleothem_type[1])

int <- mis_periods 
int_new <- int %>% filter(min_corr_age <= 0 & max_corr_age >= 1100 | entity_id %in% c(237, 305, 312, 320)) 
int_depths <- dating_tb %>% filter(entity_id %in% int_new$entity_id | entity_id %in% c(237, 305, 312, 320)) %>% group_by(entity_id) %>%
  filter(corr_age < 1100) %>% filter(date_used == 'yes') %>%select(entity_id, corr_age, depth_dating, date_type)
i <- int_depths %>% filter(!(date_type %in% c('C14', 'Event; end of laminations', 'Event; start of laminations')))
#i_new <- i %>% group_by(entity_id) %>% summarize(max_dft = max(depth_dating))
#sample_nr <- left_join(sample_tb,i_new, by='entity_id') %>%filter(entity_id %in% i_new$entity_id)# %>% group_by(entity_id) %>% 
#  filter(depth_sample < max_dft) %>% count() %>% filter(n>50)

s_nr <- original_tb %>% filter(entity_id %in% i$entity_id) %>% filter(interp_age < 1100) %>% group_by(entity_id) %>% count()
  
test <- dating_tb %>% group_by(entity_id) %>% arrange(depth_dating, .by_group = T) %>% filter(entity_id %in% int$entity_id) %>% filter(corr_age <= 10000) %>% 
  group_by(entity_id) %>% count() %>% filter(n >= 3)
test2 <- dating_tb %>% filter(entity_id %in% test$entity_id) %>% filter(corr_age <= 29000 & corr_age >= 19000) %>% group_by(entity_id) %>% count() %>% filter(n >= 3)
very_int <- int %>% filter(entity_id %in% class2$entity_id)
  