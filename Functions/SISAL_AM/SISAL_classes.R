library(rJava)
library(xlsx)
#library(dplyr)
library(tidyverse)
library(janitor)

#wd <- '/home/ariana/SISAL Data/sisalv1b_csv'
#wd <- "~/SISAL Data/interim_v2_csv_20190812"
wd <- "~/SISAL Data/sisalv2"
setwd(wd)

dating <- read.csv('dating.csv', header = T, stringsAsFactors = F) %>%
  mutate_at(vars(dating_id, depth_dating, dating_thickness, X14C_correction, corr_age, corr_age_uncert_pos, corr_age_uncert_neg), as.numeric)
hiatus <- read.csv('hiatus.csv', header = T, stringsAsFactors = F)
sample <- read.csv('sample.csv', header = T, stringsAsFactors = F)
entity <- read.csv('entity.csv', header = T, stringsAsFactors = F)
sisal <- read.csv('sisal_chronology.csv', header = T, stringsAsFactors = F)
orig <- read.csv('original_chronology.csv', header = T, stringsAsFactors = F)

dating_tb <- dating %>% left_join(.,entity, by = "entity_id") %>% filter(entity_status == 'current' | entity_status == 'current partially modified') %>%
  mutate_at(vars(dating_id, depth_dating, dating_thickness, X14C_correction, corr_age, corr_age_uncert_pos, corr_age_uncert_neg), as.numeric)
sample_tb <- sample %>% left_join(., hiatus, by = "sample_id") %>% left_join(., entity, by = "entity_id") %>% filter(entity_status == 'current' | entity_status == 'current partially modified') %>%
  mutate_at(vars(entity_id, sample_id, sample_thickness, depth_sample), as.numeric)

#dating_tb <- dating %>% left_join(.,entity, by = "entity_id") %>% filter(entity_status == 'current partially modified') %>%
#  mutate_at(vars(dating_id, depth_dating, dating_thickness, X14C_correction, corr_age, corr_age_uncert_pos, corr_age_uncert_neg), as.numeric)
#sample_tb <- sample %>% left_join(., hiatus, by = "sample_id") %>% left_join(., entity, by = "entity_id") %>% filter(entity_status == 'current partially modified') %>%
#  mutate_at(vars(entity_id, sample_id, sample_thickness, depth_sample), as.numeric)

no_pub_AM <- full_join(sample, sisal) %>% mutate_at(vars(COPRA_age), as.numeric) %>% filter(!is.na(COPRA_age)) %>% distinct(entity_id)

hiatus_tb <- sample_tb %>% filter(hiatus == 'H') %>% distinct(entity_id)
reversals <- dating_tb %>% filter(date_used == 'yes') %>%
  dplyr::select(entity_id, depth_dating, corr_age, corr_age_uncert_pos, corr_age_uncert_neg) %>%
  group_by(entity_id) %>%
  arrange(depth_dating, .by_group = T) %>%
  mutate(diff = lead(corr_age)-corr_age) %>%
  mutate(reversal = if_else(!is.na(diff) & diff < 0, TRUE, FALSE)) %>%
  group_by(entity_id) %>%
  arrange(depth_dating, .by_group = TRUE) %>%
  mutate(tractable = if_else(reversal & (abs(corr_age-lead(corr_age)) < (corr_age_uncert_pos + lead(corr_age_uncert_neg))), TRUE, FALSE))

reversal_tractable <- reversals %>% dplyr::select(entity_id, tractable) %>% group_by(entity_id) %>% filter(tractable) %>% distinct(entity_id)

reversal_nontractable <- reversals %>% filter(reversal) %>%filter(!(entity_id %in% reversal_tractable$entity_id)) %>% distinct(entity_id)
nr_dates <- sisalv2$dating %>% filter(date_used == 'yes' & date_type != 'Event; hiatus') %>% dplyr::select(entity_id, corr_age) %>% group_by(entity_id) %>%count() %>% filter(n>=3)
no_depth_sample <- sisalv2$sample %>% group_by(entity_id) %>% summarise(depths = if_else(all(is.na(depth_sample)), FALSE, TRUE)) %>% filter(!depths)

#UTh_dates <- dating_tb %>% filter(date_used == 'yes') %>% filter(date_type == "MC-ICP-MS U/Th" | date_type == "TIMS" | date_type ==  "ICP-MS U/Th Other" | date_type == "Alpha U/Th" | date_type ==  "U/Th unspecified")  %>% distinct(entity_id)
not_UTh_dates <- sisalv2$dating %>% filter(date_used == 'yes') %>% filter(date_type == 'Event; start of laminations' | date_type == 'Event; end of laminations' | date_type == 'C14' | date_type =='Multiple methods' | date_type =='other') %>%
  distinct(entity_id) 
UTh_dates <- sisalv2$dating %>% filter(date_used == 'yes') %>% filter(date_type == "MC-ICP-MS U/Th" | date_type == "TIMS" | date_type ==  "ICP-MS U/Th Other" | date_type == "Alpha U/Th" | date_type ==  "U/Th unspecified")  %>% distinct(entity_id) %>%
  filter(!(entity_id %in% not_UTh_dates$entity_id))

eId <- sisalv2$entity %>% filter(entity_status != 'superseded')

run <- sisalv2$dating %>% distinct(entity_id) %>%
    filter(entity_id %in% eId$entity_id) %>%
    filter(!(entity_id %in% no_depth_sample$entity_id)) %>%
    filter(!(entity_id %in% not_UTh_dates$entity_id)) %>%
    filter(entity_id %in% nr_dates$entity_id) %>%
    arrange(., entity_id)

write.csv(run, "~/Documents/partially.current.v2/run.csv", row.names = F)


no_chrono <- left_join(sample, orig, by = "sample_id") %>% filter(interp_age_uncert_pos == 'NULL') %>% distinct(entity_id)

dating_tb <- dating %>% left_join(.,entity, by = "entity_id") %>% filter(entity_status == 'current') %>%
  mutate_at(vars(dating_id, depth_dating, dating_thickness, X14C_correction, corr_age, corr_age_uncert_pos, corr_age_uncert_neg), as.numeric)
sample_tb <- sample %>% left_join(., hiatus, by = "sample_id") %>% left_join(., entity, by = "entity_id") %>% filter(entity_status == 'current') %>%
  mutate_at(vars(entity_id, sample_id, sample_thickness, depth_sample), as.numeric)

write.csv(class2, 'class2_test_run.csv', row.names = F)

hiatus <- sample_tb %>% filter(hiatus == 'H') %>% group_by(entity_id) %>% count()
reversals <- dating_tb %>% select(entity_id, depth_dating, corr_age, corr_age_uncert_pos, corr_age_uncert_neg) %>% group_by(entity_id) %>% 
  arrange(depth_dating, .by_group = T)
reversals_abs <- reversals %>% mutate(diff = lead(corr_age)-corr_age) %>% mutate(reversal = if_else(!is.na(diff) & diff < 0, TRUE, FALSE)) %>% group_by(entity_id) %>% #inner_join(., count(.), by = 'entity_id') %>% 
  arrange(depth_dating, .by_group = TRUE) %>% mutate(tractable = if_else(reversal & (abs(corr_age-lead(corr_age)) < (corr_age_uncert_pos + lead(corr_age_uncert_neg))), TRUE, FALSE))

r_reversal <- reversals_abs %>% select(entity_id, reversal) %>% group_by(entity_id) %>% filter(reversal) %>% count()
r_tractable <- reversals_abs %>% select(entity_id, tractable) %>% group_by(entity_id) %>% filter(tractable) %>% count()

r_nontractable <- reversals_abs %>% filter(reversal) %>%filter(!(entity_id %in% r_tractable$entity_id)) %>% distinct(entity_id)

two_dates <- dating_tb %>% select(entity_id, corr_age) %>% group_by(entity_id) %>%count() %>% filter(n==2)
#three_dates <- dating_tb %>% select(entity_id, corr_age) %>% group_by(entity_id) %>%count() %>% filter(n==3)
nr_dates <- dating_tb %>% select(entity_id, corr_age) %>% group_by(entity_id) %>%count() %>% filter(n>=3)

no_depth_sample <- sample_tb %>% group_by(entity_id) %>% summarise(depths = if_else(all(is.na(depth_sample)), FALSE, TRUE)) %>% filter(!depths)

run1 <- dating_tb%>% filter(date_used == 'yes')  %>%filter(entity_id %in% nr_dates$entity_id)  %>% filter(date_type == "MC-ICP-MS U/Th" | date_type == "TIMS" | date_type ==  "ICP-MS U/Th Other" | date_type == "Alpha U/Th" |
                                date_type ==  "U/Th unspecified") %>% filter(!(entity_id %in% r_nontractable$entity_id)) %>% filter(!(entity_id %in% no_depth_sample$entity_id)) %>% distinct(entity_id) %>% arrange(., entity_id)

run2 <- dating_tb%>% filter(date_used == 'yes') %>% filter((entity_id %in% hiatus$entity_id)) %>%filter(entity_id %in% nr_dates$entity_id) %>% filter(date_type == "MC-ICP-MS U/Th" | date_type == "TIMS" | date_type ==  "ICP-MS U/Th Other" | date_type == "Alpha U/Th" |
                                date_type ==  "U/Th unspecified") %>% filter(!(entity_id %in% r_nontractable$entity_id)) %>% filter(!(entity_id %in% no_depth_sample$entity_id)) %>% distinct(entity_id) %>% arrange(., entity_id)

r2 <- run1 %>% filter(entity_id %in% hiatus$entity_id)
class2 <- dating_tb%>% filter(date_used == 'yes') %>%filter(entity_id %in% nr_dates$entity_id) %>% filter(date_type == "MC-ICP-MS U/Th" | date_type == "TIMS" | date_type ==  "ICP-MS U/Th Other" | date_type == "Alpha U/Th" |
                               date_type ==  "U/Th unspecified") %>% filter(!(entity_id %in% hiatus$entity_id)) %>% filter(!(entity_id %in% r_nontractable$entity_id)) %>% distinct(entity_id)

class1 <- dating_tb %>% filter(date_type == 'C14' | date_type == 'Multiple methods' | (entity_id %in% hiatus$entity_id) | (entity_id %in% r_nontractable$entity_id)) %>% filter(!(entity_id %in% two_dates$entity_id))%>% distinct(entity_id)

class0 <- dating_tb %>% filter(entity_id %in% two_dates$entity_id) %>% distinct(entity_id)

duplicated <- class2 %>% filter(entity_id %in% class1$entity_id) 
dup <- dating_tb %>% filter(entity_id %in% duplicated$entity_id) %>% select(entity_id, date_type, date_used, corr_age)

d <- dating_tb %>% filter(entity_id %in% duplicated$entity_id)   


af <- dating %>% filter(date_type == 'Event; actively forming') %>% select(dating_id, entity_id, date_used, corr_age, date_type, depth_dating)
af2 <- dating %>% filter(entity_id %in% af$entity_id) %>% filter(depth_dating == 0) %>% group_by(entity_id) %>% count() %>% filter(n==2)
t <- dating %>% filter(entity_id %in% af2$entity_id) %>% filter(depth_dating == 0) %>% 
  select(entity_id, dating_id, date_used, date_type, corr_age, depth_dating) %>% filter(date_used == 'yes')
