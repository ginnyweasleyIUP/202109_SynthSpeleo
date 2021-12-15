library(tidyverse)

dating <- read.csv("~/SISAL Data/interim_v2_csv_20190812/dating.csv", header = T, stringsAsFactors = F)
lamina <- read.csv("~/SISAL Data/interim_v2_csv_20190812/dating_lamina.csv", header = T, stringsAsFactors = F)
composit <- read.csv("~/SISAL Data/interim_v2_csv_20190812/composite_link_entity.csv", header = T, stringsAsFactors = F)
entity <- read.csv("~/SISAL Data/interim_v2_csv_20190812/entity.csv", header = T, stringsAsFactors = F)
sample <- read.csv("~/SISAL Data/interim_v2_csv_20190812/sample.csv", header = T, stringsAsFactors = F)
orig <- read.csv("~/SISAL Data/interim_v2_csv_20190812/original_chronology.csv", header = T, stringsAsFactors = F)


sample_tb <- left_join(sample, orig)

holo <- left_join(entity, dating, by = "entity_id") %>% left_join(., lamina, by = 'entity_id') %>% mutate_at(vars(corr_age),as.numeric) %>% 
  filter(entity_status =='current' & !(entity_id %in% composit$composite_entity_id)) %>%
  group_by(entity_id) %>% filter(corr_age < 10000) %>% count() %>% filter(n> 3) %>% rename(n_holo = n)

lgm <- left_join(entity, dating, by = "entity_id") %>% left_join(., lamina, by = 'entity_id') %>% mutate_at(vars(corr_age),as.numeric) %>% 
  filter(entity_status =='current' & !(entity_id %in% composit$composite_entity_id)) %>%
  group_by(entity_id) %>% filter(corr_age > 19000 & corr_age<29000) %>% count() %>% filter(n>3) %>% rename(n_lgm = n)
  
all <- full_join(holo,lgm) %>% filter(!is.na(n_holo) & !is.na(n_lgm))
  
new <- all %>% full_join(., gr_full_hl)


summarise(min_age = min(corr_age, na.rm = T),
          max_age = max(corr_age, na.rm = T)) %>% filter(min_age != Inf) %>%
  filter(min_age < 10000 & max_age > 19000)


dating <- read.csv("~/SISAL Data/sisalv2/dating.csv", header = T, stringsAsFactors = F) 
entity <- read.csv("~/SISAL Data/sisalv2/entity.csv", header = T, stringsAsFactors = F)
site <- read.csv("~/SISAL Data/sisalv2/site.csv", header = T, stringsAsFactors = F)

site_tb <- left_join(site, entity, by = "site_id")
dating_tb <- left_join(dating, entity, by = "entity_id") %>% filter(entity_status == 'current') %>% mutate_at(vars(corr_age),as.numeric) %>% filter(date_used == 'yes')
p1 <- dating_tb %>% group_by(entity_id) %>% filter(corr_age <= 10000) %>% summarise(diff1 = max(corr_age, na.rm = T) - min(corr_age, na.rm = T),
                                                                                    n1 = n(), res1 = 1) %>% filter(n1 > 2 & diff1 >= 4000)
p2 <- dating_tb %>% group_by(entity_id) %>% filter(corr_age <= 20000 & corr_age > 10000) %>% summarise(diff2 = max(corr_age, na.rm = T) - min(corr_age, na.rm = T),
                                                                                                       n2 = n(), res2 = 1) %>% filter(n2 > 2 & diff2 >= 4000)
p3 <- dating_tb %>% group_by(entity_id) %>% filter(corr_age <= 30000 & corr_age > 20000) %>% summarise(diff3 = max(corr_age, na.rm = T) - min(corr_age, na.rm = T),
                                                                                                       n3 = n(), res3 = 1) %>% filter(n3 > 2 & diff3 >= 4000)
p4 <- dating_tb %>% group_by(entity_id) %>% filter(corr_age <= 40000 & corr_age > 30000) %>% summarise(diff4 = max(corr_age, na.rm = T) - min(corr_age, na.rm = T),
                                                                                                       n4= n(), res4 = 1) %>% filter(n4 > 2 & diff4 >= 4000)
p5 <- dating_tb %>% group_by(entity_id) %>% filter(corr_age <= 50000 & corr_age > 40000) %>% summarise(diff5 = max(corr_age, na.rm = T) - min(corr_age, na.rm = T),
                                                                                                       n5 = n(), res5 = 1) %>% filter(n5 > 2 & diff5 >= 4000)
p6 <- dating_tb %>% group_by(entity_id) %>% filter(corr_age <= 60000 & corr_age > 50000) %>% summarise(diff6 = max(corr_age, na.rm = T) - min(corr_age, na.rm = T),
                                                                                                       n6 = n(), res6 = 1) %>% filter(n6 > 2 & diff6 >= 4000)
p7 <- dating_tb %>% group_by(entity_id) %>% filter(corr_age <= 70000 & corr_age > 60000) %>% summarise(diff7 = max(corr_age, na.rm = T) - min(corr_age, na.rm = T),
                                                                                                       n7 = n(), res7 = 1) %>% filter(n7 > 2 & diff7 >= 4000)
p8 <- dating_tb %>% group_by(entity_id) %>% filter(corr_age <= 80000 & corr_age > 70000) %>% summarise(diff8 = max(corr_age, na.rm = T) - min(corr_age, na.rm = T),
                                                                                                       n8 = n(), res8 = 1) %>% filter(n8 > 2 & diff8 >= 4000)
p9 <- dating_tb %>% group_by(entity_id) %>% filter(corr_age <= 90000 & corr_age > 80000) %>% summarise(diff9 = max(corr_age, na.rm = T) - min(corr_age, na.rm = T),
                                                                                                       n9 = n(), res9 = 1) %>% filter(n9 > 2 & diff9 >= 4000)
p10 <- dating_tb %>% group_by(entity_id) %>% filter(corr_age <= 100000 & corr_age > 90000) %>% summarise(diff10 = max(corr_age, na.rm = T) - min(corr_age, na.rm = T),
                                                                                                       n10 = n(), res10 = 1) %>% filter(n10 > 2& diff10 >= 4000)
p11 <- dating_tb %>% group_by(entity_id) %>% filter(corr_age <= 110000 & corr_age > 100000) %>% summarise(diff11 = max(corr_age, na.rm = T) - min(corr_age, na.rm = T),
                                                                                                       n11 = n(), res11 = 1) %>% filter(n11 > 2 & diff11 >= 4000)
p12 <- dating_tb %>% group_by(entity_id) %>% filter(corr_age <= 120000 & corr_age > 110000) %>% summarise(diff12 = max(corr_age, na.rm = T) - min(corr_age, na.rm = T),
                                                                                                       n12 = n(), res12 = 1) %>% filter(n12 > 2 & diff12 >= 4000)
p13 <- dating_tb %>% group_by(entity_id) %>% filter(corr_age <= 130000 & corr_age > 120000) %>% summarise(diff13 = max(corr_age, na.rm = T) - min(corr_age, na.rm = T),
                                                                                                       n13 = n(), res13 = 1) %>% filter(n13 > 2 & diff13 >= 4000)

full <- plyr::join_all(list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13), by = 'entity_id', type = 'full') %>% group_by(entity_id) %>% 
  mutate(res = sum(c(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10,res11,res12,res13), na.rm = T)) %>% 
  filter(res > 1) %>% arrange(entity_id)

full_sites <- site_tb %>% filter(entity_id %in% full$entity_id) %>% select(site_id, site_name, entity_id, longitude, latitude)

sites <- site %>% select(site_id, longitude, latitude)  #%>% distinct(site_id)




