library(RMySQL)
library(rJava)
library(xlsx)
#library(dplyr)
library(tidyverse)
library(janitor)

# Connect to Database -----------------------------------------------####
mydb = dbConnect(MySQL(), user='root', 
                 password='Schrank8956', 
                 dbname='sisalv1b', 
                 host='localhost')


# tables
dating_tb <- dbGetQuery(mydb,"SELECT * FROM dating LEFT JOIN entity USING (entity_id) WHERE entity_status = 'current'")
original_tb <- dbGetQuery(mydb, "SELECT * FROM sample LEFT JOIN original_chronology USING (sample_id) LEFT JOIN entity USING (entity_id) WHERE entity_status = 'current'")
#sample_tb <- dbGetQuery(mydb, "SELECT * FROM sample LEFT JOIN entity USING (entity_id) WHERE entity_status = 'current'")
isotopes <- dbGetQuery(mydb, "SELECT * FROM d13C ")
sample_tb <- dbGetQuery(mydb, "SELECT * FROM sample LEFT JOIN hiatus USING (sample_id) LEFT JOIN entity USING (entity_id) WHERE entity_status ='current'")
hiatus <- dbGetQuery(mydb, 'SELECT * FROM hiatus')
composites <- dbGetQuery(mydb, 'SELECT * FROM composite_link_entity')
site_tb <- dbGetQuery(mydb, "SELECT * FROM site LEFT JOIN entity USING (site_id) LEFT JOIN entity_link_reference USING (entity_id) LEFT JOIN reference USING (ref_id) WHERE entity_status = 'current'")
dating_lamina <- dbGetQuery(mydb, "SELECT * FROM dating_lamina LEFT JOIN entity USING (entity_id) WHERE entity_status = 'current'")

# lamina dated
nr_lamina <- dating_lamina %>% distinct(entity_id)


# dublicated sample depths
dup_s <- sample_tb %>% filter(!(entity_id %in% composites$composite_entity_id) & !(entity_id %in% nr_lamina$entity_id)) %>%select(entity_id, sample_id, depth_sample) %>%
  group_by(entity_id) %>% get_dupes(depth_sample) %>% filter(!is.na(depth_sample)) %>% select(entity_id,depth_sample,sample_id)

dup_d <- dating_tb %>% filter(!(entity_id %in% composites$single_entity_id)& !(entity_id %in% nr_lamina$entity_id)) %>% select(entity_id, depth_dating, corr_age, date_used) %>% filter(date_used == 'yes') %>%
  group_by(entity_id) %>% get_dupes(depth_dating)  #%>% distinct(entity_id)

or <- original_tb %>% filter(sample_id %in% dup_s$sample_id) %>% select(entity_id,sample_id, depth_sample,interp_age)

# not marked as active
active_dating <- dating_tb %>%  select(entity_id, drip_type, date_used ,date_type, corr_age ) %>% filter(date_used == 'yes')%>% group_by(entity_id) %>% arrange(corr_age, .by_group = T) %>%slice(.,1)
active_original <- original_tb %>%  select(entity_id, interp_age) %>% group_by(entity_id) %>% arrange(interp_age, .by_group = T) %>%slice(.,1)

active <- inner_join(active_dating, active_original, by = 'entity_id') %>% filter(interp_age < corr_age & interp_age <= 0)

event_active <- dating_tb %>% filter(date_type == 'Event; actively forming') %>% select(entity_id, drip_type, date_used ,date_type, corr_age ) %>% filter(date_used == 'yes')%>% group_by(entity_id) %>% arrange(corr_age, .by_group = T) %>%slice(.,1)
event_active_orig <- original_tb %>% filter(entity_id %in% event_active$entity_id)  %>% select(entity_id, interp_age) %>% group_by(entity_id) %>% arrange(interp_age, .by_group = T) %>%slice(.,1)
event <- inner_join(event_active, event_active_orig, by = "entity_id") %>% select(entity_id, corr_age, interp_age,date_type, drip_type, date_used)

# missing hiatuses
dft_age <- original_tb %>% select(entity_id, depth_sample, interp_age)  %>% group_by(entity_id) %>% summarise(depth_sample_max = max(depth_sample, na.rm=T), depth_sample_min = min(depth_sample, na.rm = T), interp_age_max = max(interp_age,na.rm = T),
                                                                                                              interp_age_min = min(interp_age,na.rm = T)) %>% mutate(dft = depth_sample_max - depth_sample_min, age = interp_age_max - interp_age_min) 
gr <- dft_age %>% group_by(dft, age) %>% filter(age != Inf & age != -Inf & dft != Inf & dft != -Inf)%>% mutate(gr = dft/(age/1000)) %>% group_by(gr)

mean_gr <- original_tb %>% group_by(entity_id) %>% arrange(depth_sample, .by_group = T) %>% 
  mutate(growth_rate = (lead(depth_sample)-depth_sample)/(lead(interp_age)-interp_age)) 

ave_gr <- mean_gr %>% filter (!is.infinite(growth_rate)) %>% mutate(ave_growth_rate = mean(growth_rate, na.rm = T)) %>% 
  select(entity_id, entity_name, interp_age, depth_sample, growth_rate, ave_growth_rate)
missing_hiatus_entities <- ave_gr %>% group_by(entity_id) %>% filter(!is.na(growth_rate)) %>% filter(50*growth_rate < ave_growth_rate & growth_rate > 0) %>% count() %>% filter(n == 1)#%>% summarise(avg = mean(ave))
missing_hiatus <- ave_gr %>% filter(entity_id %in% missing_hiatus_entities$entity_id) %>% group_by(entity_id) %>% 
  mutate(lead_interp_age = lead(interp_age),
         lead_depth_sample = lead(depth_sample)) %>%
  filter(!is.na(growth_rate)) %>% filter(50*growth_rate < ave_growth_rate & growth_rate > 0) #%>%
  

# wrong order in depth
wrong_order <- ave_gr %>% group_by(entity_id) %>% filter(ave_growth_rate < 0) %>% summarise(av_growth_rate = mean(ave_growth_rate))
dating_wrong_order <- dating_tb %>% group_by(entity_id) %>%arrange(depth_dating, .by_group = T) %>% select(entity_id, entity_name ,depth_dating, corr_age, date_type, date_used, depth_ref) %>% 
  filter(date_used == 'yes' & entity_id %in% wrong_order$entity_id & depth_ref == 'from top' & entity_id != 314)

setwd( "/Users/Carlii/Documents/Masterarbeit/SISAL Errors")

write.csv(or, "duplicated_sample_depths.csv", row.names = F)
write.csv(active, "not_marked_active.csv", row.names = F)
write.csv(missing_hiatus, "possibly_missing_hiatus.csv", row.names = F)
write.csv(dating_wrong_order, "wrong_depth_order.csv", row.names = F)


#unuseable orig. chronologies
unuseable_orig_chrono <- dating_tb %>% filter(!(entity_id %in% composites$single_entity_id)) %>% select(entity_id, depth_dating, corr_age, date_used, date_type, calib_used,speleothem_type) %>% 
  filter(date_used == 'yes' & date_type == 'C14')

o <- original_tb %>% filter(entity_id == 152)
d <- dating_tb %>% filter(entity_id == 152 & date_used == 'yes')
plot(o$interp_age, o$depth_sample)
points(d$corr_age, d$depth_dating, col = 'orange', cex = 2, lwd = 4)

o <- original_tb %>% filter(entity_id == 181)
d <- dating_tb %>% filter(entity_id == 181 & date_used == 'yes')
plot(o$interp_age, o$depth_sample)
points(d$corr_age, d$depth_dating, col = 'orange', cex = 2, lwd = 4)


o <- original_tb %>% filter(entity_id == 185)
d <- dating_tb %>% filter(entity_id == 185 & date_used == 'yes')
plot(o$interp_age, o$depth_sample)
points(d$corr_age, d$depth_dating, col = 'orange', cex = 2, lwd = 4)


o <- original_tb %>% filter(entity_id == 210)
d <- dating_tb %>% filter(entity_id == 210 & date_used == 'yes')
plot(o$interp_age, o$depth_sample)
points(d$corr_age, d$depth_dating, col = 'orange', cex = 2, lwd = 4)

o <- original_tb %>% filter(entity_id == 210)
d <- dating_tb %>% filter(entity_id == 210 & date_used == 'yes')
plot(o$interp_age, o$depth_sample)
points(d$corr_age, d$depth_dating, col = 'orange', cex = 2, lwd = 4)

o <- original_tb %>% filter(entity_id == 382)
d <- dating_tb %>% filter(entity_id == 382 & date_used == 'yes')
plot(o$interp_age, o$depth_sample)
points(d$corr_age, d$depth_dating, col = 'orange', cex = 2, lwd = 4)

