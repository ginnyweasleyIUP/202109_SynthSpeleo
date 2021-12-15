dropped <- exp.ranking.final.new.2 %>% filter(!is.na(DROP))

dropped$comment <- c('badly dated', 'dates missing for orig. AM', 'date_used wrong', 'reversals', 'dates missing for orig. AM', 'Hiatus wrong', 'not enough dates',
                     'Hiatus wrong', 'dates missing for orig. AM', 'not to DROP', 'Hiatus wrong', 'dates missing fro orig. AM', 'dates missinf for orig. AM', 
                     'reversals','date_used wrong','actively growing?','REDO manually', 'REDO manually', 'date_used wrong','not to DROP','not enough dates', 
                     'date_used wrong', 'badly dated','not enough dates')


refs <- sisalv2$entity_link_reference %>% filter(entity_id %in% c(127,277,336,369,525,548,555,614))
data <- 

cites <- sisalv2$reference %>% filter(ref_id %in% refs$ref_id) %>% left_join(.,refs, by = 'ref_id') %>% left_join(., dropped, by = 'entity_id') %>% 
  select(entity_id, ref_id, citation, publication_DOI, comment) %>% left_join(., sisalv2$entity %>% select(entity_id, data_DOI_URL), by = 'entity_id')


##### correct dating table ---------------------------###

sisalv2interim <- read.SISAL.files('~/SISAL Data/interimv2_20191209/', '')


ad <- sisalv2$dating %>% filter(date_type == 'Event; actively forming') %>% select(-max) %>% ungroup() %>% slice(1) %>% mutate(dating_id == -1, entity_id = 548, corr_age = -56)

dating_upd <- read.csv('~/SISAL Data/sisalv2/dating_old2.csv', header = T, stringsAsFactors = F) %>%
  mutate_at(vars(depth_dating, dating_thickness,min_weight, max_weight, uncorr_age, uncorr_age_uncert_pos, uncorr_age_uncert_neg, starts_with('X'), corr_age), as.numeric) %>%
  mutate(date_used = if_else(entity_id == 525 & date_used == 'no', 'yes', date_used)) %>% # entity_id 525 change date_used to 'yes'
  mutate(depth_dating = if_else(entity_id == 614 & depth_dating == 223.9, 233.9, depth_dating)) %>% # entity_id 614 correct wrong depth dating
  mutate(depth_dating = if_else(entity_id == 369 & depth_dating == 163,160, depth_dating)) %>%
  bind_rows(., ad) %>% # add actively forming to entity id 548
  mutate(depth_dating = if_else(entity_id == 277 & depth_dating == 457.8, 480, depth_dating)) %>%
  mutate(date_used = if_else(entity_id == 277 & date_used == 'unknown', 'yes', date_used)) %>%
  mutate(date_used = if_else(entity_id == 582 & date_used == 'no', 'yes', date_used)) %>%
  mutate(date_used = if_else(entity_id == 582 & depth_dating == 30.3, 'no', date_used)) %>%
  mutate(date_used = if_else(entity_id == 561 & date_type == 'Event; hiatus', 'yes', date_used)) 

                
test <- sisalv2interim$dating %>% filter(entity_id %in% c(277, 369,614,525,548,582))
t48 <- sisalv2interim$dating %>% filter(entity_id == 48)

write.csv(dating_upd, '~/SISAL Data/sisalv2/dating.csv', row.names = F)

hiatus_upd <- read.csv('~/SISAL Data/sisalv2/sample.csv', header = T, stringsAsFactors = F) %>% mutate_at(vars(depth_sample), as.numeric) %>% mutate(depth_sample = if_else(entity_id == 277 & sample_id == 205227, 480, depth_sample),
                                                                                                       depth_sample = if_else(entity_id == 369 & sample_id == 196190, 160, depth_sample))
                                                                                                       
write.csv(hiatus_upd, '~/SISAL Data/sisalv2/sample.csv', row.names = F)
