gr <- growth.rates(ens, b=T)

run<- read.csv('/stacywork/ariana/v2_90ci/run2.csv', header = T)
runFile <- tibble(entity_id = run$entity_id, Bacon = rep(T, length(run$entity_id)), Bchron = rep(T, length(run$entity_id)), 
                  StalAge = rep(T, length(run$entity_id)), linInterp = rep(T, length(run$entity_id)), copRa = rep(T, length(run$entity_id)),
                  linReg = rep(T, length(run$entity_id)), 
                  working_directory = rep("/stacywork/ariana/v2_new/", length(run$entity_id)), 
                  wd = rep('/home/ariana/SISAL Data/sisalv2', length(run$entity_id))) 
runFile <- runFile %>% mutate(Bchron = if_else(entity_id %in% c(115,330,331,378,432,506,510), F, T)) %>% arrange(., entity_id)  # 2: ab 115


run1<- read.csv('~/Documents/partially.current.v2/run.csv', header = T)
runFile1 <- tibble(entity_id = run1$entity_id, Bacon = rep(T, length(run1$entity_id)), Bchron = rep(T, length(run1$entity_id)), 
                  StalAge = rep(T, length(run1$entity_id)), linInterp = rep(T, length(run1$entity_id)), copRa = rep(T, length(run1$entity_id)),
                  linReg = rep(T, length(run1$entity_id)), 
                  working_directory = rep("/stacywork/ariana/v2_new/", length(run1$entity_id)), 
                  wd = rep('/home/ariana/SISAL Data/sisalv2', length(run1$entity_id))) 

runFile_new <- bind_rows(runFile,runFile1) %>% arrange(entity_id)

r1 <- merge_SISAL_ensemble_final(runFile, dating_from_base)
r <- merge_SISAL_ensemble_final(runFile, dating_from_base)
r <- merge_S

o <- oxcal_final %>% distinct(entity_id) %>% filter(!(entity_id %in% missing$entity_id))

r_new <- left_join(r, exp.lk, by = 'entity_id') %>% left_join(., k %>% select(entity_id, lR, oxcal)) %>%
  mutate(Bacon = if_else(Bacon.x & is.na(Bacon.y), TRUE, FALSE),
         Bchron = if_else(Bchron.x & is.na(Bchron.y), TRUE, FALSE),
         StalAge = if_else(StalAge.x & is.na(StalAge.y), TRUE, FALSE),
         linInterp = if_else(linInterp & is.na(LI), TRUE, FALSE),
         linReg = if_else(linReg & is.na(LR) & lR != 0 & !is.na(lR), TRUE, FALSE),
         copRa = if_else(copRa.x & is.na(copRa.y), TRUE, FALSE),
         OxCal = if_else(entity_id %in% o$entity_id & oxcal != 0 & !is.na(oxcal), TRUE, FALSE )) %>%
  select(entity_id, Bacon, Bchron, StalAge, linInterp, linReg, copRa, OxCal, working_directory, wd)

r1_new <- left_join(r1, exp.lk, by = 'entity_id') %>% left_join(., k %>% select(entity_id, lR, oxcal)) %>%
  mutate(Bacon = if_else(Bacon.x & is.na(Bacon.y), TRUE, FALSE),
         Bchron = if_else(Bchron.x & is.na(Bchron.y), TRUE, FALSE),
         StalAge = if_else(StalAge.x & is.na(StalAge.y), TRUE, FALSE),
         linInterp = if_else(linInterp & is.na(LI), TRUE, FALSE),
         linReg = if_else(linReg & is.na(LR) & lR != 0 & !is.na(lR), TRUE, FALSE),
         copRa = if_else(copRa.x & is.na(copRa.y), TRUE, FALSE),
         OxCal = if_else(entity_id %in% o$entity_id & oxcal != 0 & !is.na(oxcal), TRUE, FALSE )) %>%
  select(entity_id, Bacon, Bchron, StalAge, linInterp, linReg, copRa, OxCal, working_directory, wd)

final.plot(r_new, entity, dating_from_base, drop)
final.plot(r1_new,entity, dating_from_base, drop)


final <- merge_SISAL_chrono_final(runFile_new,se,entity)
r_f <- final[[1]]

s <- final[[2]] 
oxcal_new <- read.csv("~/SISAL/Oxcal/SISAL_v2_Oxcal_Chronology_LCB_20191103.csv", header = T, stringsAsFactors = F) 
oxcal <- read.csv("~/SISAL/Oxcal/SISAL_v2_Oxcal_Chronology.csv", header = T, stringsAsFactors = F) 
oxcal_final <- left_join(oxcal_new, oxcal %>% select(sample_id, depth), by = 'sample_id')
colnames(oxcal_final) <- c('entity_id', 'sample_id','OxCal_age', 'OxCal_age_uncert_pos', 'OxCal_age_uncert_neg', 'depth')
oxcal_dating <- read.csv("~/SISAL/Oxcal/SISAL_v2_Oxcal_Dating.csv", header = T, stringsAsFactors = F) 
s_fil <- s %>% filter(entity_id %in% runFile_new$entity_id)
full_chrono <- full_join(s_fil, oxcal_final, by = c('sample_id', 'entity_id')) %>% select(-depth)

eval <- eval(full_chrono, sisalv2$dating, sisalv2$hiatus)

SISAL_eval <- eval[[5]] #%>% select(-OxCal)
k <- SISAL_eval %>% filter(lR == 0 | oxcal == 0)
k.lr <- k %>% filter(lR == 0)
k.oxcal <- k %>% filter(oxcal == 0)
s_fil_new <- s_fil %>% mutate(lin_reg_age = if_else(entity_id %in% k.lr$entity_id, NA_real_, lin_reg_age),
                              lin_reg_age_uncert_pos = if_else(entity_id %in% k.lr$entity_id, NA_real_, lin_reg_age_uncert_pos),
                              lin_reg_age_uncert_neg = if_else(entity_id %in% k.lr$entity_id, NA_real_, lin_reg_age_uncert_neg))

full_chrono_eval <- full_chrono %>% mutate(lin_reg_age = if_else(entity_id %in% k.lr$entity_id, NA_real_, lin_reg_age),
                         lin_reg_age_uncert_pos = if_else(entity_id %in% k.lr$entity_id, NA_real_, lin_reg_age_uncert_pos),
                         lin_reg_age_uncert_neg = if_else(entity_id %in% k.lr$entity_id, NA_real_, lin_reg_age_uncert_neg),
                         OxCal_age = if_else(entity_id %in% k.oxcal$entity_id, NA_real_, OxCal_age),
                         OxCal_age_uncert_pos = if_else(entity_id %in% k.oxcal$entity_id, NA_real_, OxCal_age_uncert_pos),
                         OxCal_age_uncert_neg = if_else(entity_id %in% k.oxcal$entity_id, NA_real_, OxCal_age_uncert_neg))


write.csv(s_fil_new, '/home/ariana/SISAL/sisal_chronology_iup.csv', row.names = F)
write.csv(full_chrono, '/home/ariana/SISAL/sisal_chronology.csv', row.names = F)

dating_new <- full_join(dating_from_base, oxcal_dating, by = c('entity_id','dating_id'))

plot.oxcal(oxcal_final, entity,dating_new)

### tables for entities where NAs have to be included around hiatuses ----------------####
hiatus_tb <- sisalv2$dating %>% filter(date_type == 'Event; hiatus') %>% filter(date_used == 'yes') %>% group_by(entity_id) %>% arrange(depth_dating, .by_group = T) %>% slice(.,1) %>% select(entity_id, dating_id, depth_dating) %>% rename(hiatus_depth = depth_dating)
hiatus_top <-sisalv2$dating %>% filter(date_used == 'yes') %>% filter(entity_id %in% hiatus_tb$entity_id) %>% group_by(entity_id) %>% arrange(depth_dating, .by_group = T) %>% select(entity_id, dating_id, depth_dating, date_type) %>%
  mutate(hiatus_depth = if_else(dating_id %in% hiatus_tb$dating_id,1,0)) %>% 
  mutate(remove = cumsum(cumsum(hiatus_depth))) %>%
  filter(remove <= 1) %>%
  print()

top <- sisalv2$dating %>% filter(date_used == 'yes') %>% group_by(entity_id) %>% arrange(., depth_dating, .by_group = T)%>% mutate(h = 1, remove = cumsum(h)) %>% filter(date_type == 'Event; hiatus' & remove <=2) %>% filter(entity_id %in% eID$entity_id)
bottom <- sisalv2$dating %>% filter(date_used == 'yes') %>% group_by(entity_id) %>% arrange(., desc(depth_dating), .by_group = T)%>% mutate(h = 1, remove = cumsum(h)) %>% filter(date_type == 'Event; hiatus' & remove <=2) %>% filter(entity_id %in% eID$entity_id)

rm_hiatus <- hiatus_top %>% count() %>% filter(n < 3) %>% left_join(., hiatus_tb) #%>% mutate(dating_id = -dating_id)### entities where there is only one date before/after a hiatus
#bloed <- rm_dt %>% filter(diff < 0)

rm_hiatus_dt <- sisalv2$dating %>% filter(date_used == 'yes') %>% mutate(remove_h = if_else(dating_id %in% c(rm_hiatus$dating_id,bottom$dating_id),1,
                                                                                              if_else(lead(dating_id) %in% top$dating_id,1,
                                                                                                      if_else(lag(dating_id) %in% top$dating_id,1,
                                                                                                              if_else(lead(dating_id) %in% bottom$dating_id,1,
                                                                                                                      if_else(lag(dating_id) %in% bottom$dating_id, 1,0)))))) %>% 
  mutate(dating_id = -dating_id,
         remove = cumsum(remove_h))

rm_dt <- sisalv2$dating %>% filter(date_used == 'yes') %>% mutate(dating_id = -dating_id) %>% group_by(entity_id) %>% arrange(depth_dating, .by_group = T) %>% filter(!(dating_id %in% c(rm_hiatus$dating_id, bottom$dating_id))) %>% select(entity_id, dating_id, depth_dating, date_type, corr_age) %>%
  mutate(r = if_else(date_type == 'Event; hiatus',2,
                     if_else(lag(date_type) == 'Event; hiatus', -1,
                             if_else(lead(date_type) == 'Event; hiatus',1,0)))) %>%
  filter(r %in% c(1,-1)) %>%
  mutate(diff = lead(corr_age) - corr_age) %>%
  filter(r == 1) %>%
  filter(diff > 2000)

help <- sisalv2$dating %>% filter(date_used == 'yes') %>% mutate(dating_id = -dating_id) %>% group_by(entity_id) %>% arrange(depth_dating, .by_group = T) %>% 
  filter(!(dating_id %in% c(rm_hiatus$dating_id, bottom$dating_id))) %>% select(entity_id, dating_id, depth_dating, date_type, corr_age) %>%
  mutate(r = if_else(lag(lag(dating_id)) %in% rm_dt$dating_id, 1,0)) %>% filter(r == 1)
### tables for entities where NAs have to be included around hiatuses ----------------####




##### correct sisal chronology #######
remove <- read.csv('~/SISAL/updates_to_chrono/remove_agreed_kira_laia_2019-11-03.csv', header = T)


dt_h <- sisalv2$dating %>% filter(date_used == 'yes') %>% group_by(entity_id) %>% arrange(., depth_dating, .by_group = T) %>%
  select(entity_id, dating_id, depth_dating) %>%
  mutate(remove = if_else(dating_id %in% rm_hiatus$dating_id,1,
                          if_else(dating_id %in% bottom$dating_id, -1,
                          if_else(lead(dating_id) %in% top$dating_id,1,
                                  if_else(lag(dating_id) %in% top$dating_id,1,
                                          if_else(lead(dating_id) %in% bottom$dating_id,-1,
                                                  if_else(lag(dating_id) %in% bottom$dating_id, -1,0))))))) %>% 
  mutate(dating_id = -dating_id) %>% 
  rename(depth_sample = depth_dating, sample_id = dating_id) %>% ungroup() 
dt_h <- dt_h %>%
  mutate(lin_reg_age = rep(NA_real_, length(dt_h$sample_id)),
         lin_reg_age_uncert_pos = rep(NA_real_, length(dt_h$sample_id)),
         lin_reg_age_uncert_neg = rep(NA_real_, length(dt_h$sample_id)),
         lin_interp_age = rep(NA_real_, length(dt_h$sample_id)),
         lin_interp_age_uncert_pos = rep(NA_real_, length(dt_h$sample_id)),
         lin_interp_age_uncert_neg = rep(NA_real_, length(dt_h$sample_id)),
         copRa_age = rep(NA_real_, length(dt_h$sample_id)),
         copRa_age_uncert_pos = rep(NA_real_, length(dt_h$sample_id)),
         copRa_age_uncert_neg = rep(NA_real_, length(dt_h$sample_id)),
         bacon_age = rep(NA_real_, length(dt_h$sample_id)),
         bacon_age_uncert_pos = rep(NA_real_, length(dt_h$sample_id)),
         bacon_age_uncert_neg = rep(NA_real_, length(dt_h$sample_id)),
         bchron_age = rep(NA_real_, length(dt_h$sample_id)),
         bchron_age_uncert_pos = rep(NA_real_, length(dt_h$sample_id)),
         bchron_age_uncert_neg = rep(NA_real_, length(dt_h$sample_id)),
         StalAge_age = rep(NA_real_, length(dt_h$sample_id)),
         StalAge_age_uncert_pos = rep(NA_real_, length(dt_h$sample_id)),
         StalAge_age_uncert_neg = rep(NA_real_, length(dt_h$sample_id))) %>%
  mutate_at(vars(everything()), as.character) %>%
  mutate_at(vars(everything()), as.numeric)


full_chrono_new <- bind_rows(full_chrono_eval,dt_h) %>% group_by(entity_id) %>% arrange(depth_sample, .by_group = T) %>% filter(!(sample_id %in% hiatus$sample_id)) %>%
  mutate(remove = if_else(is.na(remove),0,remove)) %>%
  mutate(remove_new = cumsum(remove)) %>% 
  mutate(lin_reg_age = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                                 remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                                 remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                                 remove_new >= 0 & entity_id %in% bottom$entity_id ~lin_reg_age,
                                 remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ lin_reg_age,
                                 remove_new >= 3 & entity_id %in% top$entity_id ~ lin_reg_age,
                                 !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ lin_reg_age),
         lin_reg_age_uncert_pos = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                                            remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                                            remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                                            remove_new >= 0 & entity_id %in% bottom$entity_id ~lin_reg_age_uncert_pos,
                                            remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ lin_reg_age_uncert_pos,
                                            remove_new >= 3 & entity_id %in% top$entity_id ~ lin_reg_age_uncert_pos,
                                            !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ lin_reg_age_uncert_pos),
         lin_reg_age_uncert_neg = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                                            remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                                            remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                                            remove_new >= 0 & entity_id %in% bottom$entity_id ~lin_reg_age_uncert_neg,
                                            remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ lin_reg_age_uncert_neg,
                                            remove_new >= 3 & entity_id %in% top$entity_id ~ lin_reg_age_uncert_neg,
                                            !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ lin_reg_age_uncert_neg),
         lin_interp_age = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                                    remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                                    remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                                    remove_new >= 0 & entity_id %in% bottom$entity_id ~lin_interp_age,
                                    remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ lin_interp_age,
                                    remove_new >= 3 & entity_id %in% top$entity_id ~ lin_interp_age,
                                    !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ lin_interp_age),
         lin_interp_age_uncert_pos = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                                               remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                                               remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                                               remove_new >= 0 & entity_id %in% bottom$entity_id ~lin_interp_age_uncert_pos,
                                               remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ lin_interp_age_uncert_pos,
                                               remove_new >= 3 & entity_id %in% top$entity_id ~ lin_interp_age_uncert_pos,
                                               !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ lin_interp_age_uncert_pos),
         lin_interp_age_uncert_neg = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                                               remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                                               remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                                               remove_new >= 0 & entity_id %in% bottom$entity_id ~lin_interp_age_uncert_neg,
                                               remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ lin_interp_age_uncert_neg,
                                               remove_new >= 3 & entity_id %in% top$entity_id ~ lin_interp_age_uncert_neg,
                                               !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ lin_interp_age_uncert_neg),
         copRa_age = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                               remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                               remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                               remove_new >= 0 & entity_id %in% bottom$entity_id ~copRa_age,
                               remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ copRa_age,
                               remove_new >= 3 & entity_id %in% top$entity_id ~ copRa_age,
                               !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ copRa_age),
         copRa_age_uncert_pos = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                                          remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                                          remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                                          remove_new >= 0 & entity_id %in% bottom$entity_id ~copRa_age_uncert_pos,
                                          remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ copRa_age_uncert_pos,
                                          remove_new >= 3 & entity_id %in% top$entity_id ~ copRa_age_uncert_pos,
                                          !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ copRa_age_uncert_pos),
         copRa_age_uncert_neg = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                                          remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                                          remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                                          remove_new >= 0 & entity_id %in% bottom$entity_id ~copRa_age_uncert_neg,
                                          remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ copRa_age_uncert_neg,
                                          remove_new >= 3 & entity_id %in% top$entity_id ~ copRa_age_uncert_neg,
                                          !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ copRa_age_uncert_neg),
         StalAge_age = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                                 remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                                 remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                                 remove_new >= 0 & entity_id %in% bottom$entity_id ~ StalAge_age,
                                 remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ StalAge_age,
                                 remove_new >= 3 & entity_id %in% top$entity_id ~ StalAge_age,
                                 !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ StalAge_age),
         StalAge_age_uncert_pos = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                                            remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                                            remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                                            remove_new >= 0 & entity_id %in% bottom$entity_id ~ StalAge_age_uncert_pos,
                                            remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ StalAge_age_uncert_pos,
                                            remove_new >= 3 & entity_id %in% top$entity_id ~ StalAge_age_uncert_pos,
                                            !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ StalAge_age_uncert_pos),
         StalAge_age_uncert_neg = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                                            remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                                            remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                                            remove_new >= 0 & entity_id %in% bottom$entity_id ~ StalAge_age_uncert_neg,
                                            remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ StalAge_age_uncert_neg,
                                            remove_new >= 3 & entity_id %in% top$entity_id ~ StalAge_age_uncert_neg,
                                            !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ StalAge_age_uncert_neg),
         bacon_age = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                               remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                               remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                               remove_new >= 0 & entity_id %in% bottom$entity_id ~ bacon_age,
                               remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ bacon_age,
                               remove_new >= 3 & entity_id %in% top$entity_id ~ bacon_age,
                               !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ bacon_age),
         bacon_age_uncert_pos = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                                          remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                                          remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                                          remove_new >= 0 & entity_id %in% bottom$entity_id ~ bacon_age_uncert_pos,
                                          remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ bacon_age_uncert_pos,
                                          remove_new >= 3 & entity_id %in% top$entity_id ~ bacon_age_uncert_pos,
                                          !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ bacon_age_uncert_pos),
         bacon_age_uncert_neg = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                                          remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                                          remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                                          remove_new >= 0 & entity_id %in% bottom$entity_id ~ bacon_age_uncert_neg,
                                          remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ bacon_age_uncert_neg,
                                          remove_new >= 3 & entity_id %in% top$entity_id ~ bacon_age_uncert_neg,
                                          !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ bacon_age_uncert_neg),
         bchron_age = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                                remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                                remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                                remove_new >= 0 & entity_id %in% bottom$entity_id ~ bchron_age,
                                remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ bchron_age,
                                remove_new >= 3 & entity_id %in% top$entity_id ~ bchron_age,
                                !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ bchron_age),
         bchron_age_uncert_pos = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                                           remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                                           remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                                           remove_new >= 0 & entity_id %in% bottom$entity_id ~ bchron_age_uncert_pos,
                                           remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ bchron_age_uncert_pos,
                                           remove_new >= 3 & entity_id %in% top$entity_id ~ bchron_age_uncert_pos,
                                           !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ bchron_age_uncert_pos),
         bchron_age_uncert_neg = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                                           remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                                           remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                                           remove_new >= 0 & entity_id %in% bottom$entity_id ~ bchron_age_uncert_neg,
                                           remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ bchron_age_uncert_neg,
                                           remove_new >= 3 & entity_id %in% top$entity_id ~ bchron_age_uncert_neg,
                                           !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ bchron_age_uncert_neg),
         OxCal_age = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                               remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                               remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                               remove_new >= 0 & entity_id %in% bottom$entity_id ~ OxCal_age,
                               remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ OxCal_age,
                               remove_new >= 3 & entity_id %in% top$entity_id ~ OxCal_age,
                               !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ OxCal_age),
         OxCal_age_uncert_pos = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                                          remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                                          remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                                          remove_new >= 0 & entity_id %in% bottom$entity_id ~ OxCal_age_uncert_pos,
                                          remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ OxCal_age_uncert_pos,
                                          remove_new >= 3 & entity_id %in% top$entity_id ~ OxCal_age_uncert_pos,
                                          !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ OxCal_age_uncert_pos),
         OxCal_age_uncert_neg = case_when(remove_new < 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ NA_real_,
                                          remove_new < 3 & entity_id %in% top$entity_id ~ NA_real_,
                                          remove_new < 0 & entity_id %in% bottom$entity_id ~ NA_real_,
                                          remove_new >= 0 & entity_id %in% bottom$entity_id ~ OxCal_age_uncert_neg,
                                          remove_new >= 3 & entity_id %in% top$entity_id & entity_id %in% bottom$entity_id ~ OxCal_age_uncert_neg,
                                          remove_new >= 3 & entity_id %in% top$entity_id ~ OxCal_age_uncert_neg,
                                          !(entity_id %in% c(top$entity_id, bottom$entity_id)) ~ OxCal_age_uncert_neg)) %>% 
  filter(!(sample_id %in% dt_h$sample_id))

  
  
dt <- dating_from_base %>% filter(date_type != 'Event; hiatus') %>% mutate(dating_id = -dating_id) %>% filter(date_used == 'yes') %>% select(entity_id, dating_id, depth_dating) %>% rename(depth_sample = depth_dating, sample_id = dating_id) %>% ungroup() 
dt <- dt %>%
  mutate(lin_reg_age = rep(NA_real_, length(dt$sample_id)),
         lin_reg_age_uncert_pos = rep(NA_real_, length(dt$sample_id)),
         lin_reg_age_uncert_neg = rep(NA_real_, length(dt$sample_id)),
         lin_interp_age = rep(NA_real_, length(dt$sample_id)),
         lin_interp_age_uncert_pos = rep(NA_real_, length(dt$sample_id)),
         lin_interp_age_uncert_neg = rep(NA_real_, length(dt$sample_id)),
         copRa_age = rep(NA_real_, length(dt$sample_id)),
         copRa_age_uncert_pos = rep(NA_real_, length(dt$sample_id)),
         copRa_age_uncert_neg = rep(NA_real_, length(dt$sample_id)),
         bacon_age = rep(NA_real_, length(dt$sample_id)),
         bacon_age_uncert_pos = rep(NA_real_, length(dt$sample_id)),
         bacon_age_uncert_neg = rep(NA_real_, length(dt$sample_id)),
         bchron_age = rep(NA_real_, length(dt$sample_id)),
         bchron_age_uncert_pos = rep(NA_real_, length(dt$sample_id)),
         bchron_age_uncert_neg = rep(NA_real_, length(dt$sample_id)),
         StalAge_age = rep(NA_real_, length(dt$sample_id)),
         StalAge_age_uncert_pos = rep(NA_real_, length(dt$sample_id)),
         StalAge_age_uncert_neg = rep(NA_real_, length(dt$sample_id))) %>%
  mutate_at(vars(everything()), as.character) %>%
  mutate_at(vars(everything()), as.numeric)

hiatus_removed <- full_chrono_eval %>% filter(sample_id %in% hiatus$sample_id)

full_chrono_new_new <- bind_rows(full_chrono_new, dt, hiatus_removed) %>% select(-remove_new) %>% group_by(entity_id) %>% arrange(depth_sample, .by_group = T) %>%
  mutate(remove = if_else(sample_id %in% rm_dt$dating_id, 1, 
                          if_else(sample_id %in% hiatus$sample_id,2,
                                  if_else(sample_id %in% help$dating_id,-3,0))),
         r = cumsum(remove)) %>%
  mutate(lin_reg_age = if_else(r%%2 != 0, NA_real_,lin_reg_age),
         lin_reg_age_uncert_pos = if_else(r%%2 != 0, NA_real_,lin_reg_age_uncert_pos),
         lin_reg_age_uncert_neg = if_else(r%%2 != 0, NA_real_,lin_reg_age_uncert_neg),
         lin_interp_age = if_else(r%%2 != 0, NA_real_,lin_interp_age),
         lin_interp_age_uncert_pos = if_else(r%%2 != 0, NA_real_,lin_interp_age_uncert_pos),
         lin_interp_age_uncert_neg = if_else(r%%2 != 0, NA_real_,lin_interp_age_uncert_neg),
         copRa_age = if_else(r%%2 != 0, NA_real_,copRa_age),
         copRa_age_uncert_pos = if_else(r%%2 != 0, NA_real_,copRa_age_uncert_pos),
         copRa_age_uncert_neg = if_else(r%%2 != 0, NA_real_,copRa_age_uncert_neg),
         StalAge_age = if_else(r%%2 != 0, NA_real_,StalAge_age),
         StalAge_age_uncert_pos = if_else(r%%2 != 0, NA_real_,StalAge_age_uncert_pos),
         StalAge_age_uncert_neg = if_else(r%%2 != 0, NA_real_,StalAge_age_uncert_neg),
         bacon_age = if_else(r%%2 != 0, NA_real_,bacon_age),
         bacon_age_uncert_pos = if_else(r%%2 != 0, NA_real_,bacon_age_uncert_pos),
         bacon_age_uncert_neg = if_else(r%%2 != 0, NA_real_,bacon_age_uncert_neg),
         bchron_age = if_else(r%%2 != 0, NA_real_,bchron_age),
         bchron_age_uncert_pos = if_else(r%%2 != 0, NA_real_,bchron_age_uncert_pos),
         bchron_age_uncert_neg = if_else(r%%2 != 0, NA_real_,bchron_age_uncert_neg),
         OxCal_age = if_else(r%%2 != 0, NA_real_,OxCal_age),
         OxCal_age_uncert_pos = if_else(r%%2 != 0, NA_real_,OxCal_age_uncert_pos),
         OxCal_age_uncert_neg = if_else(r%%2 != 0, NA_real_,OxCal_age_uncert_neg)) %>%
  filter(!(sample_id %in% dt$sample_id)) %>% select(-c(remove, r)) 

t <- full_chrono_new_new %>% group_by(entity_id) %>% slice(.,1)

full_chrono_new_new_new <- full_chrono_new_new %>% mutate(lin_reg_age = if_else(lin_reg_age < -68, NA_real_,lin_reg_age),
                                                          lin_reg_age_uncert_pos = if_else(lin_reg_age < -68, NA_real_,lin_reg_age_uncert_pos),
                                                          lin_reg_age_uncert_neg = if_else(lin_reg_age < -68, NA_real_,lin_reg_age_uncert_neg),
                                                          lin_interp_age = if_else(lin_interp_age < -68, NA_real_,lin_interp_age),
                                                          lin_interp_age_uncert_pos = if_else(lin_interp_age < -68, NA_real_,lin_interp_age_uncert_pos),
                                                          lin_interp_age_uncert_neg = if_else(lin_interp_age < -68, NA_real_,lin_interp_age_uncert_neg),
                                                          copRa_age = if_else(copRa_age < -68, NA_real_,copRa_age),
                                                          copRa_age_uncert_pos = if_else(copRa_age < -68, NA_real_,copRa_age_uncert_pos),
                                                          copRa_age_uncert_neg = if_else(copRa_age < -68, NA_real_,copRa_age_uncert_neg),
                                                          StalAge_age = if_else(StalAge_age < -68, NA_real_,StalAge_age),
                                                          StalAge_age_uncert_pos = if_else(StalAge_age < -68, NA_real_,StalAge_age_uncert_pos),
                                                          StalAge_age_uncert_neg = if_else(StalAge_age < -68, NA_real_,StalAge_age_uncert_neg),
                                                          bacon_age = if_else(bacon_age < -68, NA_real_,bacon_age),
                                                          bacon_age_uncert_pos = if_else(bacon_age < -68, NA_real_,bacon_age_uncert_pos),
                                                          bacon_age_uncert_neg = if_else(bacon_age < -68, NA_real_,bacon_age_uncert_neg),
                                                          bchron_age = if_else(bchron_age < -68, NA_real_,bchron_age),
                                                          bchron_age_uncert_pos = if_else(bchron_age < -68, NA_real_,bchron_age_uncert_pos),
                                                          bchron_age_uncert_neg = if_else(bchron_age < -68, NA_real_,bchron_age_uncert_neg),
                                                          OxCal_age = if_else(OxCal_age < -68, NA_real_,OxCal_age),
                                                          OxCal_age_uncert_pos = if_else(OxCal_age < -68, NA_real_,OxCal_age_uncert_pos),
                                                          OxCal_age_uncert_neg = if_else(OxCal_age < -68, NA_real_,OxCal_age_uncert_neg))

te <- full_chrono_new_new_new %>% group_by(entity_id) %>% slice(.,1)


write.csv(full_chrono_new_new, '~/SISAL/sisal_chronology_new.csv', row.names = F)

### expert Laia Kira ranking -----------------------------------####

library(tidyverse)
library(readxl)

exp.kira <- read_xlsx('/home/ariana/SISAL/expert/experts/expert_eval_kr.xlsx', col_types=c(rep("numeric",8),"text"))
exp.laia <- read_xlsx('/home/ariana/SISAL/expert/experts/rankings_laia.xlsx', col_types=c(rep("numeric",7),"text"))
colnames(exp.laia) <- c('EID','LR','LI','copRa','bacon','stalage','bchron','comments')
colnames(exp.kira) <- c(colnames(exp.kira)[1:8],'comments')

exp.laia.new <- exp.laia %>% mutate(DROP = if_else(LR + LI + copRa +bacon + stalage + bchron == 6, 1, NA_real_))

exp.lk <- full_join(exp.kira, exp.laia.new, by = 'EID') %>% select(-c(comments.x, comments.y)) %>% 
  mutate(copRa = copRa.x + copRa.y,
         StalAge = stalage.x + stalage.y,
         LI = LI.x + LI.y,
         LR = LR.x + LR.y,
         Bacon = bacon.x + bacon.y,
         Bchron = bchron.x + bchron.y,
         DROP = DROP.x + DROP.y) %>% 
  select(EID, copRa, StalAge, LI, LR, Bacon, Bchron, DROP) %>% 
  mutate(Bchron = if_else(EID == 614, 2, Bchron)) %>% 
  rename(entity_id = EID)

missing <- oxcal %>% filter(!(entity_id %in% s_fil$entity_id)) %>% distinct(entity_id)
sisal.ext <- read_xlsx('/home/ariana/SISAL/expert/experts/entities_with_SISAL_chronology.xlsx')
colnames(sisal.ext) <- c('entity_id', 'entity_name', 'age_model')

exp.ranking.final.new.2 <- exp.ranking.final.new %>% replace(. == 0, NA_real_) %>% 
  mutate_at(vars(LR, LI, copRa, StalAge, Bacon, Bchron, OxCal, DROP), list(~replace(.,which(.<2), NA_real_))) %>% rename(entity_id = EntityID)

exp.chrono <- left_join(full_chrono_new_new_new, exp.ranking.final.new.2, by = 'entity_id') %>% #group_by(entity_id) %>%
  mutate(lin_reg_age = if_else(is.na(LR) & is.na(DROP), lin_reg_age, NA_real_),
         lin_reg_age_uncert_pos = if_else(is.na(LR) & is.na(DROP),lin_reg_age_uncert_pos, NA_real_),
         lin_reg_age_uncert_neg = if_else(is.na(LR) & is.na(DROP),lin_reg_age_uncert_neg, NA_real_),
         lin_interp_age = if_else(is.na(LI) & is.na(DROP),lin_interp_age, NA_real_),
         lin_interp_age_uncert_pos = if_else(is.na(LI) & is.na(DROP),lin_interp_age_uncert_pos, NA_real_),
         lin_interp_age_uncert_neg = if_else(is.na(LI) & is.na(DROP),lin_interp_age_uncert_neg, NA_real_),
         copRa_age = if_else(is.na(copRa) & is.na(DROP),copRa_age, NA_real_),
         copRa_age_uncert_pos = if_else(is.na(copRa) & is.na(DROP),copRa_age_uncert_pos, NA_real_),
         copRa_age_uncert_neg = if_else(is.na(copRa) & is.na(DROP),copRa_age_uncert_neg, NA_real_),
         StalAge_age = if_else(is.na(StalAge) & is.na(DROP),StalAge_age, NA_real_),
         StalAge_age_uncert_pos = if_else(is.na(StalAge) & is.na(DROP), StalAge_age_uncert_pos, NA_real_),
         StalAge_age_uncert_neg = if_else(is.na(StalAge) & is.na(DROP),StalAge_age_uncert_neg, NA_real_),
         bacon_age = if_else(is.na(Bacon) & is.na(DROP), bacon_age, NA_real_),
         bacon_age_uncert_pos = if_else(is.na(Bacon) & is.na(DROP),bacon_age_uncert_pos, NA_real_),
         bacon_age_uncert_neg = if_else(is.na(Bacon) & is.na(DROP), bacon_age_uncert_neg, NA_real_),
         bchron_age = if_else(is.na(Bchron) & is.na(DROP),bchron_age, NA_real_),
         bchron_age_uncert_pos = if_else(is.na(Bchron) & is.na(DROP),bchron_age_uncert_pos, NA_real_),
         bchron_age_uncert_neg = if_else(is.na(Bchron) & is.na(DROP),bchron_age_uncert_neg, NA_real_),
         OxCal_age = if_else(is.na(DROP) & !(entity_id %in% missing$entity_id), OxCal_age, NA_real_),
         OxCal_age_uncert_pos = if_else(is.na(DROP) & !(entity_id %in% missing$entity_id),OxCal_age_uncert_pos, NA_real_),
         OxCal_age_uncert_neg = if_else(is.na(DROP) & !(entity_id %in% missing$entity_id),OxCal_age_uncert_neg, NA_real_)) %>% 
  filter(!(entity_id %in% missing$entity_id)) %>% 
  select(entity_id, sample_id, depth_sample,
         lin_reg_age, lin_reg_age_uncert_pos, lin_reg_age_uncert_neg,
         lin_interp_age, lin_interp_age_uncert_pos, lin_interp_age_uncert_neg,
         copRa_age, copRa_age_uncert_pos, copRa_age_uncert_neg, 
         StalAge_age, StalAge_age_uncert_pos, StalAge_age_uncert_neg, 
         bacon_age, bacon_age_uncert_pos, bacon_age_uncert_neg, 
         bchron_age, bchron_age_uncert_pos, bchron_age_uncert_neg, 
         OxCal_age, OxCal_age_uncert_pos, OxCal_age_uncert_neg)

#exp.chrono <- exp.chrono %>% filter(!(sample_id %in% sisalv2$hiatus$sample_id))

date_file <- sisalv2$dating %>% rename (date_used_lin_reg = date_used_linear_regress, date_used_copRa = date_used_COPRA, date_used_lin_interp = date_used_linear) %>% mutate(date_used_lin_reg = date_used,
                                       date_used_lin_interp = date_used,
                                       date_used_copRa = date_used,
                                       date_used_StalAge = date_used,
                                       date_used_Bacon = date_used,
                                       date_used_Bchron = date_used) %>%
  select(entity_id, dating_id, starts_with('date_used')) %>% replace(. == 'unknown', 'no') %>% left_join(., exp.ranking.final.new.2, by = 'entity_id') %>%
  mutate(date_used_lin_reg = if_else(is.na(LR) & is.na(DROP) & entity_id %in% SISAL_eval$entity_id, date_used, 'NULL'),
         date_used_lin_interp = if_else(is.na(LI) & is.na(DROP) & entity_id %in% SISAL_eval$entity_id, date_used, 'NULL'),
         date_used_copRa = if_else(is.na(copRa) & is.na(DROP) & entity_id %in% SISAL_eval$entity_id, date_used, 'NULL'),
         date_used_StalAge = if_else(is.na(StalAge) & is.na(DROP) & entity_id %in% SISAL_eval$entity_id, date_used, 'NULL'),
         date_used_Bacon = if_else(is.na(Bacon) & is.na(DROP) & entity_id %in% SISAL_eval$entity_id, date_used, 'NULL'),
         date_used_Bchron = if_else(is.na(Bchron) & is.na(DROP) & entity_id %in% SISAL_eval$entity_id, date_used, 'NULL')) %>% left_join(., oxcal_dating, by = c('entity_id', 'dating_id')) %>% 
  select(entity_id, dating_id, starts_with('date_used')) %>% replace(is.na(.), 'NULL') %>% left_join(., SISAL_eval, by = 'entity_id') %>%
  mutate(date_used_lin_reg = if_else(lR > 0, date_used_lin_reg, 'NULL'),
         date_used_lin_interp = if_else(lI > 0, date_used_lin_interp, 'NULL'),
         date_used_copRa = if_else(copRa > 0, date_used_copRa, 'NULL'),
         date_used_StalAge = if_else(stalage > 0, date_used_StalAge, 'NULL'),
         date_used_Bacon = if_else(bacon > 0, date_used_Bacon, 'NULL'),
         date_used_Bchron = if_else(bchron > 0, date_used_Bchron, 'NULL'),
         date_used_OxCal = if_else(oxcal > 0, date_used_OxCal, 'NULL')) %>%
  select(entity_id, dating_id, starts_with('date_used')) %>% select(-date_used)

write.csv(date_file,'~/SISAL/v2_upd_2/date_file.csv', row.names = F)


final_eval <- exp.chrono %>% ungroup() %>% group_by(entity_id) %>% 
  mutate(., diff_lr = lead(lin_reg_age) - lin_reg_age,
         diff_li = lead(lin_interp_age) - lin_interp_age,
         diff_copRa = lead(copRa_age) - copRa_age,
         diff_bacon = lead(bacon_age) - bacon_age,
         diff_bchron = lead(bchron_age) - bchron_age,
         diff_stalage = lead(StalAge_age) - StalAge_age,
         diff_oxcal = lead(OxCal_age) - OxCal_age) %>%
  mutate(rev_lr = if_else( diff_lr< 0, FALSE, TRUE),
         rev_li = if_else( diff_li< 0, FALSE, TRUE),
         rev_copRa = if_else( diff_copRa< 0, FALSE, TRUE),
         rev_bacon = if_else( diff_bacon< 0, FALSE, TRUE),
         rev_bchron = if_else( diff_bchron< 0, FALSE, TRUE),
         rev_StalAge = if_else( diff_stalage< 0, FALSE, TRUE),
         rev_oxcal = if_else( diff_oxcal< 0, FALSE, TRUE)) %>%
  group_by(entity_id) %>% 
  summarise(lR_reversal = if_else(all(is.na(rev_lr)),-1, if_else(any(!(rev_lr), na.rm = T), 0, 1)),
            lI_reversal = if_else(all(is.na(rev_li)),-1, if_else(any(!(rev_li), na.rm = T), 0, 1)),
            copRa_reversal = if_else(all(is.na(rev_copRa)),-1, if_else(any(!(rev_copRa), na.rm = T), 0, 1)),
            bacon_reversal = if_else(all(is.na(rev_bacon)),-1, if_else(any(!(rev_bacon), na.rm = T), 0, 1)),
            bchron_reversal = if_else(all(is.na(rev_bchron)),-1, if_else(any(!(rev_bchron), na.rm = T), 0, 1)),
            stalage_reversal = if_else(all(is.na(rev_StalAge)),-1, if_else(any(!(rev_StalAge), na.rm = T), 0, 1)),
            oxcal_reversal = if_else(all(is.na(rev_oxcal)),-1, if_else(any(!(rev_oxcal), na.rm = T), 0, 1)))

t <- exp.chrono %>% filter(entity_id == 186)

exp.chrono.final <- exp.chrono %>% ungroup() %>% mutate_at(vars(everything()), list(~replace(.,which(is.na(.)), 'NULL'))) %>% filter(!(sample_id %in% hiatus$sample_id))
  
exp.chrono.fil <- exp.chrono.final %>% filter(!(entity_id %in% sisal.ext$entity_id)) 

exp.chrono.sisal <- exp.chrono.final %>% filter(entity_id %in% sisal.ext$entity_id)
write.csv(exp.chrono.fil, '~/SISAL/v2_upd_2/sisal_chrono_final.csv', row.names = F)
write.csv(exp.chrono.sisal, '~/SISAL/v2_upd_2/sisal_chrono_exists.csv', row.names = F)

correct_am(exp.chrono, '/home/ariana/SISAL/v2_upd/', '/home/ariana/SISAL/v2_upd_2/', entity)

correct_am <- function(exp.chrono, file_path, file_path_new, entity.) {
  eID <- exp.chrono %>% distinct(entity_id)
  
  for(i in eID$entity_id) {
    print(i)
    entity_name <- (entity. %>% filter(entity_id == i))$entity_name
    file_name <- get(load(file.path(file_path, paste(i, '-', entity_name, '.RData', sep = ''))))
    exp.chrono.file <- exp.chrono %>% ungroup() %>% filter(entity_id == i)
    file_name$linReg$AM <- exp.chrono.file %>% select(sample_id, lin_reg_age, lin_reg_age_uncert_pos, lin_reg_age_uncert_neg)
    file_name$linInterp$AM <- exp.chrono.file %>% select(sample_id, lin_interp_age, lin_interp_age_uncert_pos, lin_interp_age_uncert_neg)
    file_name$copRa$AM <- exp.chrono.file %>% select(sample_id, copRa_age, copRa_age_uncert_pos, copRa_age_uncert_neg)
    file_name$StalAge$AM <- exp.chrono.file %>% select(sample_id, StalAge_age, StalAge_age_uncert_pos, StalAge_age_uncert_neg)
    file_name$Bacon$AM <- exp.chrono.file %>% select(sample_id, bacon_age, bacon_age_uncert_pos, bacon_age_uncert_neg)
    file_name$Bchron$AM <- exp.chrono.file %>% select(sample_id, bchron_age, bchron_age_uncert_pos, bchron_age_uncert_neg)
    file_name$OxCal$AM <- exp.chrono.file %>% select(sample_id, OxCal_age, OxCal_age_uncert_pos, OxCal_age_uncert_neg)
    
    if(all(is.na(file_name$linReg$AM$lin_reg_age))){file_name$linReg$AM <- list()}
    if(all(is.na(file_name$linInterp$AM$lin_interp_age))){file_name$linInterp$AM <- list()}
    if(all(is.na(file_name$copRa$AM$copRa_age))){file_name$copRa$AM <- list()}
    if(all(is.na(file_name$StalAge$AM$StalAge_age))){file_name$StalAge$AM <- list()}
    if(all(is.na(file_name$Bacon$AM$bacon_age))){file_name$Bacon$AM <- list()}
    if(all(is.na(file_name$Bchron$AM$bchron_age))){file_name$Bchron$AM <- list()}
    if(all(is.na(file_name$OxCal$AM$OxCal_age))){file_name$OxCal$AM <- list()}
    
    save(file_name, file = file.path(file_path_new, paste(i, '-', entity_name, '.RData', sep = '')))
  }
}

correct_gr(exp.chrono, entity, '/home/ariana/SISAL/v2_upd_2/')

correct_gr <- function(exp, entity., file_path){
  eID <- exp %>% distinct(entity_id)
  
  for(i in eID$entity_id) {
    print(i)
    entity_name <- (entity. %>% filter(entity_id == i))$entity_name
    file_name <- get(load(file.path(file_path, paste(i, '-', entity_name, '.RData', sep = ''))))
    exp.chrono.file <- exp %>% ungroup() %>% filter(entity_id == i)
    
    if(!plyr::empty(data.frame(file_name$linReg$gr))){
      if(plyr::empty(data.frame(file_name$linReg$AM))){
        file_name$linReg$gr <- list()
      } else{
        file_name$linReg$gr <- cbind(exp.chrono.file[-dim(exp.chrono.file)[1],], file_name$linReg$gr) %>% mutate(lin_reg_gr = if_else(is.na(lin_reg_age), NA_real_, lin_reg_gr),
                                                                                                                 lin_reg_gr_uncert_pos = if_else(is.na(lin_reg_age_uncert_pos), NA_real_, lin_reg_gr_uncert_pos),
                                                                                                                 lin_reg_gr_uncert_neg = if_else(is.na(lin_reg_age_uncert_neg), NA_real_, lin_reg_gr_uncert_neg)) %>%
          select(lin_reg_gr, lin_reg_gr_uncert_pos, lin_reg_gr_uncert_neg)
      }
      
    }
    
    if(!plyr::empty(data.frame(file_name$linInterp$gr))){
      if(plyr::empty(data.frame(file_name$linInterp$AM))){
        file_name$linInterp$gr <- list()
      } else{
        file_name$linInterp$gr <- cbind(exp.chrono.file[-dim(exp.chrono.file)[1],], file_name$linInterp$gr) %>% mutate(lin_interp_gr = if_else(is.na(lin_interp_age), NA_real_, lin_interp_gr),
                                                                                                                    lin_interp_gr_uncert_pos = if_else(is.na(lin_interp_age_uncert_pos), NA_real_, lin_interp_gr_uncert_pos),
                                                                                                                    lin_interp_gr_uncert_neg = if_else(is.na(lin_interp_age_uncert_neg), NA_real_, lin_interp_gr_uncert_neg)) %>%
          select(lin_interp_gr, lin_interp_gr_uncert_pos, lin_interp_gr_uncert_neg)
      }
    }
    
    if(!plyr::empty(data.frame(file_name$copRa$gr))){
      if(plyr::empty(data.frame(file_name$copRa$AM))){
        file_name$copRa$gr <- list()
      } else{
        file_name$copRa$gr <- cbind(exp.chrono.file[-dim(exp.chrono.file)[1],], file_name$copRa$gr) %>% mutate(copRa_gr = if_else(is.na(copRa_age), NA_real_, copRa_gr),
                                                                                                          copRa_gr_uncert_pos = if_else(is.na(copRa_age_uncert_pos), NA_real_, copRa_gr_uncert_pos),
                                                                                                          copRa_gr_uncert_neg = if_else(is.na(copRa_age_uncert_neg), NA_real_, copRa_gr_uncert_neg)) %>%
          select(copRa_gr, copRa_gr_uncert_pos, copRa_gr_uncert_neg)
      }
    }
    
    if(!plyr::empty(data.frame(file_name$Bacon$gr))){
      if(plyr::empty(data.frame(file_name$Bacon$AM))){
        file_name$Bacon$gr <- list()
      } else{
        file_name$Bacon$gr <- cbind(exp.chrono.file[-dim(exp.chrono.file)[1],], file_name$Bacon$gr) %>% mutate(bacon_gr = if_else(is.na(bacon_age), NA_real_, bacon_gr),
                                                                                                          bacon_gr_uncert_pos = if_else(is.na(bacon_age_uncert_pos), NA_real_, bacon_gr_uncert_pos),
                                                                                                          bacon_gr_uncert_neg = if_else(is.na(bacon_age_uncert_neg), NA_real_, bacon_gr_uncert_neg)) %>%
          select(bacon_gr, bacon_gr_uncert_pos, bacon_gr_uncert_neg)
      }
    }
    
    if(!plyr::empty(data.frame(file_name$Bchron$gr))){
      if(plyr::empty(data.frame(file_name$Bchron$AM))){
        file_name$Bchron$gr <- list()
      } else{
        file_name$Bchron$gr <- cbind(exp.chrono.file[-dim(exp.chrono.file)[1],], file_name$Bchron$gr) %>% mutate(bchron_gr = if_else(is.na(bchron_age), NA_real_, bchron_gr),
                                                                                                              bchron_gr_uncert_pos = if_else(is.na(bchron_age_uncert_pos), NA_real_, bchron_gr_uncert_pos),
                                                                                                              bchron_gr_uncert_neg = if_else(is.na(bchron_age_uncert_neg), NA_real_, bchron_gr_uncert_neg)) %>%
          select(bchron_gr, bchron_gr_uncert_pos, bchron_gr_uncert_neg)
      }
    }

    save(file_name, file = file.path(file_path, paste(i, '-', entity_name, '.RData', sep = '')))
  }
}

correct_ens(sc_final, sisalv2$entity, '/home/ariana/SISAL/v2_upd_2/')

correct_ens <- function(exp, entity., file_path){
  eID <- exp %>% distinct(entity_id)
  
  for(i in eID$entity_id) {
    print(i)
    entity_name <- (entity. %>% filter(entity_id == i))$entity_name
    file_name <- get(load(file.path(file_path, paste(i, '-', entity_name, '.RData', sep = ''))))
    #exp.chrono.file <- exp %>% ungroup() %>% filter(entity_id == i)
    
    if(!plyr::empty(data.frame(file_name$linReg$ens))){
      if(plyr::empty(data.frame(file_name$linReg$AM))){
        file_name$linReg$ens <- list()
      } else{
        file_name$linReg$ens[which(is.na(file_name$linReg$AM$lin_reg_age)),] <- NA_real_
      }
      
    }
    
    if(!plyr::empty(data.frame(file_name$linInterp$gr))){
      if(plyr::empty(data.frame(file_name$linInterp$AM))){
        file_name$linInterp$ens <- list()
      } else{
        file_name$linInterp$ens[which(is.na(file_name$linInterp$AM$lin_interp_age)),] <- NA_real_
      }
    }
    
    if(!plyr::empty(data.frame(file_name$copRa$gr))){
      if(plyr::empty(data.frame(file_name$copRa$AM))){
        file_name$copRa$ens <- list()
      } else{
        file_name$copRa$ens[which(is.na(file_name$copRa$AM$copRa_age)),] <- NA_real_
      }
    }
    
    if(!plyr::empty(data.frame(file_name$Bacon$gr))){
      if(plyr::empty(data.frame(file_name$Bacon$AM))){
        file_name$Bacon$ens <- list()
      } else{
        file_name$Bacon$ens[which(is.na(file_name$Bacon$AM$bacon_age)),] <- NA_real_
      }
    }
    
    if(!plyr::empty(data.frame(file_name$Bchron$gr))){
      if(plyr::empty(data.frame(file_name$Bchron$AM))){
        file_name$Bchron$ens <- list()
      } else{
        file_name$Bchron$ens[which(is.na(file_name$Bchron$AM$bchron_age)),] <- NA_real_
      }
    }
    
    save(file_name, file = file.path(file_path, paste(i, '-', entity_name, '.RData', sep = '')))
  }
}

drop <- exp.ranking.final.new.2 %>% select(entity_id, DROP) 
drop <- left_join(r_f %>% select(entity_id), drop)

final.plot(exp.chrono, sisalv2$entity, sisalv2$dating, drop)


date_file <- read.csv('~/SISAL/date_file.csv', header = T, stringsAsFactors = F) %>% filter(entity_id %in% su$entity_id)


write.csv(date_file, '~/SISAL/date_file.csv', row.names = F)

sisal.plot<- function(s){
  eId <- (s$entity %>% distinct(entity_id) %>% arrange(entity_id))$entity_id
  
  for(i in eId) {
    print(i)
    
    entity_name <- (.entity %>% filter(entity_id == i))$entity_name
    graphics.off()
    df_fil <- get(load(paste('/home/ariana/SISAL/v2_upd_2/',i,'-',entity_name, '.RData' ,sep = '')))
    
    cairo_pdf(paste("~/SISAL/v2_upd_2/",i,'-',entity_name,'.pdf',sep = ''), 12, 12)
    layout(matrix(c(1,1,2,3,3,4,5,5,6,7,7,7), nrow=4, ncol=3, byrow = T))
    par(mar = c(4,3,3,4), oma = c(1,3,2,1))
    x.lim <- plot_am(df_fil, i, .dating, entity_name)
    if(any(c(!is_empty(data.frame(df_fil$Bacon$AM)),!is_empty(data.frame(df_fil$Bchron$AM)),!is_empty(data.frame(df_fil$copRa$AM)),!is_empty(data.frame(df_fil$StalAge$AM)),!is_empty(data.frame(df_fil$linInterp$AM)),
             !is_empty(data.frame(df_fil$linReg$AM)),!is_empty(data.frame(df_fil$OxCal$AM)))) & is.na(DROP)){
      #if(any(c(m$Bacon,m$Bchron, m$StalAge, m$linInterp, m$copRa, m$linReg, m$OxCal)) & is.na(DROP)){
      plot_iqr(df_fil, i, .dating, entity_name)
      mtext(paste(i,'-', entity_name, sep = ''), side = 3, line = -2, outer = TRUE)
      if(all(c(is_empty(data.frame(df_fil$Bacon$gr)),is_empty(data.frame(df_fil$Bchron$gr)),is_empty(data.frame(df_fil$copRa$gr)),is_empty(data.frame(df_fil$linInterp$gr)),
               is_empty(data.frame(df_fil$linReg$gr))))){
        plot.new()
        plot.new()
      } else {
        plot_gr(df_fil,i, .dating, entity_name, x.lim)
        plot_gr_iqr(df_fil)
      }
    } else {
      plot.new()
      mtext(paste(i,'-', entity_name, sep = ''), side = 3, line = -2, outer = TRUE)
      plot.new()
      plot.new()
    }
    plot_isotopes(df_fil, i, entity_name, sisal.chrono, x.lim)
    plot_dating(df_fil, i, .dating)
    plot.new()
    
    orig_am <- df_fil$origAM %>% filter(age_model_type != 'NA') %>% distinct(age_model_type)
    
    legend('center',legend=c(paste('original AM:', orig_am),'copRa median age','StalAge mean age','lin. interp. median age','lin. reg. median age','Bacon median age','Bchron median age', 'OxCal median age', "Event; actively forming", "Event; hiatus"),
           lty=c(1,1,1,1,1,1,1,1,3,3),cex=1.5,
           col = c('black','royalblue3','forestgreen','skyblue','springgreen2','purple', 'sienna1', 'brown2',"orange","grey"),title=paste("Age model Information"),bty = "n",  ncol = 3, lwd = 4)
    
    dev.off()
    
    
  }
  
}








########## old stuff -----------------------####

correct_am <- function(eID, am, rm_hiatus, rm_dt, dating, hiatus){
  orig_names <- colnames(am)
  colnames(am) <- c('sample_id','age','uncert_pos','uncert_neg')
  
  dt <- dating_from_base %>% filter(entity_id == eID) %>% filter(date_type != 'Event; hiatus') %>% mutate(dating_id = -dating_id) %>% filter(date_used == 'yes') %>% select(entity_id, dating_id, depth_dating) %>% rename(depth_sample = depth_dating, sample_id = dating_id)
  dt <- dt %>% ungroup() %>%
    mutate(age= rep(NA_real_, length(dt$sample_id)),
           uncert_pos = rep(NA_real_, length(dt$sample_id)),
           uncert_neg = rep(NA_real_, length(dt$sample_id))) %>%
    mutate_at(vars(everything()), as.character) %>%
    mutate_at(vars(everything()), as.numeric)
  
  am_new <- am %>% mutate(remove = if_else(sample_id %in% hiatus$sample_id,1,0)) %>% 
    mutate(remove_new = cumsum(remove)) %>% 
    mutate(age = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,age),
           uncert_pos = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,uncert_pos),
           uncert_neg = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,uncert_neg))
  
  am_new_new <- bind_rows(am_new, dt) %>% select(-remove_new) %>% arrange(depth_sample, .by_group = T) %>%
    mutate(remove = if_else(sample_id %in% rm_dt$dating_id, 1, 
                            if_else(sample_id %in% hiatus$sample_id,2,
                                    if_else(sample_id %in% help$dating_id,-3,0))),
           r = cumsum(remove)) %>%
    mutate(age = if_else(r%%2 != 0, NA_real_,age),
           uncert_pos = if_else(r%%2 != 0, NA_real_,uncert_pos),
           uncert_neg = if_else(r%%2 != 0, NA_real_,uncert_neg)) %>%
    filter(!(sample_id %in% dt$sample_id)) %>% select(-c(remove, r)) %>% filter(!(sample_id %in% hiatus$sample_id))
  
  colnames(am) <- orig_names
  
  return(am_new_new)
  
}


s_fil_new_new <- s_fil_new %>% mutate(remove = if_else(sample_id %in% hiatus$sample_id,1,0)) %>% group_by(entity_id) %>%
  mutate(remove_new = cumsum(remove)) %>% 
  mutate(lin_reg_age = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,lin_reg_age),
         lin_reg_age_uncert_pos = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,lin_reg_age_uncert_pos),
         lin_reg_age_uncert_neg = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,lin_reg_age_uncert_neg),
         lin_interp_age = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,lin_interp_age),
         lin_interp_age_uncert_pos = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,lin_interp_age_uncert_pos),
         lin_interp_age_uncert_neg = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,lin_interp_age_uncert_neg),
         copRa_age = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,copRa_age),
         copRa_age_uncert_pos = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,copRa_age_uncert_pos),
         copRa_age_uncert_neg = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,copRa_age_uncert_neg),
         StalAge_age = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,StalAge_age),
         StalAge_age_uncert_pos = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,StalAge_age_uncert_pos),
         StalAge_age_uncert_neg = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,StalAge_age_uncert_neg),
         bacon_age = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,bacon_age),
         bacon_age_uncert_pos = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,bacon_age_uncert_pos),
         bacon_age_uncert_neg = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,bacon_age_uncert_neg),
         bchron_age = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,bchron_age),
         bchron_age_uncert_pos = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,bchron_age_uncert_pos),
         bchron_age_uncert_neg = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,bchron_age_uncert_neg))

s_fil_new_new_new <- bind_rows(s_fil_new_new, dt) %>% select(-remove_new) %>%group_by(entity_id) %>% arrange(depth_sample, .by_group = T) %>%
  mutate(remove = if_else(sample_id %in% rm_dt$dating_id, 1, 
                          if_else(sample_id %in% hiatus$sample_id,2,
                                  if_else(sample_id %in% help$dating_id,-3,0))),
         r = cumsum(remove)) %>%
  mutate(lin_reg_age = if_else(r%%2 != 0, NA_real_,lin_reg_age),
         lin_reg_age_uncert_pos = if_else(r%%2 != 0, NA_real_,lin_reg_age_uncert_pos),
         lin_reg_age_uncert_neg = if_else(r%%2 != 0, NA_real_,lin_reg_age_uncert_neg),
         lin_interp_age = if_else(r%%2 != 0, NA_real_,lin_interp_age),
         lin_interp_age_uncert_pos = if_else(r%%2 != 0, NA_real_,lin_interp_age_uncert_pos),
         lin_interp_age_uncert_neg = if_else(r%%2 != 0, NA_real_,lin_interp_age_uncert_neg),
         copRa_age = if_else(r%%2 != 0, NA_real_,copRa_age),
         copRa_age_uncert_pos = if_else(r%%2 != 0, NA_real_,copRa_age_uncert_pos),
         copRa_age_uncert_neg = if_else(r%%2 != 0, NA_real_,copRa_age_uncert_neg),
         StalAge_age = if_else(r%%2 != 0, NA_real_,StalAge_age),
         StalAge_age_uncert_pos = if_else(r%%2 != 0, NA_real_,StalAge_age_uncert_pos),
         StalAge_age_uncert_neg = if_else(r%%2 != 0, NA_real_,StalAge_age_uncert_neg),
         bacon_age = if_else(r%%2 != 0, NA_real_,bacon_age),
         bacon_age_uncert_pos = if_else(r%%2 != 0, NA_real_,bacon_age_uncert_pos),
         bacon_age_uncert_neg = if_else(r%%2 != 0, NA_real_,bacon_age_uncert_neg),
         bchron_age = if_else(r%%2 != 0, NA_real_,bchron_age),
         bchron_age_uncert_pos = if_else(r%%2 != 0, NA_real_,bchron_age_uncert_pos),
         bchron_age_uncert_neg = if_else(r%%2 != 0, NA_real_,bchron_age_uncert_neg)) %>%
  filter(!(sample_id %in% dt$sample_id)) %>% select(-c(remove, r)) %>% filter(!(sample_id %in% hiatus$sample_id))



write.csv(s_fil_new_new_new, '~/SISAL/sisal_chronology_iup.csv', row.names = F)



full_chrono_new <- full_chrono_eval %>% mutate(remove = if_else(sample_id %in% hiatus$sample_id,1,0)) %>% group_by(entity_id) %>%
  mutate(remove_new = cumsum(remove)) %>% 
  mutate(lin_reg_age = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,lin_reg_age),
         lin_reg_age_uncert_pos = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,lin_reg_age_uncert_pos),
         lin_reg_age_uncert_neg = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,lin_reg_age_uncert_neg),
         lin_interp_age = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,lin_interp_age),
         lin_interp_age_uncert_pos = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,lin_interp_age_uncert_pos),
         lin_interp_age_uncert_neg = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,lin_interp_age_uncert_neg),
         copRa_age = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,copRa_age),
         copRa_age_uncert_pos = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,copRa_age_uncert_pos),
         copRa_age_uncert_neg = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,copRa_age_uncert_neg),
         StalAge_age = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,StalAge_age),
         StalAge_age_uncert_pos = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,StalAge_age_uncert_pos),
         StalAge_age_uncert_neg = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,StalAge_age_uncert_neg),
         bacon_age = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,bacon_age),
         bacon_age_uncert_pos = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,bacon_age_uncert_pos),
         bacon_age_uncert_neg = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,bacon_age_uncert_neg),
         bchron_age = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,bchron_age),
         bchron_age_uncert_pos = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,bchron_age_uncert_pos),
         bchron_age_uncert_neg = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,bchron_age_uncert_neg),
         OxCal_age = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,OxCal_age),
         OxCal_age_uncert_pos = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,OxCal_age_uncert_pos),
         OxCal_age_uncert_neg = if_else(remove_new < 1 & entity_id %in% rm_hiatus$entity_id, NA_real_,OxCal_age_uncert_neg))


full_chrono_new <- bind_rows(full_chrono_eval,dt_h) %>% group_by(entity_id) %>% arrange(depth_sample, .by_group = T) %>% filter(!(sample_id %in% hiatus$sample_id)) %>%
  mutate(remove = if_else(is.na(remove),0,remove)) %>%
  mutate(remove_new = cumsum(remove)) %>% 
  mutate(lin_reg_age = if_else(remove_new < 3  & entity_id %in% c(top$entity_id, bottom$entity_id) , NA_real_,lin_reg_age),
lin_reg_age_uncert_pos = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,lin_reg_age_uncert_pos),
lin_reg_age_uncert_neg = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,lin_reg_age_uncert_neg),
lin_interp_age = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,lin_interp_age),
lin_interp_age_uncert_pos = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,lin_interp_age_uncert_pos),
lin_interp_age_uncert_neg = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,lin_interp_age_uncert_neg),
copRa_age = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,copRa_age),
copRa_age_uncert_pos = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,copRa_age_uncert_pos),
copRa_age_uncert_neg = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,copRa_age_uncert_neg),
StalAge_age = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,StalAge_age),
StalAge_age_uncert_pos = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,StalAge_age_uncert_pos),
StalAge_age_uncert_neg = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,StalAge_age_uncert_neg),
bacon_age = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,bacon_age),
bacon_age_uncert_pos = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,bacon_age_uncert_pos),
bacon_age_uncert_neg = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,bacon_age_uncert_neg),
bchron_age = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,bchron_age),
bchron_age_uncert_pos = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,bchron_age_uncert_pos),
bchron_age_uncert_neg = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,bchron_age_uncert_neg),
OxCal_age = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,OxCal_age),
OxCal_age_uncert_pos = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,OxCal_age_uncert_pos),
OxCal_age_uncert_neg = if_else(remove_new < 3 & entity_id %in% c(top$entity_id, bottom$entity_id), NA_real_,OxCal_age_uncert_neg)) %>% 
  filter(!(sample_id %in% dt_h$sample_id))

exp.chrono <- left_join(full_chrono_new_new_new, exp.lk, by = 'entity_id') %>% #group_by(entity_id) %>%
  mutate(lin_reg_age = if_else(is.na(LR) & is.na(DROP), lin_reg_age, NA_real_),
         lin_reg_age_uncert_pos = if_else(is.na(LR) & is.na(DROP),lin_reg_age_uncert_pos, NA_real_),
         lin_reg_age_uncert_neg = if_else(is.na(LR) & is.na(DROP),lin_reg_age_uncert_neg, NA_real_),
         lin_interp_age = if_else(is.na(LI) & is.na(DROP),lin_interp_age, NA_real_),
         lin_interp_age_uncert_pos = if_else(is.na(LI) & is.na(DROP),lin_interp_age_uncert_pos, NA_real_),
         lin_interp_age_uncert_neg = if_else(is.na(LI) & is.na(DROP),lin_interp_age_uncert_neg, NA_real_),
         copRa_age = if_else(is.na(copRa) & is.na(DROP),copRa_age, NA_real_),
         copRa_age_uncert_pos = if_else(is.na(copRa) & is.na(DROP),copRa_age_uncert_pos, NA_real_),
         copRa_age_uncert_neg = if_else(is.na(copRa) & is.na(DROP),copRa_age_uncert_neg, NA_real_),
         StalAge_age = if_else(is.na(StalAge) & is.na(DROP),StalAge_age, NA_real_),
         StalAge_age_uncert_pos = if_else(is.na(StalAge) & is.na(DROP), StalAge_age_uncert_pos, NA_real_),
         StalAge_age_uncert_neg = if_else(is.na(StalAge) & is.na(DROP),StalAge_age_uncert_neg, NA_real_),
         bacon_age = if_else(is.na(Bacon) & is.na(DROP),bacon_age, NA_real_),
         bacon_age_uncert_pos = if_else(is.na(Bacon) & is.na(DROP),bacon_age_uncert_pos, NA_real_),
         bacon_age_uncert_neg = if_else(is.na(Bacon) & is.na(DROP),bacon_age_uncert_neg, NA_real_),
         bchron_age = if_else(is.na(Bchron) & is.na(DROP),bchron_age, NA_real_),
         bchron_age_uncert_pos = if_else(is.na(Bchron) & is.na(DROP),bchron_age_uncert_pos, NA_real_),
         bchron_age_uncert_neg = if_else(is.na(Bchron) & is.na(DROP),bchron_age_uncert_neg, NA_real_),
         OxCal_age = if_else(is.na(DROP) & !(entity_id %in% missing$entity_id), OxCal_age, NA_real_),
         OxCal_age_uncert_pos = if_else(is.na(DROP) & !(entity_id %in% missing$entity_id),OxCal_age_uncert_pos, NA_real_),
         OxCal_age_uncert_neg = if_else(is.na(DROP) & !(entity_id %in% missing$entity_id),OxCal_age_uncert_neg, NA_real_)) %>% 
  filter(!(entity_id %in% missing$entity_id)) %>% 
  select(entity_id, sample_id, depth_sample,
         lin_reg_age, lin_reg_age_uncert_pos, lin_reg_age_uncert_neg,
         lin_interp_age, lin_interp_age_uncert_pos, lin_interp_age_uncert_neg,
         copRa_age, copRa_age_uncert_pos, copRa_age_uncert_neg, 
         StalAge_age, StalAge_age_uncert_pos, StalAge_age_uncert_neg, 
         bacon_age, bacon_age_uncert_pos, bacon_age_uncert_neg, 
         bchron_age, bchron_age_uncert_pos, bchron_age_uncert_neg, 
         OxCal_age, OxCal_age_uncert_pos, OxCal_age_uncert_neg)

