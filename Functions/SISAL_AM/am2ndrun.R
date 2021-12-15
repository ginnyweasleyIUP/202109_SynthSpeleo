run_interim<- tibble(entity_id = c(277,367,369,525,548,561,582,614))
runFile_interim <- tibble(entity_id = run$entity_id, Bacon = rep(T, length(run$entity_id)), Bchron = rep(T, length(run$entity_id)), 
                  StalAge = rep(T, length(run$entity_id)), linInterp = rep(T, length(run$entity_id)), copRa = rep(T, length(run$entity_id)),
                  linReg = rep(T, length(run$entity_id)), 
                  working_directory = rep("/stacywork/ariana/v2_new/", length(run$entity_id)), 
                  wd = rep('/home/ariana/SISAL Data/sisalv2', length(run$entity_id))) 

sisalv2 <- read.SISAL.files('~/SISAL Data/sisalv2/', '')

r_interim <- merge_SISAL_chrono_final(runFile_interim, left_join(sisalv2$sample, sisalv2$entity) ,sisalv2$entity)
r_2 <- r_interim[[1]]
s_2 <- r_interim[[2]]
s_2_full <- full_join(s_2, oxcal_final) %>% filter(entity_id %in% run_interim$entity_id) %>% select(-depth)

eval_2 <- eval(s_2_full, sisalv2$dating, sisalv2$hiatus)

SISAL_eval_2 <- eval_2[[5]]

s_2_eval <- s_2_full %>% mutate(lin_reg_age = if_else(entity_id == 277, NA_real_, lin_reg_age),
                              lin_reg_age_uncert_pos = if_else(entity_id ==277, NA_real_, lin_reg_age_uncert_pos),
                              lin_reg_age_uncert_neg = if_else(entity_id == 277, NA_real_, lin_reg_age_uncert_neg))


interim_rank <- tibble(entity_id = c(277,369,525,548,582,614),
                       LI =  c(2,2,0,2,0,0),
                       LR = c(2,2,2,2,2,2),
                       copRa = c(2,2,0,2,0,2),
                       StalAge = c(2,2,0,2,0,2),
                       Bacon = c(2,2,0,0,0,2),
                       Bchron = c(0,0,0,0,0,0),
                       OxCal = c(0,0,0,0,0,0),
                       DROP = c(0,0,0,0,0,0)) %>% replace(. == 0, NA_real_)

dt_h_interim <- dt_h %>% filter(entity_id %in% run_interim$entity_id)

s_2_new <- bind_rows(s_2_eval,dt_h) %>% group_by(entity_id) %>% arrange(depth_sample, .by_group = T) %>% filter(!(sample_id %in% sisalv2$hiatus$sample_id)) %>%
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



dt <- sisalv2$dating %>% filter(date_type != 'Event; hiatus') %>% mutate(dating_id = -dating_id) %>% filter(date_used == 'yes') %>% select(entity_id, dating_id, depth_dating) %>% rename(depth_sample = depth_dating, sample_id = dating_id) %>% ungroup() 
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

hiatus_removed <- s_2_new %>% filter(sample_id %in% hiatus$sample_id)

s_2_new_new <- bind_rows(s_2_new, dt) %>% select(-remove_new) %>% group_by(entity_id) %>% arrange(depth_sample, .by_group = T) %>%
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

t <- s_2_new_new %>% group_by(entity_id) %>% slice(.,1)

hi <- left_join(sisalv2$sample, sisalv2$hiatus, by = 'sample_id') %>% filter(entity_id %in% run_interim$entity_id & hiatus == 'H') %>% select(entity_id, sample_id, depth_sample) %>% ungroup() %>%
  mutate(lin_reg_age = rep(NA_real_, 7),
         lin_reg_age_uncert_pos = rep(NA_real_, 7),
         lin_reg_age_uncert_neg = rep(NA_real_, 7),
         lin_interp_age = rep(NA_real_, 7),
         lin_interp_age_uncert_pos = rep(NA_real_, 7),
         lin_interp_age_uncert_neg = rep(NA_real_, 7),
         copRa_age = rep(NA_real_, 7),
         copRa_age_uncert_pos = rep(NA_real_, 7),
         copRa_age_uncert_neg = rep(NA_real_, 7),
         bacon_age = rep(NA_real_, 7),
         bacon_age_uncert_pos = rep(NA_real_, 7),
         bacon_age_uncert_neg = rep(NA_real_, 7),
         bchron_age = rep(NA_real_, 7),
         bchron_age_uncert_pos = rep(NA_real_, 7),
         bchron_age_uncert_neg = rep(NA_real_, 7),
         StalAge_age = rep(NA_real_,7),
         StalAge_age_uncert_pos = rep(NA_real_, 7),
         StalAge_age_uncert_neg = rep(NA_real_, 7),
         OxCal_age = rep(NA_real_, 7),
         OxCal_age_uncert_pos = rep(NA_real_, 7),
         OxCal_age_uncert_neg = rep(NA_real_, 7))%>%
  mutate_at(vars(everything()), as.character) %>%
  mutate_at(vars(everything()), as.numeric)

s_2_new_new_new <- s_2_new_new %>% mutate(lin_reg_age = if_else(lin_reg_age < -68, NA_real_,lin_reg_age),
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
                                                          OxCal_age_uncert_neg = if_else(OxCal_age < -68, NA_real_,OxCal_age_uncert_neg)) %>% filter(entity_id %in% run_interim$entity_id) %>%
  bind_rows(.,hi) %>% group_by(entity_id) %>% arrange(., depth_sample, .by_group = T)

chrono.interim <- left_join(s_2_new_new_new, interim_rank, by = 'entity_id') %>% #group_by(entity_id) %>%
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


te <- s_2_new_new_new %>% group_by(entity_id) %>% slice(.,1)

drop <- tibble(entity_id = run_interim$entity_id, DROP = NA_real_)

merge_SISAL_ensemble_final(r_2, sisalv2$dating)

correct_am(chrono.interim, '/home/ariana/SISAL/v2_upd_2/', '/home/ariana/SISAL/v2_upd_2/', entity)
correct_gr(chrono.interim, sisalv2$entity, '/home/ariana/SISAL/v2_upd_2/')

final.plot(chrono.interim, sisalv2$entity, sisalv2$dating, drop)

write.csv(chrono.interim %>% filter(!(sample_id %in% hi$sample_id)), '~/SISAL/v2_upd_2/sisal_chrono_updated_entites.csv', row.names = F)



date_file_upd <- sisalv2$dating %>% rename (date_used_lin_reg = date_used_linear_regress, date_used_copRa = date_used_COPRA, date_used_lin_interp = date_used_linear) %>% 
  mutate(date_used_lin_reg = date_used,
         date_used_lin_interp = date_used,
         date_used_copRa = date_used,
         date_used_StalAge = date_used,
         date_used_Bacon = date_used,
         date_used_Bchron = date_used) %>%
  select(entity_id, dating_id, starts_with('date_used')) %>% replace(. == 'unknown', 'no') %>% left_join(., interim_rank, by = 'entity_id') %>%
  mutate(date_used_lin_reg = if_else(is.na(LR) & is.na(DROP), date_used, 'NULL'),
         date_used_lin_interp = if_else(is.na(LI) & is.na(DROP), date_used, 'NULL'),
         date_used_copRa = if_else(is.na(copRa) & is.na(DROP), date_used, 'NULL'),
         date_used_StalAge = if_else(is.na(StalAge) & is.na(DROP), date_used, 'NULL'),
         date_used_Bacon = if_else(is.na(Bacon) & is.na(DROP), date_used, 'NULL'),
         date_used_Bchron = if_else(is.na(Bchron) & is.na(DROP), date_used, 'NULL')) %>% left_join(., oxcal_dating, by = c('entity_id', 'dating_id')) %>% 
  select(entity_id, dating_id, starts_with('date_used')) %>% replace(is.na(.), 'NULL') %>% left_join(., SISAL_eval_2, by = 'entity_id') %>%
  mutate(date_used_lin_reg = if_else(lR > 0, date_used_lin_reg, 'NULL'),
         date_used_lin_interp = if_else(lI > 0, date_used_lin_interp, 'NULL'),
         date_used_copRa = if_else(copRa > 0, date_used_copRa, 'NULL'),
         date_used_StalAge = if_else(stalage > 0, date_used_StalAge, 'NULL'),
         date_used_Bacon = if_else(bacon > 0, date_used_Bacon, 'NULL'),
         date_used_Bchron = if_else(bchron > 0, date_used_Bchron, 'NULL'),
         date_used_OxCal = if_else(oxcal > 0, date_used_OxCal, 'NULL'))  %>% 
  select(entity_id, dating_id, starts_with('date_used')) %>% select(-date_used) %>% filter(entity_id %in% run_interim$entity_id)

write.csv(date_file_upd, '~/SISAL/v2_upd_2/date_file_updated_entites.csv', row.names = F)


date_file_final <- left_join(date_file, date_file_upd, by = c('entity_id', 'dating_id')) %>%
  mutate(date_used_lin_reg.x = if_else(entity_id %in% date_file_upd$entity_id, date_used_lin_reg.y, date_used_lin_reg.x),
         date_used_lin_interp.x = if_else(entity_id %in% date_file_upd$entity_id, date_used_lin_interp.y, date_used_lin_interp.x),
         date_used_copRa.x = if_else(entity_id %in% date_file_upd$entity_id, date_used_copRa.y, date_used_copRa.x),
         date_used_StalAge.x = if_else(entity_id %in% date_file_upd$entity_id, date_used_StalAge.y, date_used_StalAge.x),
         date_used_Bacon.x = if_else(entity_id %in% date_file_upd$entity_id, date_used_Bacon.y, date_used_Bacon.x),
         date_used_Bchron.x = if_else(entity_id %in% date_file_upd$entity_id, date_used_Bchron.y, date_used_Bchron.x),
         date_used_OxCal.x = if_else(entity_id %in% date_file_upd$entity_id, date_used_OxCal.y, date_used_OxCal.x)) %>%
  select(entity_id, dating_id, date_used_lin_reg.x, date_used_lin_interp.x, date_used_copRa.x, date_used_StalAge.x, date_used_Bacon.x, date_used_Bchron.x, date_used_OxCal.x) %>%
  rename(date_used_lin_reg = date_used_lin_reg.x, date_used_lin_interp = date_used_lin_interp.x, date_used_copRa = date_used_copRa.x, date_used_StalAge = date_used_StalAge.x, 
         date_used_Bacon = date_used_Bacon.x, date_used_Bchron = date_used_Bchron.x, date_used_OxCal = date_used_OxCal.x) %>% replace(is.na(.),'NULL')
  

write.csv(date_file_final, '~/SISAL/v2_upd_2/date_file_final.csv', row.names = F)














