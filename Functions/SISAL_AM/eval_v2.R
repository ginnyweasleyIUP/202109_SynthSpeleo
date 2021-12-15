######-----------EVALUATION --------------------------#####


setwd("~/SISAL Data/sisalv2/")
sample <- read.csv('sample.csv', header = T, stringsAsFactors = F) %>% mutate_at(vars(depth_sample), as.numeric)
entity <- read.csv('entity.csv', header = T, stringsAsFactors = F)
site <- read.csv('site.csv', header = T, stringsAsFactors = F)
dating <- read.csv('dating.csv', header = T, stringsAsFactors = F) 
dating <- dating %>% mutate_at(vars(dating_id, entity_id, depth_dating, dating_thickness, X14C_correction, corr_age, corr_age_uncert_pos, corr_age_uncert_neg), as.numeric)
orig2 <- read.csv('original_chronology.csv', header = T, stringsAsFactors = F)
site_tb <- left_join(site, entity, by = 'site_id')
sisal2 <- read.csv('sisal_chronology.csv', header = T, stringsAsFactors = F)

entity_from_base <-  entity %>% filter(depth_ref == 'from base') %>% distinct(entity_id)
sample_from_base <- sample %>% filter(entity_id %in% entity_from_base$entity_id) %>% dplyr::select(entity_id,depth_sample) %>% group_by(entity_id) %>% summarise(max = max(depth_sample))
sampling_from_base <- full_join(sample, sample_from_base, by = 'entity_id') %>% group_by(entity_id) %>% 
  mutate(depth_conv = if_else(entity_id %in% entity_from_base$entity_id, max-depth_sample, NA_real_)) %>% 
  mutate(depth_sample = if_else(!is.na(depth_conv), depth_conv, depth_sample)) %>%
  dplyr::select(-depth_conv) %>% arrange(., depth_sample, .by_group = T)

se <- left_join(sampling_from_base, entity, by = 'entity_id')

setwd("~/Documents/v2")
r1 <- read.csv('runFile.csv', header = T) %>% filter(entity_id < 331)
r2 <- read.csv('runFile2.csv', header = T) %>% filter(entity_id < 378 & entity_id >= 331)
r3 <- read.csv('runFile3.csv', header = T) %>% filter(entity_id >=378 & entity_id < 506)
r4 <- read.csv('runFile4.csv', header = T) %>% filter(entity_id >= 507 & entity_id  <510)
r5 <- read.csv('runFile5.csv', header = T) %>% filter(entity_id > 510)
r <- rbind(rbind(rbind(rbind(r1, r2), r3),r4),r5) %>% arrange(., entity_id)

r <- r %>% add_row(., entity_id = 378, Bacon = FALSE, Bchron = FALSE,StalAge = FALSE,linInterp = FALSE,copRa = FALSE,linReg = FALSE,
                   working_directory = '/home/ariana/Documents/v2',wd =  '/home/ariana/SISAL Data/interim_v2_csv_20190812', .before=1)
r <- r %>% add_row(., entity_id = 331, Bacon = FALSE, Bchron = FALSE,StalAge = FALSE,linInterp = FALSE,copRa = FALSE,linReg = FALSE,
                   working_directory = '/home/ariana/Documents/v2',wd =  '/home/ariana/SISAL Data/interim_v2_csv_20190812', .before=1)
r <- r %>% add_row(., entity_id = 506, Bacon = FALSE, Bchron = FALSE,StalAge = FALSE,linInterp = FALSE,copRa = FALSE,linReg = FALSE,
                   working_directory = '/home/ariana/Documents/v2',wd =  '/home/ariana/SISAL Data/interim_v2_csv_20190812', .before=1)
r <- r %>% add_row(., entity_id = 510, Bacon = FALSE, Bchron = FALSE,StalAge = FALSE,linInterp = FALSE,copRa = FALSE,linReg = FALSE,
                   working_directory = '/home/ariana/Documents/v2',wd =  '/home/ariana/SISAL Data/interim_v2_csv_20190812', .before=1)
r <- r %>% arrange(., entity_id) 

SISAL_chronology_new <- data.frame(entity_id = se$entity_id, sample_id = se$sample_id, depth_sample = se$depth_sample,
                                   lin_reg_age = rep(NA_real_, length(sample$sample_id)),
                                   lin_reg_age_uncert_pos = rep(NA_real_, length(sample$sample_id)),
                                   lin_reg_age_uncert_neg = rep(NA_real_, length(sample$sample_id)),
                                   lin_interp_age = rep(NA_real_, length(sample$sample_id)),
                                   lin_interp_age_uncert_pos = rep(NA_real_, length(sample$sample_id)),
                                   lin_interp_age_uncert_neg = rep(NA_real_, length(sample$sample_id)),
                                   bacon_age = rep(NA_real_, length(sample$sample_id)),
                                   bacon_age_uncert_pos = rep(NA_real_, length(sample$sample_id)),
                                   bacon_age_uncert_neg = rep(NA_real_, length(sample$sample_id)),
                                   bchron_age = rep(NA_real_, length(sample$sample_id)),
                                   bchron_age_uncert_pos = rep(NA_real_, length(sample$sample_id)),
                                   bchron_age_uncert_neg = rep(NA_real_, length(sample$sample_id)),
                                   StalAge_age = rep(NA_real_, length(sample$sample_id)),
                                   StalAge_age_uncert_pos = rep(NA_real_, length(sample$sample_id)),
                                   StalAge_age_uncert_neg = rep(NA_real_, length(sample$sample_id)),
                                   copRa_age = rep(NA_real_, length(sample$sample_id)),
                                   copRa_age_uncert_pos = rep(NA_real_, length(sample$sample_id)),
                                   copRa_age_uncert_neg = rep(NA_real_, length(sample$sample_id)))


merge_SISAL_chrono <- function(runFile, SISAL_chronology_new, entity){
  
  for(i in runFile$entity_id) {
    print(i)
    y <- runFile %>% filter(entity_id == i) 
    entity_name <- (entity %>% filter(entity_id == i))$entity_name
    file_name <- paste(i,"-",entity_name, sep = '')
    
    if(y$linReg){
      setwd(file.path(y$working_directory, file_name, '/linReg'))
      lR_chrono <- read.csv('linReg_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) 
      names(lR_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
      SISAL_chronology_new <- full_join(lR_chrono, SISAL_chronology_new, by = 'sample_id') %>% 
        mutate(lin_reg_age = if_else(sample_id %in% lR_chrono$sample_id, age, lin_reg_age),
               lin_reg_age_uncert_pos = if_else(sample_id %in% lR_chrono$sample_id, uncert_pos, lin_reg_age_uncert_pos),
               lin_reg_age_uncert_neg = if_else(sample_id %in% lR_chrono$sample_id, uncert_neg, lin_reg_age_uncert_neg)) %>% 
        select(-age, -uncert_pos, -uncert_neg)
      #print('linreg')
    }
    
    if(y$linInterp){
      setwd(file.path(y$working_directory, file_name, '/linInterp'))
      lI_chrono <- read.csv('linInt_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) 
      names(lI_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
      SISAL_chronology_new <- full_join(lI_chrono, SISAL_chronology_new, by = 'sample_id')  %>% 
        mutate(lin_interp_age = if_else(sample_id %in% lI_chrono$sample_id, age, lin_interp_age),
               lin_interp_age_uncert_pos = if_else(sample_id %in% lI_chrono$sample_id, uncert_pos, lin_interp_age_uncert_pos),
               lin_interp_age_uncert_neg = if_else(sample_id %in% lI_chrono$sample_id, uncert_neg, lin_interp_age_uncert_neg)) %>% 
        select(-age, -uncert_pos, -uncert_neg)
      #print('lininterp')
    }
    
    if(y$copRa){
      setwd(file.path(y$working_directory, file_name, '/copRa'))
      copRa_chrono <- read.csv('copRa_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) 
      names(copRa_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
      SISAL_chronology_new <- full_join(copRa_chrono, SISAL_chronology_new, by = 'sample_id')  %>% 
        mutate(copRa_age = if_else(sample_id %in% copRa_chrono$sample_id, age, copRa_age),
               copRa_age_uncert_pos = if_else(sample_id %in% copRa_chrono$sample_id, uncert_pos, copRa_age_uncert_pos),
               copRa_age_uncert_neg = if_else(sample_id %in% copRa_chrono$sample_id, uncert_neg, copRa_age_uncert_neg)) %>% 
        select(-age, -uncert_pos, -uncert_neg)
      #print('lininterp')
    }
    
    if(y$Bchron){
      setwd(file.path(y$working_directory, file_name, '/Bchron'))
      bchron_chrono <- read.csv('bchron_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric'))
      names(bchron_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
      SISAL_chronology_new <- full_join(bchron_chrono, SISAL_chronology_new, by = 'sample_id') %>% 
        mutate(bchron_age = if_else(sample_id %in% bchron_chrono$sample_id, age, bchron_age),
               bchron_age_uncert_pos = if_else(sample_id %in% bchron_chrono$sample_id, uncert_pos, bchron_age_uncert_pos),
               bchron_age_uncert_neg = if_else(sample_id %in% bchron_chrono$sample_id, uncert_neg, bchron_age_uncert_neg)) %>% 
        select(-age, -uncert_pos, -uncert_neg)
    }
    
    if(y$StalAge){
      setwd(file.path(y$working_directory, file_name, '/StalAge'))
      stalage_chrono <- read.csv('StalAge_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) 
      names(stalage_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
      SISAL_chronology_new <- full_join(stalage_chrono, SISAL_chronology_new, by = 'sample_id')  %>% 
        mutate(StalAge_age = if_else(sample_id %in% stalage_chrono$sample_id, age, StalAge_age),
               StalAge_age_uncert_pos = if_else(sample_id %in% stalage_chrono$sample_id, uncert_pos, StalAge_age_uncert_pos),
               StalAge_age_uncert_neg = if_else(sample_id %in% stalage_chrono$sample_id, uncert_neg, StalAge_age_uncert_neg)) %>% 
        select(-age, -uncert_pos, -uncert_neg)
      #print('bacon')
    }
    
    if(y$Bacon){
      setwd(file.path(y$working_directory, file_name, '/Bacon_runs'))
      bacon_chrono <- read.csv('bacon_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) #
      names(bacon_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
      SISAL_chronology_new <- full_join(bacon_chrono, SISAL_chronology_new, by = 'sample_id') %>% 
        mutate(bacon_age = if_else(sample_id %in% bacon_chrono$sample_id, age, bacon_age),
               bacon_age_uncert_pos = if_else(sample_id %in% bacon_chrono$sample_id, uncert_pos, bacon_age_uncert_pos),
               bacon_age_uncert_neg = if_else(sample_id %in% bacon_chrono$sample_id, uncert_neg, bacon_age_uncert_neg)) %>% 
        select(-age, -uncert_pos, -uncert_neg)
      #print('bchron')
    }
    
  }
  return(SISAL_chronology_new)
}


merge_SISAL_chrono_2 <- function(runFile, SISAL_chronology_new, entity){
  
  for(i in runFile$entity_id) {
    print(i)
    y <- runFile %>% filter(entity_id == i) 
    entity_name <- (entity %>% filter(entity_id == i))$entity_name
    file_name <- paste(i,"-",entity_name, sep = '')
    
    if(y$linReg){
      setwd(file.path(y$working_directory, file_name, '/linReg'))
      if(file.exists('linReg_chronology.csv')){
        lR_chrono <- read.csv('linReg_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) 
        names(lR_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        SISAL_chronology_new <- full_join(lR_chrono, SISAL_chronology_new, by = 'sample_id') %>% 
          mutate(lin_reg_age = if_else(sample_id %in% lR_chrono$sample_id, age, lin_reg_age),
                 lin_reg_age_uncert_pos = if_else(sample_id %in% lR_chrono$sample_id, uncert_pos, lin_reg_age_uncert_pos),
                 lin_reg_age_uncert_neg = if_else(sample_id %in% lR_chrono$sample_id, uncert_neg, lin_reg_age_uncert_neg)) %>% 
          select(-age, -uncert_pos, -uncert_neg)
      }
    }
    
    if(y$linInterp){
      setwd(file.path(y$working_directory, file_name, '/linInterp'))
      if(file.exists('linInt_chronology.csv')){
        lI_chrono <- read.csv('linInt_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) 
        names(lI_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        SISAL_chronology_new <- full_join(lI_chrono, SISAL_chronology_new, by = 'sample_id')  %>% 
          mutate(lin_interp_age = if_else(sample_id %in% lI_chrono$sample_id, age, lin_interp_age),
                 lin_interp_age_uncert_pos = if_else(sample_id %in% lI_chrono$sample_id, uncert_pos, lin_interp_age_uncert_pos),
                 lin_interp_age_uncert_neg = if_else(sample_id %in% lI_chrono$sample_id, uncert_neg, lin_interp_age_uncert_neg)) %>% 
          select(-age, -uncert_pos, -uncert_neg)
      }
    }
    
    if(y$copRa){
      setwd(file.path(y$working_directory, file_name, '/copRa'))
      if(file.exists('copRa_chronology.csv')){
        copRa_chrono <- read.csv('copRa_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) 
        names(copRa_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        SISAL_chronology_new <- full_join(copRa_chrono, SISAL_chronology_new, by = 'sample_id')  %>% 
          mutate(copRa_age = if_else(sample_id %in% copRa_chrono$sample_id, age, copRa_age),
                 copRa_age_uncert_pos = if_else(sample_id %in% copRa_chrono$sample_id, uncert_pos, copRa_age_uncert_pos),
                 copRa_age_uncert_neg = if_else(sample_id %in% copRa_chrono$sample_id, uncert_neg, copRa_age_uncert_neg)) %>% 
          select(-age, -uncert_pos, -uncert_neg)
      }
    }
    
    if(y$Bchron){
      setwd(file.path(y$working_directory, file_name, '/Bchron'))
      if(file.exists('bchron_chronology.csv')){
        bchron_chrono <- read.csv('bchron_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric'))
        names(bchron_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        SISAL_chronology_new <- full_join(bchron_chrono, SISAL_chronology_new, by = 'sample_id') %>% 
          mutate(bchron_age = if_else(sample_id %in% bchron_chrono$sample_id, age, bchron_age),
                 bchron_age_uncert_pos = if_else(sample_id %in% bchron_chrono$sample_id, uncert_pos, bchron_age_uncert_pos),
                 bchron_age_uncert_neg = if_else(sample_id %in% bchron_chrono$sample_id, uncert_neg, bchron_age_uncert_neg)) %>% 
          select(-age, -uncert_pos, -uncert_neg)
      }
    }
    
    if(y$StalAge){
      setwd(file.path(y$working_directory, file_name, '/StalAge'))
      if(file.exists('StalAge_chronology.csv')){
        stalage_chrono <- read.csv('StalAge_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) 
        names(stalage_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        SISAL_chronology_new <- full_join(stalage_chrono, SISAL_chronology_new, by = 'sample_id')  %>% 
          mutate(StalAge_age = if_else(sample_id %in% stalage_chrono$sample_id, age, StalAge_age),
                 StalAge_age_uncert_pos = if_else(sample_id %in% stalage_chrono$sample_id, uncert_pos, StalAge_age_uncert_pos),
                 StalAge_age_uncert_neg = if_else(sample_id %in% stalage_chrono$sample_id, uncert_neg, StalAge_age_uncert_neg)) %>% 
          select(-age, -uncert_pos, -uncert_neg)
      }
    }
    
    if(y$Bacon){
      setwd(file.path(y$working_directory, file_name, '/Bacon_runs'))
      if(file.exists('bacon_chronology.csv')){
        bacon_chrono <- read.csv('bacon_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) #
        names(bacon_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        SISAL_chronology_new <- full_join(bacon_chrono, SISAL_chronology_new, by = 'sample_id') %>% 
          mutate(bacon_age = if_else(sample_id %in% bacon_chrono$sample_id, age, bacon_age),
                 bacon_age_uncert_pos = if_else(sample_id %in% bacon_chrono$sample_id, uncert_pos, bacon_age_uncert_pos),
                 bacon_age_uncert_neg = if_else(sample_id %in% bacon_chrono$sample_id, uncert_neg, bacon_age_uncert_neg)) %>% 
          select(-age, -uncert_pos, -uncert_neg)
      }
    }
    
  }
  return(SISAL_chronology_new)
}

get_growth_rates <- function(SISAL_chronology){
  SISAL_chronology_GR <- SISAL_chronology %>% mutate_at(vars(depth_sample, bacon_age, bacon_age_uncert_pos, bacon_age_uncert_neg, 
                                                             bchron_age, bchron_age_uncert_pos, bchron_age_uncert_neg, 
                                                             StalAge_age, StalAge_age_uncert_pos, StalAge_age_uncert_neg, 
                                                             lin_interp_age, lin_interp_age_uncert_pos, lin_interp_age_uncert_pos, 
                                                             lin_reg_age, lin_reg_age_uncert_pos, lin_reg_age_uncert_neg,
                                                             copRa_age, copRa_age_uncert_pos, copRa_age_uncert_neg), as.numeric) %>% 
    mutate(bacon_gr = (lead(depth_sample) - depth_sample)/(lead(bacon_age) -bacon_age),
           bchron_gr = (lead(depth_sample)-depth_sample)/(lead(bchron_age)-bchron_age),
           StalAge_gr = (lead(depth_sample)-depth_sample)/(lead(StalAge_age)-StalAge_age),
           copRa_gr = (lead(depth_sample)-depth_sample)/(lead(copRa_age)-copRa_age),
           lin_interp_gr = (lead(depth_sample)-depth_sample)/(lead(lin_interp_age)-lin_interp_age),
           lin_reg_gr = (lead(depth_sample)-depth_sample)/(lead(lin_reg_age)-lin_reg_age)) #%>% 
  #select(entity_id, sample_id, lin_reg_gr, lin_interp_gr, copRa_gr, bacon_gr, bchron_gr, StalAge_gr)
  
  return(SISAL_chronology_GR)
} 

s <- merge_SISAL_chrono_2(co1, SISAL_chronology_new, entity) 
#s_rr <-tryCatch({merge_SISAL_chrono(rr, SISAL_chronology_new, entity)})
s_dup <- s %>% group_by(entity_id) %>% get_dupes(sample_id) %>% distinct(entity_id) # duplicated sample_depths 
s_new <- s %>% group_by(entity_id) %>% distinct()
#s_new$depth_sample <- as.numeric(levels(s_new$depth_sample))[s_new$depth_sample]
s_new_fil <- s_new %>% filter(entity_id %in% r$entity_id) #%>% mutate_at(vars(depth_sample), as.numeric)
gr_v2 <- get_growth_rates(s_new_fil)
save(gr_v2, file = 'gr_v2.RData')
write.csv(r, 'runFile_corrected.csv', row.names = F)

SISAL_eval_rev <- s_new_fil %>% ungroup() %>% group_by(entity_id) %>% 
  mutate(., diff_lr = lead(lin_reg_age) - lin_reg_age,
         diff_li = lead(lin_interp_age) - lin_interp_age,
         diff_copRa = lead(copRa_age) - copRa_age,
         diff_bacon = lead(bacon_age) - bacon_age,
         diff_bchron = lead(bchron_age) - bchron_age,
         diff_stalage = lead(StalAge_age) - StalAge_age) %>%
  mutate(rev_lr = if_else( diff_lr< 0, FALSE, TRUE),
         rev_li = if_else( diff_li< 0, FALSE, TRUE),
         rev_copRa = if_else( diff_copRa< 0, FALSE, TRUE),
         rev_bacon = if_else( diff_bacon< 0, FALSE, TRUE),
         rev_bchron = if_else( diff_bchron< 0, FALSE, TRUE),
         rev_StalAge = if_else( diff_stalage< 0, FALSE, TRUE)) %>%
  group_by(entity_id) %>% 
  summarise(lR_reversal = if_else(all(is.na(rev_lr)),-1, if_else(any(!(rev_lr), na.rm = T), 0, 1)),
            lI_reversal = if_else(all(is.na(rev_li)),-1, if_else(any(!(rev_li), na.rm = T), 0, 1)),
            copRa_reversal = if_else(all(is.na(rev_copRa)),-1, if_else(any(!(rev_copRa), na.rm = T), 0, 1)),
            bacon_reversal = if_else(all(is.na(rev_bacon)),-1, if_else(any(!(rev_bacon), na.rm = T), 0, 1)),
            bchron_reversal = if_else(all(is.na(rev_bchron)),-1, if_else(any(!(rev_bchron), na.rm = T), 0, 1)),
            stalage_reversal = if_else(all(is.na(rev_StalAge)),-1, if_else(any(!(rev_StalAge), na.rm = T), 0, 1))) 

names(SISAL_eval_rev) <- c('entity_id', 'lin. Reg.', 'lin. Interp.', 'copRa', 'Bacon', 'Bchron', 'StalAge')
mine.data2 <- as.data.frame(gather(data = SISAL_eval_rev, key = Class, value = Fit, c(-1)))
mine.heatmap2 <- ggplot(data = mine.data2, mapping = aes(x = Class, y = factor(entity_id), fill = as.character(Fit))) + 
  geom_tile()  + 
  facet_grid( ~ Class,switch = 'x',scales = 'free_x', space = 'free_x') + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        strip.placement = 'outside', 
        strip.background = element_rect(fill = 'white',color = 'black')) + 
  xlab(label = 'Age Model') + 
  ylab(label = 'Entity ID') + 
  ggtitle(label = 'Filter for Reversals') +
  scale_fill_manual(guide = guide_legend(title = 'Fit yes-no'),values = c('0'='indianred', '1'='skyblue', '-1' = 'black'))
mine.heatmap2

#rev <- SISAL_chronology_eval_reversal %>% group_by(entity_id) %>% filter(any(`lin. Reg.` == 0 | `lin. Interp.` == 0 | copRa == 0 | Bacon == 0 | Bchron == 0 | StalAge == 0))

dt <- dating %>% filter(entity_id %in% s_new_fil$entity_id & date_type != 'Event; hiatus') %>% mutate_at(vars(depth_dating), as.numeric)

d <- data.frame(entity_id = dt$entity_id, 
                sample_id = dt$dating_id, 
                depth_sample = dt$depth_dating,
                lin_reg_age = dt$corr_age,
                lin_reg_age_uncert_pos = dt$corr_age_uncert_pos,
                lin_reg_age_uncert_neg = dt$corr_age_uncert_neg,
                lin_interp_age = dt$corr_age,
                lin_interp_age_uncert_pos = dt$corr_age_uncert_pos,
                lin_interp_age_uncert_neg = dt$corr_age_uncert_neg,
                copRa_age = dt$corr_age,
                copRa_age_uncert_pos = dt$corr_age_uncert_pos,
                copRa_age_uncert_neg = dt$corr_age_uncert_neg,
                bacon_age = dt$corr_age,
                bacon_age_uncert_pos = dt$corr_age_uncert_pos,
                bacon_age_uncert_neg = dt$corr_age_uncert_neg,
                bchron_age = dt$corr_age,
                bchron_age_uncert_pos = dt$corr_age_uncert_pos,
                bchron_age_uncert_neg = dt$corr_age_uncert_neg,
                StalAge_age = dt$corr_age,
                StalAge_age_uncert_pos = dt$corr_age_uncert_pos,
                StalAge_age_uncert_neg = dt$corr_age_uncert_neg)

dt_count <- dt %>% group_by(entity_id) %>% count() 
d_fil <- d %>% filter(entity_id %in% s_new_fil$entity_id) %>% mutate_at(vars(depth_sample), as.numeric)

SISAL_eval_fit <- bind_rows(s_new_fil, d_fil) %>% ungroup() %>% group_by(entity_id) %>% dplyr::arrange(., depth_sample, .by_group =T) %>%
  mutate(rev_lr = if_else(abs(lin_reg_age-lead(lin_reg_age)) < (2*lin_reg_age_uncert_pos + 2*lead(lin_reg_age_uncert_neg)), TRUE, FALSE),
         rev_li = if_else(abs(lin_interp_age-lead(lin_interp_age)) < (2*lin_interp_age_uncert_pos + 2*lead(lin_interp_age_uncert_neg)), TRUE, FALSE),
         rev_copRa = if_else(abs(copRa_age-lead(copRa_age)) < (2*copRa_age_uncert_pos + 2*lead(copRa_age_uncert_neg)), TRUE, FALSE),
         rev_bacon = if_else(abs(bacon_age-lead(bacon_age)) < (2*bacon_age_uncert_pos + 2*lead(bacon_age_uncert_neg)), TRUE, FALSE),
         rev_bchron = if_else(abs(bchron_age-lead(bchron_age)) < (2*bchron_age_uncert_pos + 2*lead(bchron_age_uncert_neg)), TRUE, FALSE),
         rev_StalAge = if_else(abs(StalAge_age-lead(StalAge_age)) < (2*StalAge_age_uncert_pos + 2*lead(StalAge_age_uncert_neg)), TRUE, FALSE)) %>%
  summarise(lR_count = sum(!rev_lr, na.rm = T),
            lI_count = sum(!rev_li, na.rm = T),
            copRa_count = sum(!rev_copRa, na.rm = T),
            bacon_count = sum(!rev_bacon, na.rm = T),
            bchron_count = sum(!rev_bchron, na.rm = T),
            stalage_count = sum(!rev_StalAge, na.rm = T)) %>%
  left_join(., dt_count, by = "entity_id") %>% group_by(entity_id) %>% summarise(lR_fitness = if_else(lR_count <= 0.3*n, 1, 0),
                                                                                 lI_fitness = if_else(lI_count <= 0.3*n, 1, 0),
                                                                                 copRa_fitness = if_else(copRa_count <= 0.3*n, 1, 0),
                                                                                 bacon_fitness = if_else(bacon_count <= 0.3*n, 1, 0),
                                                                                 bchron_fitness = if_else(bchron_count <= 0.3*n, 1, 0),
                                                                                 stalage_fitness = if_else(stalage_count <= 0.3*n, 1, 0))


names(SISAL_eval_fit) <- c('entity_id', 'lin. Reg.', 'lin. Interp.', 'copRa', 'Bacon', 'Bchron', 'StalAge')
mine.data <- gather(data = SISAL_eval_fit, key = Class, value = Fit, c(-1))
mine.data <- as.data.frame(mine.data)

mine.heatmap1 <- ggplot(data = mine.data, mapping = aes(x = Class, y = factor(entity_id), fill = as.character(Fit))) + 
  geom_tile()  + 
  facet_grid( ~ Class,switch = 'x',scales = 'free_x', space = 'free_x') + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        strip.placement = 'outside', 
        strip.background = element_rect(fill = 'white',color = 'black')) + 
  xlab(label = 'Age Model') + 
  ylab(label = 'Entity ID') + 
  ggtitle(label = 'Filter for flexibility/fitness') +
  scale_fill_manual(guide = guide_legend(title = 'Fit yes-no'),values = c('0'='indianred', '1'='skyblue'))
mine.heatmap1

dt_new <- dt %>% as_tibble() %>% mutate(depth_new = (lead(depth_dating)-depth_dating)/2+depth_dating)
d_new <- data.frame(sample_id = dt_new$dating_id, depth_new = dt_new$depth_new)

SISAL_IQR_check <- bind_rows(s_new_fil, d_fil) %>% ungroup() %>% group_by(entity_id) %>% dplyr::arrange(., depth_sample, .by_group =T) %>% left_join(., d_new, by = "sample_id") %>% 
  mutate(depth_new_d = sapply(depth_new, function(x) depth_sample[order(abs(x - depth_sample))][1])) %>% 
  mutate(depth_new_d = if_else(!is.na(depth_new), depth_new_d,  NA_real_)) %>% nest() %>% 
  mutate(filter = purrr::map(.$data, function(x) filter(x, !is.na(depth_new_d))$depth_new_d)) %>% 
  mutate(proxies = purrr::map2(.$data, .$filter, function(x, y) filter(x, depth_sample %in% y & is.na(depth_new) & is.na(depth_new_d)) %>% mutate( 
    lR_iqr = (lin_reg_age_uncert_pos+lin_reg_age_uncert_neg)/2 , 
    lI_iqr =  (lin_interp_age_uncert_pos+lin_interp_age_uncert_neg)/2,
    copRa_iqr = (copRa_age_uncert_pos+copRa_age_uncert_neg)/2,
    bacon_iqr =  (bacon_age_uncert_pos+bacon_age_uncert_neg)/2, 
    bchron_iqr =  (bchron_age_uncert_pos+bchron_age_uncert_neg)/2, 
    StalAge_iqr =  (StalAge_age_uncert_pos+StalAge_age_uncert_neg)/2) %>% select(sample_id, depth_sample, lR_iqr, lI_iqr, copRa_iqr, bacon_iqr, bchron_iqr, StalAge_iqr))) %>% 
  mutate(filter2 = purrr::map(.$data, function(x) filter(x, !is.na(depth_new))$sample_id)) %>% 
  mutate(dating =  purrr::map2(.$data, .$filter2, function(x, y) filter(x, sample_id %in% y) %>% mutate( 
    date_iqr = (lin_reg_age_uncert_pos+lin_reg_age_uncert_neg)/2) %>% select(sample_id, depth_new_d, date_iqr))) %>%
  mutate(final = purrr::map2(.$proxies, .$dating, function(x,y) left_join(y %>% rename(depth_join = depth_new_d),x %>% rename(depth_join = depth_sample), by= 'depth_join'))) %>%
  unnest(final) %>% print() %>%
  mutate(lR = if_else((date_iqr + lead(date_iqr))/2<lR_iqr, 1,0),
         lI = if_else((date_iqr + lead(date_iqr))/2<lI_iqr , 1,0),
         copRa = if_else((date_iqr + lead(date_iqr))/2<copRa_iqr , 1,0),
         bacon = if_else((date_iqr + lead(date_iqr))/2< bacon_iqr , 1,0),
         bchron = if_else((date_iqr + lead(date_iqr))/2<bchron_iqr , 1,0),
         stalage = if_else((date_iqr + lead(date_iqr))/2<StalAge_iqr, 1,0)) %>% group_by(entity_id) %>% 
  summarise(general = n(),
            lR_iqr = if_else(sum(lR, na.rm = T)>= 0.6*general, 1,0),
            lI_iqr = if_else(sum(lI, na.rm = T)>= 0.6*general, 1, 0),
            copRa_iqr = if_else(sum(copRa, na.rm = T)>= 0.6*general , 1,0),
            bacon_iqr= if_else(sum(bacon, na.rm = T)>= 0.6*general,1,0),
            bchron_iqr = if_else(sum(bchron, na.rm = T)>= 0.6*general,1,0),
            StalAge_iqr = if_else(sum(stalage, na.rm = T)>= 0.6*general,1,0)) %>% 
  select(entity_id, lR_iqr, lI_iqr, copRa_iqr, bacon_iqr, bchron_iqr, StalAge_iqr)

names(SISAL_IQR_check) <- c('entity_id', 'lin. Reg.', 'lin. Interp.', 'copRa', 'Bacon', 'Bchron', 'StalAge')
mine.data3 <- gather(data = SISAL_IQR_check, key = Class, value = Fit, c(-1))
mine.data3 <- as.data.frame(mine.data3)

mine.heatmap3 <- ggplot(data = mine.data3, mapping = aes(x = Class, y = factor(entity_id), fill = as.character(Fit))) + 
  geom_tile()  + 
  facet_grid( ~ Class,switch = 'x',scales = 'free_x', space = 'free_x')Moment + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        strip.placement = 'outside', 
        strip.background = element_rect(fill = 'white',color = 'black')) + 
  xlab(label = 'Age Model') + 
  ylab(label = 'Entity ID') + 
  ggtitle(label = 'Filter for IQR (increase between dates)') +
  scale_fill_manual(guide = guide_legend(title = 'Fit yes-no'),values = c('0'='indianred', '1'='skyblue'))
mine.heatmap3


SISAL_IQR <- s_new_fil %>% ungroup() %>% group_by(entity_id) %>% arrange(., depth_sample, .by_group =T) %>%
  mutate(lR_iqr = (lin_reg_age_uncert_pos+lin_reg_age_uncert_neg)/2 , 
         lI_iqr =  (lin_interp_age_uncert_pos+lin_interp_age_uncert_neg)/2, 
         copRa_iqr =  (copRa_age_uncert_pos+copRa_age_uncert_neg)/2, 
         bacon_iqr =  (bacon_age_uncert_pos+bacon_age_uncert_neg)/2, 
         bchron_iqr =  (bchron_age_uncert_pos+bchron_age_uncert_neg)/2, 
         StalAge_iqr =  (StalAge_age_uncert_pos+StalAge_age_uncert_neg)/2)

SISAL_young <- SISAL_IQR %>% group_by(entity_id) %>% mutate(general = n()/3) %>%
  slice(.,1:general) %>% group_by(entity_id) %>% summarise(lR_young = mean(lR_iqr, na.rm = T),
                                                           lI_young = mean(lI_iqr, na.rm = T),
                                                           copRa_young = mean(copRa_iqr, na.rm = T),
                                                           bacon_young = mean(bacon_iqr, na.rm = T),
                                                           bchron_young = mean(bchron_iqr, na.rm = T),
                                                           stalage_young = mean(StalAge_iqr, na.rm = T))

SISAL_old <- SISAL_IQR %>% group_by(entity_id) %>% mutate(general = n()/3) %>%
  slice(.,general*2:n()) %>% group_by(entity_id) %>% summarise(lR_old = mean(lR_iqr, na.rm = T),
                                                               lI_old = mean(lI_iqr, na.rm = T),
                                                               copRa_old = mean(copRa_iqr, na.rm = T),
                                                               bacon_old = mean(bacon_iqr, na.rm = T),
                                                               bchron_old = mean(bchron_iqr, na.rm = T),
                                                               stalage_old = mean(StalAge_iqr, na.rm = T))

SISAL_IQR_check2 <- left_join(SISAL_young, SISAL_old, by = "entity_id") %>% mutate(lR_iqr = if_else(lR_young < lR_old , 1, 0),
                                                                                   lI_iqr = if_else(lI_young < lI_old, 1,0),
                                                                                   copRa_iqr = if_else(copRa_young < lI_old, 1,0),
                                                                                   bacon_iqr = if_else(bacon_young < bacon_old,1,0),
                                                                                   bchron_iqr = if_else(bchron_young < bchron_old,1,0),
                                                                                   StalAge_iqr = if_else(stalage_young < stalage_old, 1,0)) %>%
  mutate(lR_iqr_2 = if_else(is.na(lR_young) & is.na(lR_old), -1, lR_iqr),
         lI_iqr_2 = if_else(is.na(lI_young) & is.na(lI_old), -1, lI_iqr),
         copRa_iqr_2 = if_else(is.na(copRa_young) & is.na(copRa_old), -1, copRa_iqr),
         bacon_iqr_2 = if_else(is.na(bacon_young) & is.na(bacon_old), -1, bacon_iqr),
         bchron_iqr_2 = if_else(is.na(bchron_young) & is.na(bchron_old),-1,bchron_iqr),
         StalAge_iqr_2 = if_else(is.na(stalage_young) & is.na(stalage_old),-1,StalAge_iqr)) %>%
  select(entity_id, lR_iqr_2, lI_iqr_2, copRa_iqr_2, bacon_iqr_2, bchron_iqr_2, StalAge_iqr_2)

names(SISAL_IQR_check2) <- c('entity_id', 'lin. Reg.', 'lin. Interp.', 'copRa', 'Bacon', 'Bchron', 'StalAge')
data4 <- SISAL_IQR_check2
mine.data4 <- gather(data = data4, key = Class, value = Fit, c(-1))
mine.data4 <- as.data.frame(mine.data4)

mine.heatmap4 <- ggplot(data = mine.data4, mapping = aes(x = Class, y = factor(entity_id), fill = as.character(Fit))) + 
  geom_tile()  + 
  facet_grid( ~ Class,switch = 'x',scales = 'free_x', space = 'free_x') + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        strip.placement = 'outside', 
        strip.background = element_rect(fill = 'white',color = 'black')) + 
  xlab(label = 'Age Model') + 
  ylab(label = 'Entity ID') + 
  ggtitle(label = 'Filter for IQR (increase with depth)') +
  scale_fill_manual(guide = guide_legend(title = 'Fit yes-no'),values = c('0'='indianred', '1'='skyblue', '-1' = 'black'))
mine.heatmap4

SISAL_eval <- left_join(SISAL_eval_rev, SISAL_eval_fit, by = "entity_id") %>% left_join(., SISAL_IQR_check, by = "entity_id") %>%
  left_join(., SISAL_IQR_check2, by = "entity_id") %>%
  mutate(lR = if_else(lR_reversal == -1, -1, if_else(lR_reversal == 0, 0, 1 + lR_fitness + lR_iqr +lR_iqr_2)),
         lI = if_else(lI_reversal == -1,-1,if_else(lI_reversal == 0, 0, 1 + lI_fitness + lI_iqr + lI_iqr_2)),
         copRa = if_else(copRa_reversal == -1,-1,if_else(copRa_reversal == 0, 0, 1 + copRa_fitness + copRa_iqr + copRa_iqr_2)),
         bacon = if_else(bacon_reversal == -1, -1, if_else(bacon_reversal == 0, 0, 1 + bacon_fitness + bacon_iqr + bacon_iqr_2)),
         bchron = if_else(bchron_reversal == -1, -1, if_else(bchron_reversal == 0, 0, 1+ bchron_fitness + bchron_iqr + bchron_iqr_2)),
         stalage = if_else(stalage_reversal == -1,-1, if_else(stalage_reversal == 0, 0, 1 + stalage_fitness + StalAge_iqr + StalAge_iqr_2))) %>% 
  mutate(success = if_else(lR>0 | lI >0 | copRa >0 | bacon >0 | bchron>0 | stalage >0, 'success', '-1')) %>% 
  select(entity_id, lR, lI,copRa, bacon, bchron, stalage, success) 

data5 <- SISAL_eval %>% select(entity_id, lR, lI,copRa, bacon, bchron, stalage) %>% filter(entity_id %in% r2$entity_id)
#data5 <- SISAL_eval %>% select(entity_id, lR, lI,bacon, bchron, stalage) %>% filter(entity_id %in% r2$entity_id)
not <- data5 %>% filter(lR ==-1 & lI == -1 & bacon == -1 & bchron == -1 & stalage == -1)
#names(data5) <- c('entity_id', 'lin. Reg.', 'lin. Interp.','copRa', 'Bacon', 'Bchron', 'StalAge')
names(SISAL_eval) <- c('entity_id', 'lin. Reg.', 'lin. Interp.', 'copRa','Bacon', 'Bchron', 'StalAge', 'Success')
mine.data5 <- gather(data = SISAL_eval, key = Class, value = Fit, c(-1))
mine.data5 <- as.data.frame(mine.data5)
mine.data5$Fit <- factor(mine.data5$Fit, levels = c('-1', '0', '1', '2', '3', '4', 'NA', 'success'))

mine.heatmap5 <- ggplot(data = mine.data5, mapping = aes(x = Class, y = factor(entity_id), fill = Fit, color = Fit)) + 
  geom_tile()  + 
  facet_grid( ~ Class,switch = 'x',scales = 'free_x', space = 'free_x') + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        strip.placement = 'outside', 
        strip.background = element_rect(fill = 'white',color = 'black')) + 
  xlab(label = 'Age Model') + 
  ylab(label = 'Entity ID') + 
  scale_fill_manual(guide = guide_legend(title = 'Score'),values = c('-1' = 'black', '0'='indianred', '1'='dodgerblue4', 
                                                                     '2' = 'deepskyblue3', '3' = 'deepskyblue1', '4' = 'darkslategray1', 'success' =  "#43BF71FF")) +
  scale_color_manual(guide = guide_legend(title = 'Score'),values = c('-1' = 'black', '0'='indianred', '1'='dodgerblue4', 
                                                                     '2' = 'deepskyblue3', '3' = 'deepskyblue1', '4' = 'darkslategray1', 'success' =  "#43BF71FF"))
#ggtitle(label = 'Age Model Check') +    
mine.heatmap5

pdf('~/Documents/plots Msc Midterm/age-model-check.pdf', 6,4)
mine.heatmap5
dev.off()

system("cdo copy -timmean /stacydata/data/P2Fvar/curated_ESGF_replica/MPI-ESM-P/lgm/pr_Amon_MPI-ESM-P_lgm_r1i1p1_185001-194912.nc /home/ariana/Documents/pr_mean_bg/pr_Amon_MPI-ESM-P_lgm_r1i1p1_timmean.nc")
system("cdo copy -timmean /stacydata/data/P2Fvar/curated_ESGF_replica/MPI-ESM-P/piControl/pr_Amon_MPI-ESM-P_piControl_r1i1p1_185001-199912.nc /home/ariana/Documents/pr_mean_bg/pr_Amon_MPI-ESM-P_piControl_r1i1p1_timmean.nc")
#system("cdo copy -timmean /stacydata/data/P2Fvar/curated_ESGF_replica/MPI-ESM-P/piControl/preprocessed/pr_Amon_MPI-ESM-P_piControl_preprocessed_final50yrs.nc /home/ariana/Documents/pr_mean_bg/pr_Amon_MPI-ESM-P_piControl_r1i1p1_timmean2.nc")
system("cdo sub /home/ariana/Documents/pr_mean_bg/pr_Amon_MPI-ESM-P_piControl_r1i1p1_timmean.nc /home/ariana/Documents/pr_mean_bg/pr_Amon_MPI-ESM-P_lgm_r1i1p1_timmean.nc /home/ariana/Documents/pr_mean_bg/pr_Amon_MPI-ESM-P_r1i1p1_anomaly.nc")
system("cdo div /home/ariana/Documents/pr_mean_bg/pr_Amon_MPI-ESM-P_r1i1p1_anomaly.nc /home/ariana/Documents/pr_mean_bg/pr_Amon_MPI-ESM-P_lgm_r1i1p1_timmean.nc /home/ariana/Documents/pr_mean_bg/pr_Amon_MPI-ESM-P_r1i1p1_anomaly_relative.nc")

gr_v2 <- get_growth_rates(s_new_fil)
gr_copRa <- gr_v2 %>% select(entity_id, copRa_age, copRa_gr) %>% group_by(entity_id) %>% mutate(copRa_gr = abs(copRa_gr)) %>% rename(gr = copRa_gr, age = copRa_age)
int <- gr_copRa %>% group_by(entity_id) %>%filter(median(age/1000, na.rm = T) < 30 | entity_id %in% c(367,115,53,107,416,54,312,124,1)) %>% distinct(entity_id)
gr_v2_lgm_holo <- gr_v2 %>% filter(entity_id %in% int$entity_id)
gr_copRa_lgm_holo <- gr_copRa %>% filter(entity_id %in% int$entity_id)
#int <- gr_copRa %>% filter(entity_id %in% c(367,115,53,107,416,54,312,124,1)) %>% distinct(entity_id)
  

get_mean <- function(gr, p1, p2, site_tb){
  gr_holo_mean <- gr %>% group_by(entity_id) %>% filter(age <= p1[2] & age >= p1[1]) %>% 
    summarize(holo_mean = mean(gr, na.rm = T))
  gr_lgm_mean <- gr %>% group_by(entity_id) %>% filter(age <= p2[2] & age >= p2[1]) %>% 
    summarize(lgm_mean = mean(gr, na.rm = T))
  gr_holo_lgm <- gr_holo_mean %>% filter(entity_id %in% gr_lgm_mean$entity_id) %>% 
    left_join(., gr_lgm_mean %>% filter(entity_id %in% gr_holo_mean$entity_id), by = 'entity_id') %>%
    mutate(anomaly = holo_mean - lgm_mean) 
  long_lat_lgm_holo <- site_tb %>% filter(entity_id %in% gr_holo_lgm$entity_id) %>% select(entity_id, latitude, longitude)
  gr_full<- left_join(gr_holo_lgm, long_lat_lgm_holo, by = 'entity_id') %>% mutate(long = if_else(longitude <0, longitude +360, longitude))
  
  return(gr_full)
}

plotBasicMap <- function(gr, nc){
  library(ncdf4)
  library(fields)
  lon <- ncvar_get(nc, 'lon')
  lat <- ncvar_get(nc, 'lat')
  pr <- ncvar_get(nc, 'pr')
  
  sisal_range <- c(-1,1)*max(abs(gr$anomaly))
  assign("col_anomaly",colorRampPalette(c("blue3","white","red3"),space="rgb")(33),envir=globalenv()) #Anomalies from reference
  
  image.plot(lon,lat,pr,col=col_anomaly,zlim=c(-1,1)*max(abs(pr)));map(add=T,wrap=c(0,360))
  points(gr$long,gr$latitude,col="black",
         bg=col_anomaly[round((gr$anomaly-sisal_range[1])/(sisal_range[2]-sisal_range[1])*(length(col_anomaly)-1)+1,0)],pch=21,cex=1.25)
  
}


setwd("~/Documents/growth_rates_plots")
plotfct(gr_copRa_lgm_holo, 5, F) %>% 
  ggsave(filename = 'gr_SISAL_v2_copRa.pdf', width = 30, height = 20, units = 'cm')

copRa_eID <- gr_copRa_lgm_holo %>% distinct(entity_id)
long_lat_lgm_holo <- site_tb %>% filter(entity_id %in% copRa_eID$entity_id) %>% select(entity_id, latitude, longitude)
gr_copRa_holo_mean <- gr_copRa_lgm_holo %>% group_by(entity_id) %>% filter(age <= 10000) %>% 
  summarize(holo_mean = mean(gr, na.rm = T))
gr_copRa_lgm_mean <- gr_copRa_lgm_holo %>% group_by(entity_id) %>% filter(age <= 29000 & age >= 15000) %>% 
  summarize(lgm_mean = mean(gr, na.rm = T))


gr_copRa_holo_lgm <- gr_copRa_holo_mean %>% filter(entity_id %in% gr_copRa_lgm_mean$entity_id) %>% 
  left_join(., gr_copRa_lgm_mean %>% filter(entity_id %in% gr_copRa_holo_mean$entity_id), by = 'entity_id') %>%
  mutate(anomaly = holo_mean - lgm_mean) %>% filter(!(entity_id %in% c(16,563))) 
gr_copRa <- gr_v2 %>% select(entity_id, copRa_age, copRa_gr) %>% group_by(entity_id) %>% mutate(copRa_gr = abs(copRa_gr)) %>% rename(gr = copRa_gr, age = copRa_age)
gr_copRa_full <- get_mean(gr_copRa, c(-50,10000),c(19000,29000), site_tb) %>% filter(entity_id != 115)

gr_copRa_full<- left_join(gr_copRa_holo_lgm, long_lat_lgm_holo %>% filter(entity_id %in% gr_copRa_holo_lgm$entity_id), by = 'entity_id') %>% 
  mutate(long = if_else(longitude <0, longitude +360, longitude))


gr_bacon <- gr_v2 %>% select(entity_id, bacon_age, bacon_gr) %>% group_by(entity_id) %>%
  mutate(bacon_gr = abs(bacon_gr)) %>% rename(gr = bacon_gr, age = bacon_age)
gr_bacon_full <- get_mean(gr_bacon, c(-50,10000),c(19000,29000), site_tb) %>% filter(entity_id != 115)

gr_lR <- gr_v2 %>% select(entity_id, lin_reg_age, lin_reg_gr) %>% group_by(entity_id) %>%
  mutate(lin_reg_gr = abs(lin_reg_gr)) %>% rename(gr = lin_reg_gr, age = lin_reg_age)
gr_lR_full <- get_mean(gr_lR, c(-50,10000),c(19000,29000), site_tb) %>% filter(entity_id != 115)

gr_bchron <- gr_v2 %>% select(entity_id, bchron_age, bchron_gr) %>% group_by(entity_id) %>%
  mutate(bchron_gr = abs(bchron_gr)) %>% rename(gr = bchron_gr, age = bchron_age)
gr_bchron_full <- get_mean(gr_bchron, c(-50,10000),c(19000,29000), site_tb) %>% filter(entity_id != 115)

gr_StalAge <- gr_v2 %>% select(entity_id, StalAge_age, StalAge_gr) %>% group_by(entity_id) %>%
  mutate(StalAge_gr = abs(StalAge_gr)) %>% rename(gr = StalAge_gr, age = StalAge_age)
gr_StalAge_full <- get_mean(gr_StalAge, c(-50,10000),c(19000,29000), site_tb) %>% filter(entity_id != 115)

gr_hl <- gr_v2 %>% filter(entity_id %in% gr_bacon_full$entity_id | entity_id %in% gr_bchron_full$entity_id | entity_id %in% gr_copRa_full$entity_id |
                            entity_id %in% gr_StalAge_full$entity_id | entity_id %in% gr_lR_full$entity_id) %>% distinct(entity_id)

gr_am_fit <- SISAL_eval_fit %>% filter(entity_id %in% gr_hl$entity_id)

gr_am <- SISAL_eval %>% filter(entity_id %in% gr_hl$entity_id) %>% 
  mutate(bacon_new =  if_else(gr_am_fit$Bacon == 0, 0,bacon),
         copRa_new =  if_else(gr_am_fit$copRa == 0, 0,copRa),
         stalage_new =  if_else(gr_am_fit$StalAge == 0, 0,stalage),
         lR_new =  if_else(gr_am_fit$`lin. Reg.` == 0, 0,lR),
         bchron_new =  if_else(gr_am_fit$Bchron == 0, 0,bchron)) %>% 
  filter(entity_id %in% c(587,523,512,500,353,320,319,312,305,26,219,175,124,1,54,62)) 

gr_full <- gr_v2 %>% filter(entity_id %in% c(587,523,512,500,353,320,319,312,305,26,219,175,124,1,54,62,510,237,593)) %>%
  mutate(age = if_else(entity_id %in% c(1,26,54,319,320,500,353,593,510,237), copRa_age, NA_real_),
         gr = if_else(entity_id %in% c(1,26,54,319,320,500,353,593,510,237), copRa_gr, NA_real_),
         age_model = if_else(entity_id %in% c(1,26,54,319,320,500), 'copRa', 'NA')) %>%
  mutate(age = if_else(entity_id %in% c(62,124,175,305,312,523,587), bchron_age, age),
         gr = if_else(entity_id %in% c(62,124,175,305,312,523,587), bchron_gr, gr),
         age_model = if_else(entity_id %in% c(62,124,175,305,312,353,523,587), 'Bchron', age_model)) %>%
  mutate(age = if_else(entity_id %in% c(219,512), bacon_age, age),
         gr = if_else(entity_id %in% c(219,512), bacon_gr, gr),
         age_model = if_else(entity_id %in% c(219,512), 'Bacon', age_model)) %>% select(sample_id, entity_id, age, gr, age_model)

gr_full_hl <- get_mean(gr_full, c(-50,8000),c(19000,27000), site_tb) %>% mutate(age_model = if_else(entity_id %in% c(1,26,54,319,320,500,353,593,510,237), 'copRa', 
                                                                                                     if_else(entity_id %in% c(62,124,175,305,312,523,587), 'Bchron',
                                                                                                             if_else(entity_id %in% c(219,512), 'Bacon','NA')))) %>% 
  mutate(age_model_new = if_else(age_model == 'copRa', 21, if_else(age_model == 'Bchron',22, if_else(age_model == 'Bacon',23,0))))# %>%
  mutate(ano_rel = anomaly/lgm_mean) %>% mutate(ano_norm = anomaly/max(abs(anomaly)))


gr_full_hl_new <- rbind(gr_full_hl, orig_gr_hl, gr_50)


pr_m_new <- pr_m/max(abs(pr_m))

get_mean <- function(gr, p1, p2, site_tb){
  gr_holo_mean <- gr %>% group_by(entity_id) %>% filter(age <= p1[2] & age >= p1[1]) %>% 
    summarize(holo_mean = mean(gr, na.rm = T))
  gr_lgm_mean <- gr %>% group_by(entity_id) %>% filter(age <= p2[2] & age >= p2[1]) %>% 
    summarize(lgm_mean = mean(gr, na.rm = T))
  gr_holo_lgm <- gr_holo_mean %>% filter(entity_id %in% gr_lgm_mean$entity_id) %>% 
    left_join(., gr_lgm_mean %>% filter(entity_id %in% gr_holo_mean$entity_id), by = 'entity_id') %>%
    mutate(anomaly = (holo_mean - lgm_mean)) 
  long_lat_lgm_holo <- site_tb %>% filter(entity_id %in% gr_holo_lgm$entity_id) %>% select(entity_id, latitude, longitude)
  gr_full<- left_join(gr_holo_lgm, long_lat_lgm_holo, by = 'entity_id') %>% mutate(long = if_else(longitude <0, longitude +360, longitude))
  
  return(gr_full)
}


library(ncdf4)
library(fields)
nc <- nc_open('/home/ariana/Documents/pr_mean_bg/pr_Amon_MPI-ESM-P_r1i1p1_anomaly.nc')
lon <- ncvar_get(nc, 'lon')
lat <- ncvar_get(nc, 'lat')
pr <- ncvar_get(nc, 'pr')

sisal_range <- c(-1,1)*max(abs(gr_full_hl$anomaly))
assign("col_anomaly",colorRampPalette(c("blue3","white","red3"),space="rgb")(33),envir=globalenv()) #Anomalies from reference

setwd("/home/ariana/Documents/growth_rates_plots")
pdf('/home/ariana/Documents/growth_rates_plots/gr.pdf', 6, 4)
image.plot(lon,lat,pr,col=col_anomaly,zlim=c(-1,1)*max(abs(pr)));map(add=T,wrap=c(0,360))
points(gr_full_hl$long,gr_full_hl$latitude,col="black",
       bg=col_anomaly[round((gr_full_hl$anomaly-sisal_range[1])/(sisal_range[2]-sisal_range[1])*(length(col_anomaly)-1)+1,0)],pch=gr_full_hl$age_model_new,cex=1.25)
legend("bottomleft",legend=c('copRa', 'Bacon','Bchron'), pch=c(21,23,22), cex = 0.7, bg = 'white')
dev.off()

nc_rel <- nc_open('/home/ariana/Documents/pr_mean_bg/pr_Amon_MPI-ESM-P_r1i1p1_anomaly_relative.nc')
lon_rel <- ncvar_get(nc_rel, 'lon')
lat_rel <- ncvar_get(nc_rel, 'lat')
pr_rel <- ncvar_get(nc_rel, 'pr')

sisal_range <- c(-1,1)*max(abs(gr_full_hl$ano_rel))

pdf('/home/ariana/Documents/growth_rates_plots/gr_rel.pdf', 6, 4)
image.plot(lon_rel,lat_rel,pr_rel,col=col_anomaly,zlim=c(-1,1)*max(abs(pr_rel)));map(add=T,wrap=c(0,360))
points(gr_full_hl$long,gr_full_hl$latitude,col="black",
       bg=col_anomaly[round((gr_full_hl$ano_rel-sisal_range[1])/(sisal_range[2]-sisal_range[1])*(length(col_anomaly)-1)+1,0)],pch=gr_full_hl$age_model_new,cex=1.25)
legend("bottomleft",legend=c('copRa', 'Bacon','Bchron'), pch=c(21,23,22), cex = 0.7, bg = 'white')
dev.off()


setwd("/home/ariana/Documents/growth_rates_plots")
pdf('bacon_gr.pdf', 6, 4)
plotBasicMap(gr_bacon_full, nc)
dev.off()

setwd("/home/ariana/Documents/growth_rates_plots")
pdf('lR_gr.pdf', 6, 4)
plotBasicMap(gr_lR_full, nc)
dev.off()

setwd("/home/ariana/Documents/growth_rates_plots")
pdf('bchron_gr.pdf', 6, 4)
plotBasicMap(gr_bchron_full, nc)
dev.off()


save(gr_copRa_full, file = 'gr_copRa_full.RData')


orig_gr <- sample_tb %>% filter(entity_id %in% c(163)) %>% mutate(gr = (lead(depth_sample)-depth_sample)/(lead(interp_age)-interp_age),
                                                                  age = interp_age,
                                                                  age_model = rep('original', length(gr))) %>%select(sample_id, entity_id, age, gr, age_model)
orig_gr_hl = get_mean(orig_gr, c(-50,8000),c(19000,27000), site_tb) %>% mutate(age_model = 'original',age_model_new = 24)



chrono510 <- read.csv("~/Documents/v2/510-Itzamna/copRa/copRa_chronology.csv", header = T, stringsAsFactors = F) %>% 
  left_join(., sample_tb %>% filter(entity_id == 510) %>% select(sample_id, depth_sample, entity_id)) %>% 
  mutate(gr = (lead(depth_sample)-depth_sample)/(lead(copRa_age)-copRa_age)) %>% rename(age = copRa_age)

gr_50 <- get_mean(chrono510, c(-50,8000),c(19000,27000), site_tb) %>% mutate(age_model = 'copRa', age_model_new=21)






gr_h_l <- full_join(gr_bacon_full %>% rename(bacon_holo = holo_mean, bacon_lgm = lgm_mean, bacon_ano = anomaly), 
                    gr_bchron_full %>% rename(bchron_holo = holo_mean, bchron_lgm = lgm_mean, bchron_ano = anomaly), by = c('longitude', 'latitude', 'long' ,'entity_id')) %>% 
  full_join(., gr_copRa_full %>% rename(copRa_holo = holo_mean, copRa_lgm = lgm_mean, copRa_ano = anomaly), by = c('longitude', 'latitude', 'long' ,'entity_id')) %>% 
  select(entity_id, long, longitude, latitude, everything())

pdf('comparison.pdf', 6,6)
layout(matrix(c(1,2,3), nrow=3, ncol=1, byrow = T))#, widths = c(rep(lcm(9),3)), heights = c(rep(lcm(6.5),3)))
par(mar = c(1,0,2,0), oma = c(2,1,2,1))
image.plot(lon,lat,pr,col=col_anomaly,zlim=c(-1,1)*max(abs(pr)));map(add=T,wrap=c(0,360))
points((gr_h_l%>% filter(!is.na(bacon_ano)))$long,(gr_h_l%>% filter(!is.na(bacon_ano)))$latitude,col="black",
       bg=col_anomaly[round(((gr_h_l%>% filter(!is.na(bacon_ano)))$bacon_ano-sisal_range[1])/(sisal_range[2]-sisal_range[1])*(length(col_anomaly)-1)+1,0)],pch=21,cex=1.25)
mtext('Bacon', side = 3)
image.plot(lon,lat,pr,col=col_anomaly,zlim=c(-1,1)*max(abs(pr)));map(add=T,wrap=c(0,360))
points((gr_h_l%>% filter(!is.na(bchron_ano)))$long,(gr_h_l%>% filter(!is.na(bchron_ano)))$latitude,col="black",
       bg=col_anomaly[round(((gr_h_l%>% filter(!is.na(bchron_ano)))$bchron_ano-sisal_range[1])/(sisal_range[2]-sisal_range[1])*(length(col_anomaly)-1)+1,0)],pch=21,cex=1.25)
mtext('Bchron', side = 3)
image.plot(lon,lat,pr,col=col_anomaly,zlim=c(-1,1)*max(abs(pr)));map(add=T,wrap=c(0,360))
points((gr_h_l %>% filter(!is.na(copRa_ano)))$long,(gr_h_l %>% filter(!is.na(copRa_ano)))$latitude,col="black",
       bg=col_anomaly[round(((gr_h_l %>% filter(!is.na(copRa_ano)))$copRa_ano-sisal_range[1])/(sisal_range[2]-sisal_range[1])*(length(col_anomaly)-1)+1,0)],pch=21,cex=1.25)
mtext('copRa', side = 3)
dev.off()  

names(SISAL_eval_rev) <- c('entity_id', 'lin. Reg.', 'lin. Interp.', 'copRa', 'Bacon', 'Bchron', 'StalAge')
mine.data <- as.data.frame(gather(data = SISAL_eval_rev, key = Class, value = Fit, c(-1)))
mine.heatmap <- ggplot(data = mine.data, mapping = aes(x = Class, y = factor(entity_id), fill = as.character(Fit))) + 
  geom_tile()  + 
  facet_grid( ~ Class,switch = 'x',scales = 'free_x', space = 'free_x') + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        strip.placement = 'outside', 
        strip.background = element_rect(fill = 'white',color = 'black')) + 
  xlab(label = 'Age Model') + 
  ylab(label = 'Entity ID') + 
  ggtitle(label = 'Filter for Reversals') +
  scale_fill_manual(guide = guide_legend(title = 'Fit yes-no'),values = c('0'='indianred', '1'='skyblue', '-1' = 'black'))
mine.heatmap

names(SISAL_eval_fit) <- c('entity_id', 'lin. Reg.', 'lin. Interp.', 'copRa', 'Bacon', 'Bchron', 'StalAge')
mine.data2 <- gather(data = SISAL_eval_fit, key = Class, value = Fit, c(-1))
mine.data2 <- as.data.frame(mine.data2)

mine.heatmap2 <- ggplot(data = mine.data2, mapping = aes(x = Class, y = factor(entity_id), fill = as.character(Fit))) + 
  geom_tile()  + 
  facet_grid( ~ Class,switch = 'x',scales = 'free_x', space = 'free_x') + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        strip.placement = 'outside', 
        strip.background = element_rect(fill = 'white',color = 'black')) + 
  xlab(label = 'Age Model') + 
  ylab(label = 'Entity ID') + 
  ggtitle(label = 'Filter for flexibility/fitness') +
  scale_fill_manual(guide = guide_legend(title = 'Fit yes-no'),values = c('0'='indianred', '1'='skyblue'))
mine.heatmap2

names(SISAL_IQR_check) <- c('entity_id', 'lin. Reg.', 'lin. Interp.', 'copRa', 'Bacon', 'Bchron', 'StalAge')
mine.data3 <- gather(data = SISAL_IQR_check, key = Class, value = Fit, c(-1))
mine.data3 <- as.data.frame(mine.data3)

mine.heatmap3 <- ggplot(data = mine.data3, mapping = aes(x = Class, y = factor(entity_id), fill = as.character(Fit))) + 
  geom_tile()  + 
  facet_grid( ~ Class,switch = 'x',scales = 'free_x', space = 'free_x') + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        strip.placement = 'outside', 
        strip.background = element_rect(fill = 'white',color = 'black')) + 
  xlab(label = 'Age Model') + 
  ylab(label = 'Entity ID') + 
  ggtitle(label = 'Filter for IQR (increase between dates)') +
  scale_fill_manual(guide = guide_legend(title = 'Fit yes-no'),values = c('0'='indianred', '1'='skyblue'))
mine.heatmap3

names(SISAL_IQR_check2) <- c('entity_id', 'lin. Reg.', 'lin. Interp.', 'copRa', 'Bacon', 'Bchron', 'StalAge')
mine.data4 <- gather(data = SISAL_IQR_check2, key = Class, value = Fit, c(-1))
mine.data4 <- as.data.frame(mine.data4)

mine.heatmap4 <- ggplot(data = mine.data4, mapping = aes(x = Class, y = factor(entity_id), fill = as.character(Fit))) + 
  geom_tile()  + 
  facet_grid( ~ Class,switch = 'x',scales = 'free_x', space = 'free_x') + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        strip.placement = 'outside', 
        strip.background = element_rect(fill = 'white',color = 'black')) + 
  xlab(label = 'Age Model') + 
  ylab(label = 'Entity ID') + 
  ggtitle(label = 'Filter for IQR (increase with depth)') +
  scale_fill_manual(guide = guide_legend(title = 'Fit yes-no'),values = c('0'='indianred', '1'='skyblue', '-1' = 'black'))
mine.heatmap4







