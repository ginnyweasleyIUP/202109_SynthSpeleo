merge_SISAL_chrono_final <- function(runFile, se, entity){
  
  SISAL_chronology_new <- data.frame(entity_id = se$entity_id, 
                                     sample_id = se$sample_id, 
                                     depth_sample = se$depth_sample,
                                     lin_reg_age = rep(NA_real_, length(se$sample_id)),
                                     lin_reg_age_uncert_pos = rep(NA_real_, length(se$sample_id)),
                                     lin_reg_age_uncert_neg = rep(NA_real_, length(se$sample_id)),
                                     lin_interp_age = rep(NA_real_, length(se$sample_id)),
                                     lin_interp_age_uncert_pos = rep(NA_real_, length(se$sample_id)),
                                     lin_interp_age_uncert_neg = rep(NA_real_, length(se$sample_id)),
                                     copRa_age = rep(NA_real_, length(se$sample_id)),
                                     copRa_age_uncert_pos = rep(NA_real_, length(se$sample_id)),
                                     copRa_age_uncert_neg = rep(NA_real_, length(se$sample_id)),
                                     StalAge_age = rep(NA_real_, length(se$sample_id)),
                                     StalAge_age_uncert_pos = rep(NA_real_, length(se$sample_id)),
                                     StalAge_age_uncert_neg = rep(NA_real_, length(se$sample_id)),
                                     bacon_age = rep(NA_real_, length(se$sample_id)),
                                     bacon_age_uncert_pos = rep(NA_real_, length(se$sample_id)),
                                     bacon_age_uncert_neg = rep(NA_real_, length(se$sample_id)),
                                     bchron_age = rep(NA_real_, length(se$sample_id)),
                                     bchron_age_uncert_pos = rep(NA_real_, length(se$sample_id)),
                                     bchron_age_uncert_neg = rep(NA_real_, length(se$sample_id))
  )
  
  for(i in runFile$entity_id) {
    print(i)
    y <- runFile %>% filter(entity_id == i) 
    entity_name <- (entity %>% filter(entity_id == i))$entity_name
    file_name <- paste(i,"-",entity_name, sep = '')
    
    if(y$linReg){
      setwd(file.path(y$working_directory, file_name, '/linReg'))
      if(file.exists('linReg_chronology.csv')){
        lR_chrono <- read.csv('linReg_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% distinct(sample_id, .keep_all = T)
        names(lR_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        SISAL_chronology_new <- full_join(lR_chrono, SISAL_chronology_new, by = 'sample_id') %>% 
          mutate(lin_reg_age = if_else(sample_id %in% lR_chrono$sample_id, age, lin_reg_age),
                 lin_reg_age_uncert_pos = if_else(sample_id %in% lR_chrono$sample_id, uncert_pos, lin_reg_age_uncert_pos),
                 lin_reg_age_uncert_neg = if_else(sample_id %in% lR_chrono$sample_id, uncert_neg, lin_reg_age_uncert_neg)) %>% 
          select(-age, -uncert_pos, -uncert_neg)
      } else {y$linReg <- F}
    }
    
    if(y$linInterp){
      setwd(file.path(y$working_directory, file_name, '/linInterp'))
      if(file.exists('linInt_chronology.csv')){
        lI_chrono <- read.csv('linInt_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% distinct(sample_id, .keep_all = T)
        names(lI_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        SISAL_chronology_new <- full_join(lI_chrono, SISAL_chronology_new, by = 'sample_id')  %>% 
          mutate(lin_interp_age = if_else(sample_id %in% lI_chrono$sample_id, age, lin_interp_age),
                 lin_interp_age_uncert_pos = if_else(sample_id %in% lI_chrono$sample_id, uncert_pos, lin_interp_age_uncert_pos),
                 lin_interp_age_uncert_neg = if_else(sample_id %in% lI_chrono$sample_id, uncert_neg, lin_interp_age_uncert_neg)) %>% 
          select(-age, -uncert_pos, -uncert_neg)
      } else { y$linInterp <- F}
    }
    
    if(y$copRa){
      setwd(file.path(y$working_directory, file_name, '/copRa'))
      if(file.exists('copRa_chronology.csv')){
        copRa_chrono <- read.csv('copRa_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% distinct(sample_id, .keep_all = T)
        names(copRa_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        SISAL_chronology_new <- full_join(copRa_chrono, SISAL_chronology_new, by = 'sample_id')  %>% 
          mutate(copRa_age = if_else(sample_id %in% copRa_chrono$sample_id, age, copRa_age),
                 copRa_age_uncert_pos = if_else(sample_id %in% copRa_chrono$sample_id, uncert_pos, copRa_age_uncert_pos),
                 copRa_age_uncert_neg = if_else(sample_id %in% copRa_chrono$sample_id, uncert_neg, copRa_age_uncert_neg)) %>% 
          select(-age, -uncert_pos, -uncert_neg)
      } else {y$copRa <- F}
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
      } else {y$Bchron <- F}
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
      } else {y$StalAge <- F}
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
      } else {y$Bacon <- F}
    }
    
  }
  return(list(runFile, SISAL_chronology_new))
}

merge_dating_chrono <- function(dating, runFile){
  
  dating_file <- dating %>% filter(entity_id %in% runFile$entity_id) %>%
    mutate(date_used_copRa = date_used,
           date_used_Bacon = date_used,
           date_used_Bchron = date_used,
           date_used_lin_reg = date_used,
           date_used_lin_interp = date_used,
           date_used_StalAge = date_used) %>% 
    select(entity_id, dating_id, date_used_lin_reg, date_used_lin_interp, date_used_copRa, date_used_StalAge, date_used_Bacon, date_used_Bchron)
  
}

setwd("~/SISAL Data/sisalv2/")
sample <- read.csv('sample.csv', header = T, stringsAsFactors = F) %>% mutate_at(vars(depth_sample), as.numeric)
entity <- read.csv('entity.csv', header = T, stringsAsFactors = F)
site <- read.csv('site.csv', header = T, stringsAsFactors = F)
dating <- read.csv('dating.csv', header = T, stringsAsFactors = F) 
dating <- dating %>% mutate_at(vars(dating_id, entity_id, depth_dating, dating_thickness, X14C_correction, corr_age, corr_age_uncert_pos, corr_age_uncert_neg), as.numeric)
orig2 <- read.csv('original_chronology.csv', header = T, stringsAsFactors = F)
site_tb <- left_join(site, entity, by = 'site_id')
sisal2 <- read.csv('sisal_chronology.csv', header = T, stringsAsFactors = F)
hiatus_tb <- read.csv("hiatus.csv", header = T, stringsAsFactors = F)

sample_tb <- sample %>% 
  mutate_at(vars(entity_id, sample_id, sample_thickness, depth_sample), as.numeric) 

entity_from_base <-  entity %>% filter(depth_ref == 'from base') %>% distinct(entity_id)
sample_from_base <- sample %>% filter(entity_id %in% entity_from_base$entity_id) %>% dplyr::select(entity_id,depth_sample) %>% group_by(entity_id) %>% summarise(max = max(depth_sample))
sampling_from_base <- full_join(sample, sample_from_base, by = 'entity_id') %>% 
  mutate_at(vars(max),as.numeric) %>%
  group_by(entity_id) %>% 
  mutate(depth_conv = if_else(entity_id %in% entity_from_base$entity_id, max-depth_sample, NA_real_)) %>% 
  mutate(depth_sample = if_else(!is.na(depth_conv), depth_conv, depth_sample)) %>%
  dplyr::select(-depth_conv) %>% arrange(., depth_sample, .by_group = T)
dating_from_base <- full_join(dating, sample_from_base, by = 'entity_id') %>% group_by(entity_id) %>% 
  mutate(depth_conv = if_else(entity_id %in% entity_from_base$entity_id, max-depth_dating, NA_real_)) %>% 
  mutate(depth_dating = if_else(!is.na(depth_conv), depth_conv, depth_dating)) %>%
  select(-depth_conv) %>% arrange(., depth_dating, .by_group = T) %>% ungroup()

se <- left_join(sampling_from_base, entity, by = 'entity_id')

run<- read.csv('/home/ariana/Documents/runv2.csv', header = T)
runFile <- tibble(entity_id = run$entity_id, Bacon = rep(T, length(run$entity_id)), Bchron = rep(T, length(run$entity_id)), 
                  StalAge = rep(T, length(run$entity_id)), linInterp = rep(T, length(run$entity_id)), copRa = rep(T, length(run$entity_id)),
                  linReg = rep(T, length(run$entity_id)), 
                  working_directory = rep("/stacywork/ariana/v2_90ci/", length(run$entity_id)), 
                  wd = rep('/home/ariana/SISAL Data/sisalv2', length(run$entity_id))) %>% 
  mutate(Bchron = if_else(entity_id %in% c(330,378,432,506,510), F, T)) %>% filter(!(entity_id %in% c(112,130,131))) %>% arrange(., entity_id)

#s <- merge_SISAL_chrono_final(runFile, se, entity)

#run_new <- as.data.frame(s[1])
#chrono <- as.data.frame(s[2])
#s_dup <- s %>% group_by(entity_id) %>% get_dupes(sample_id) %>% distinct(entity_id) # duplicated sample_depths 
#chrono<- chrono %>% group_by(entity_id) %>% distinct()
#chronology_table <- chrono %>% filter(!(sample_id %in% hiatus_tb$sample_id))

dating_new <- merge_dating_chrono(dating, run_new)

write.csv(chronology_table, "~/SISAL.chronology/isotopes_file.csv", row.names = F)
write.csv(dating_new, "~/SISAL.chronology/date_file.csv", row.names = F)
write.csv(run_new, "~/SISAL.chronology/runFile.csv", row.names = F)


runTest2 <- runFile %>% filter(entity_id > 329)

r <- merge_SISAL_ensemble_final(runFile, dating_from_base)
write.csv(r %>% select(-wd, -working_directory), "/home/ariana/SISAL.chronology/logFile.csv", row.names = F)

r_new <- merge_SISAL_ensemble_final(runFile, dating_from_base)
write.csv(r_new %>% select(-wd, -working_directory), "/home/ariana/SISAL/MC.ensemble.partially.current/logFile.csv", row.names = F)

merge_SISAL_ensemble_final <- function(runFile, dating){
  j <- 0
  entity <- read.csv('/home/ariana/SISAL Data/sisalv2/entity.csv', header = T, stringsAsFactors = F)
  
  for(i in runFile$entity_id) {
    print(i)
    j <- j+1
    y <- runFile %>% filter(entity_id == i) 
    entity_name <- (entity %>% filter(entity_id == i))$entity_name
    file_name <- paste(i,"-",entity_name, sep = '')
    file_name_wd <- file_name
    file_name <- list()
    
    dates <- dating %>% filter(entity_id == i) %>% select(dating_id, depth_dating, date_used, date_type, corr_age, corr_age_uncert_pos, corr_age_uncert_neg)
    
    setwd(file.path(y$working_directory, file_name_wd))
    hiatus_tb <- read.csv("hiatus.csv", header = T, colClasses = c('numeric', 'numeric'))
    info_tb <- read.csv("info.csv", header = T, stringsAsFactors = F)
    orig_tb <- read.csv('original_chronology.csv', header = T, stringsAsFactors = F)
    proxy_tb <- read.csv('proxy_data.csv', header = T, stringsAsFactors = F)
    
    file_name$hiatus <- hiatus_tb
    file_name$info <- info_tb
    file_name$origAM <- orig_tb
    file_name$proxy <- proxy_tb
    file_name$dating <- dates
    
    setwd(file.path(y$working_directory, file_name_wd, '/linReg'))
    if(file.exists('mc_linReg_ensemble.txt')){
      AM <- read.csv('linReg_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) 
      ens <- read.table('mc_linReg_ensemble.txt', header = F, colClasses = 'numeric', row.names = NULL)
      colnames(ens) <- NULL
      file_name$linReg$AM <- AM
      file_name$linReg$ens <- ens
    } else {
      file_name$linReg$AM <- list()
      file_name$linReg$ens <- list()
      runFile[j,7] <- F
    }
    
    setwd(file.path(y$working_directory, file_name_wd, '/linInterp'))
    if(file.exists('mc_linInt_ensemble.txt')){
      AM <- read.csv('linInt_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) 
      ens <- read.table('mc_linInt_ensemble.txt', header = F, colClasses = 'numeric', row.names = NULL)
      colnames(ens) <- NULL
      file_name$linInterp$AM <- AM
      file_name$linInterp$ens <- ens
    } else {
      file_name$linInterp$AM <- list()
      file_name$linInterp$ens <- list()
      runFile[j,5] <- F
    }
    
    setwd(file.path(y$working_directory, file_name_wd, '/copRa'))
    if(file.exists('mc_copRa_ensemble.txt')){
      AM <- read.csv('copRa_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) 
      ens <- read.table('mc_copRa_ensemble.txt', header = F, colClasses = 'numeric', row.names = NULL)
      colnames(ens) <- NULL
      file_name$copRa$AM <- AM
      file_name$copRa$ens <- ens
    } else {
      file_name$copRa$AM <- list()
      file_name$copRa$ens <- list()
      runFile[j,6] <- F
    }
    
    setwd(file.path(y$working_directory, file_name_wd, '/StalAge'))
    if(file.exists('StalAge_chronology.csv')){
      AM <- read.csv('StalAge_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) 
      file_name$StalAge$AM <- AM
      file_name$StalAge$ens <- list()
    } else {
      file_name$StalAge$AM <- list()
      file_name$StalAge$ens <- list()
      runFile[j, 4] <- F
    }
    
    if(y$Bchron){
      setwd(file.path(y$working_directory, file_name_wd, '/Bchron'))
      if(file.exists('bchron_ensemble.txt')){
        AM <- read.csv('bchron_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric'))
        ens <- t(read.table('bchron_ensemble.txt', header = F, colClasses = 'numeric', row.names = NULL))
        colnames(ens) <- NULL
        file_name$Bchron$AM <- AM
        file_name$Bchron$ens <- as.data.frame(ens)
      } else {
        file_name$Bchron$AM <- list()
        file_name$Bchron$AM <- list()
        runFile[j,3] <- F
      }
    } else {
      file_name$Bchron$AM <- list()
      file_name$Bchron$AM <- list()
      runFile[j,3] <- F
    }
    
    setwd(file.path(y$working_directory, file_name_wd, '/Bacon_runs'))
    if(file.exists('mc_bacon_ensemble.txt')){
      AM <- read.csv('bacon_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) 
      ens <- read.table('mc_bacon_ensemble.txt', header = F, colClasses = 'numeric', row.names = NULL)
      colnames(ens) <- NULL
      file_name$Bacon$AM <- AM
      file_name$Bacon$ens <- ens
    } else {
      file_name$Bacon$AM <- list()
      file_name$Bacon$ens <- list()
      runFile[j,2] <- F
    }
    
    save(file_name, file = file.path("/home/ariana/SISAL/MC.ensemble.partially.current/",paste(file_name_wd,".RData", sep = '')))
  }
  return(runFile)
}
