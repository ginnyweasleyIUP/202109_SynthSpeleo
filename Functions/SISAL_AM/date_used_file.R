existing <- left_join(sisalv2interim$sisal_chronology, sisalv2interim$sample, by = 'sample_id') %>% distinct(entity_id)
#
original <- sisalv2interim$dating %>% select(entity_id, dating_id, starts_with('date_used')) %>% left_join(., date_file_final, by = c('entity_id', 'dating_id')) %>%
  mutate(date_used_copRa = if_else(entity_id %in% existing$entity_id & date_used_copRa == 'NULL', date_used_COPRA, date_used_copRa),
         date_used_lin_interp = if_else(entity_id %in% existing$entity_id & date_used_lin_interp == 'NULL', date_used_linear, date_used_lin_interp),
         date_used_lin_reg = if_else(entity_id %in% existing$entity_id & date_used_lin_reg == 'NULL', date_used_linear_regress, date_used_lin_reg),
         date_used_Bacon = if_else(entity_id %in% existing$entity_id & date_used_Bacon.y == 'NULL', date_used_Bacon.x, date_used_Bacon.y),
         date_used_Bchron = if_else(entity_id %in% existing$entity_id & date_used_Bchron.y == 'NULL', date_used_Bchron.x, date_used_Bchron.y)) %>% 
  select(entity_id, dating_id, date_used_lin_reg, date_used_lin_interp, date_used_copRa, date_used_StalAge, date_used_Bacon, date_used_Bchron, date_used_OxCal)

original.f <- original.new %>% filter(entity_id == 691)

write.csv(original,'~/SISAL/date_file_merged.csv', row.names = F)
x <- original %>% distinct(entity_id)



original.new <- original %>% mutate(date_used_Bacon = case_when(entity_id == 163 & dating_id %in% e163.new$labID ~ 'yes',
                                                                entity_id == 163 & !(dating_id %in% e163.new$labID) ~ 'no',
                                                                TRUE ~ date_used_Bacon))
original.f <- original.new %>% filter(entity_id == 163)
write.csv(original.new, '~/SISAL/date_file_merged.csv')
save(sisal_chrono_191222, file = '~/SISAL/sisal_chrono_191222.RData')
save(date_file_191222, file = '~/SISAL/date_file_191222.RData')

save_ens(r_f, '~/SISAL/v2_upd_2/', sisalv2interim$entity)

save_ens <- function(r, fp, entity) {
  eID <- (r %>% distinct(entity_id))$entity_id
  
  for(i in eID) {
    print(i)
    
    entity_name <- (entity %>% filter(entity_id == i))$entity_name
    file_name <- get(load(file.path(fp, paste(i, '-', entity_name, '.RData', sep = ''))))
    
    paste(i, '-', entity_name, '-', sep = '')
    if(!plyr::empty(data.frame(file_name$linReg$ens))){
      ens <- file_name$linReg$ens
      names(ens) <- c('sample_id', 'depth_sample', paste(seq(dim(ens)[2]-2)))
      assign(paste(i, '-', entity_name, '-linReg', sep = ''),ens)
      save(list = paste(i, '-', entity_name, '-linReg', sep = ''), file = file.path(paste('/stacywork/ariana/Kira/',i, '-', entity_name, '-linReg','.RData', sep = '')))
    }
    
    if(!plyr::empty(data.frame(file_name$linInterp$gr))){
      ens <- file_name$linInterp$ens
      names(ens) <- c('sample_id', 'depth_sample', paste(seq(dim(ens)[2]-2)))
      assign(paste(i, '-', entity_name, '-linInterp', sep = ''),ens)
      save(list = paste(i, '-', entity_name, '-linInterp', sep = ''), file = file.path(paste('/stacywork/ariana/Kira/',i, '-', entity_name, '-linInterp','.RData', sep = '')))
    }
    
    if(!plyr::empty(data.frame(file_name$copRa$gr))){
      ens <- file_name$copRa$ens
      names(ens) <- c('sample_id', 'depth_sample', paste(seq(dim(ens)[2]-2)))
      assign(paste(i, '-', entity_name, '-copRa', sep = ''),ens)
      save(list = paste(i, '-', entity_name, '-copRa', sep = ''), file = file.path(paste('/stacywork/ariana/Kira/',i, '-', entity_name, '-copRa','.RData', sep = '')))
    }
    
    if(!plyr::empty(data.frame(file_name$Bacon$gr))){
      ens <- file_name$Bacon$ens
      names(ens) <- c('sample_id', 'depth_sample', paste(seq(dim(ens)[2]-2)))
      assign(paste(i, '-', entity_name, '-Bacon', sep = ''),ens)
      save(list = paste(i, '-', entity_name, '-Bacon', sep = ''), file = file.path(paste('/stacywork/ariana/Kira/',i, '-', entity_name, '-Bacon','.RData', sep = '')))
    }
    
    if(!plyr::empty(data.frame(file_name$Bchron$gr))){
      ens <- file_name$Bchron$ens
      names(ens) <- c('sample_id', 'depth_sample', paste(seq(dim(ens)[2]-2)))
      assign(paste(i, '-', entity_name, '-Bchron', sep = ''),ens)
      save(list = paste(i, '-', entity_name, '-Bchron', sep = ''), file = file.path(paste('/stacywork/ariana/Kira/',i, '-', entity_name, '-Bchron','.RData', sep = '')))
    }
  }
    
}
