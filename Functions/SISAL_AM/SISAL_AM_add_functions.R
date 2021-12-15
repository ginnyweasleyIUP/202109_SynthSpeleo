load_data <- function(prefix, wd) {
  setwd(wd)
  
  composite_link_entity <- read.csv(paste(prefix, 'composite_link_entity.csv',sep = ''), header = T,stringsAsFactors = F)
  d13C <- read.csv(paste(prefix, 'd13C.csv',sep='') ,header = T, stringsAsFactors = F)
  d13C <- rename(d13C, iso_std_d13C = iso_std )
  d18O <- read.csv(paste(prefix, 'd18O.csv', sep =''),header = T, stringsAsFactors = F)
  d18O <- rename(d18O, iso_std_d18O = iso_std)
  dating_lamina <- read.csv(paste(prefix, 'dating_lamina.csv', sep = ''), header = T, stringsAsFactors = F)
  dating <- read.csv(paste(prefix, 'dating.csv',sep = ''), header = T, stringsAsFactors = F)
  entity_link_reference <- read.csv(paste(prefix, 'entity_link_reference.csv', sep = ''), header =T, stringsAsFactors = F)
  entity <- read.csv(paste(prefix, 'entity.csv', sep = ''), header = T, stringsAsFactors = F)
  gap <- read.csv(paste(prefix, 'gap.csv', sep = ''), header = T, stringsAsFactors = F)
  hiatus <- read.csv(paste(prefix, 'hiatus.csv', sep =''), header = T, stringsAsFactors = F)
  notes <- read.csv(paste(prefix, 'notes.csv', sep = ''), header = T, stringsAsFactors = F)
  original_chronology <- read.csv(paste(prefix, 'original_chronology.csv', sep = ''), header = T, stringsAsFactors = F)
  reference <- read.csv(paste(prefix, 'reference.csv', sep = ''), header = T, stringsAsFactors = F)
  sample <- read.csv(paste(prefix, 'sample.csv', sep = ''), header = T, stringsAsFactors = F)
  sisal_chronology <- read.csv(paste(prefix, 'sisal_chronology.csv', sep = ''), header = T, stringsAsFactors = F)
  site <- read.csv(paste(prefix, 'site.csv', sep = ''), header = T, stringsAsFactors = F)
  
  site_tb <- left_join(site, entity, by = 'site_id') %>% left_join(., entity_link_reference, by = 'entity_id') %>% 
    left_join(., reference, by = 'ref_id') %>% left_join(., notes, by = 'site_id') %>% mutate_at(vars(site_id, entity_id), as.numeric)
  dating_tb <- dating %>% group_by(entity_id) %>%mutate(laminar_dated = if_else((entity_id %in% dating_lamina$entity_id), 'yes', 'no')) %>% 
    mutate_at(vars(dating_id, depth_dating, dating_thickness, X14C_correction, corr_age, corr_age_uncert_pos, corr_age_uncert_neg), as.numeric) %>%ungroup()
  sample_tb <- plyr::join_all(list(sample,hiatus, gap, original_chronology, sisal_chronology, d13C, d18O), by = 'sample_id', type = 'left', match = 'all') %>% 
    mutate_at(vars(entity_id, sample_id, sample_thickness, depth_sample, interp_age, interp_age_uncert_pos, interp_age_uncert_neg, COPRA_age,
                  COPRA_age_uncert_pos, COPRA_age_uncert_neg, linear_age, linear_age_uncert_pos, linear_age_uncert_neg, d13C_measurement,
                  d13C_precision, d18O_measurement, d18O_precision), as.numeric)
  
  return(list(site_tb, dating_tb, sample_tb))
}

write_files_old <- function(entid, site_tb, dating_tb, sample_tb, bacon = F, bchron = F, stalage = F, linInterp = F, copRa = F, linReg = F, working_directory, du) {
  site_tb <- as.data.frame(site_tb) 
  dating_tb <- as.data.frame(dating_tb)
  sample_tb <- as.data.frame(sample_tb)
  
  entity_from_base <- site_tb %>% filter(depth_ref == 'from base') %>% distinct(entity_id)
  sample_from_base <- sample_tb %>% filter(entity_id %in% entity_from_base$entity_id) %>% select(entity_id,depth_sample) %>% group_by(entity_id) %>% summarise(max = max(depth_sample))
  
  dating_from_base <- full_join(dating_tb, sample_from_base, by = 'entity_id') %>% group_by(entity_id) %>% 
    mutate(depth_conv = if_else(entity_id %in% entity_from_base$entity_id, max-depth_dating, NA_real_)) %>% 
    mutate(depth_dating = if_else(!is.na(depth_conv), depth_conv, depth_dating)) %>%
    select(-depth_conv) %>% arrange(., depth_dating, .by_group = T)
  
  sampling_from_base <- full_join(sample_tb, sample_from_base, by = 'entity_id') %>% group_by(entity_id) %>% 
    mutate(depth_conv = if_else(entity_id %in% entity_from_base$entity_id, max-depth_sample, NA_real_)) %>% 
    mutate(depth_sample = if_else(!is.na(depth_conv), depth_conv, depth_sample)) %>%
    select(-depth_conv) %>% arrange(., depth_sample, .by_group = T)
  
  dating_tb_new <- dating_from_base %>% group_by(entity_id) %>%arrange(depth_dating) %>% mutate(depth_dating_new = depth_dating/10, 
                                                              corr_age_uncert = (corr_age_uncert_pos + corr_age_uncert_neg)/4,
                                                              thickness_new = if_else(is.na(dating_thickness), 0.5, dating_thickness/20)) %>% ungroup() %>% 
   mutate(calib_curve_new = case_when(
     calib_used == 'NULL' ~ 'normal',
     calib_used == 'unknown' ~ 'unknown',
     calib_used == "not calibrated" ~ 'ask again',
     calib_used == "INTCAL13 NH" ~ 'intcal13'
   ))%>% 
    mutate(cc = case_when(
      date_type == 'C14' ~ 1,
      date_type != 'C14' ~ 0
    )) 
            
 sample_tb_new <- sampling_from_base %>% group_by(entity_id) %>% mutate(depth_sample_new = depth_sample/10) %>% ungroup()
 
 
 entity_name <- unique((site_tb %>% filter(entity_id == entid))$entity_name)
 file_name <- paste(entid,"-",entity_name, sep = '')
 setwd(working_directory)
 dir.create(file.path(working_directory,file_name))
 setwd(file.path(working_directory, file_name))
 #write.csv(unique((site_tb %>% filter(entity_id == entid))$notes), "notes.txt", row.names = F, col.names = F)
 write.csv(site_tb %>% filter(entity_id == entid) %>% select(site_id, site_name, latitude, longitude, entity_id, entity_name, depth_ref, speleothem_type, contact, data_DOI_URL, ref_id, citation, publication_DOI, notes), 'info.csv', row.names = F, col.names = T)
 write.csv(sample_tb_new %>% filter(entity_id == entid & hiatus =='H') %>% select(sample_id, depth_sample) %>% arrange(., depth_sample), 'hiatus.csv', row.names = F, col.names = T)
 write.csv(dating_tb_new %>% filter(entity_id == entid & get(du) == 'yes' & date_type != 'Event; hiatus') %>% select(depth_dating, corr_age, corr_age_uncert_pos, corr_age_uncert_neg, date_type) %>% arrange(., depth_dating), 'used_dates.csv', row.names = F, col.names = T)
 write.csv(dating_tb_new %>% filter(entity_id == entid  & date_type != 'Event; hiatus') %>% filter(get(du) == 'no' | get(du) == 'unknown') %>% select(depth_dating, corr_age) %>% arrange(., depth_dating), 'not_used_dates.csv', row.names = F, col.names = T)
 write.csv(sample_tb_new %>% filter(entity_id == entid) %>% select(sample_id, depth_sample, interp_age, d13C_measurement, d18O_measurement) %>% arrange(., depth_sample), "proxy_data.csv", row.names = F, col.names = T)
 write.csv(sample_tb_new %>% filter(entity_id == entid) %>% select(sample_id, depth_sample, interp_age, interp_age_uncert_pos, interp_age_uncert_neg, age_model_type) %>% arrange(., depth_sample), "original_chronology.csv", row.names = F, col.names = T)
 
 if (bacon) {
   setwd(file.path(working_directory, file_name))
   dir.create('Bacon_runs')
   setwd(file.path(getwd(),'/Bacon_runs'))
   dir.create(file_name)
   setwd(file.path(working_directory, file_name,'Bacon_runs', file_name))
   
   hiatus_bacon <- sample_tb_new %>% filter(entity_id == entid & hiatus =='H') %>% arrange(depth_sample) %>% mutate(depth_sample_bacon = depth_sample/10)  %>% select(sample_id, depth_sample_bacon)
   id <- sample_tb_new %>% filter(entity_id == entid & is.na(hiatus) & is.na(gap)) %>% select(sample_id, depth_sample_new) %>% arrange(., depth_sample_new)
   
   write.csv(id$sample_id, 'sample_id.csv', row.names = F, col.names = F)
   write.table(id$depth_sample_new, paste(entid,'-',entity_name,'_depths.txt', sep = ''), row.names = F, col.names = F)
   write.csv(dating_tb_new %>% filter(entity_id == entid & get(du) == 'yes' & date_type != 'Event; hiatus') %>% select(dating_id, corr_age, corr_age_uncert, depth_dating_new, cc) %>% arrange(., depth_dating_new), paste(entid,'-',entity_name,'.csv', sep = ''), row.names = F, col.names = T)
   write.csv(hiatus_bacon,'hiatus_bacon.csv', row.names = F, col.names = T)
 }
 
 if (bchron){
   setwd(file.path(working_directory, file_name))
   dir.create('Bchron')
   setwd(file.path(getwd(),'/Bchron'))
   
   write.csv(dating_tb_new %>% filter(entity_id == entid & get(du) == 'yes' & date_type != 'Event; hiatus') %>% 
               select(dating_id, corr_age, corr_age_uncert, depth_dating_new, thickness_new, calib_curve_new) %>% arrange(., depth_dating_new), 'ages.csv', row.names = F, col.names = T)
   write.csv(sample_tb_new %>% filter(entity_id == entid) %>% select(sample_id, depth_sample_new) %>% arrange(., depth_sample_new), 'depths.csv', row.names = F, col.names = T)
 }
 
 if (stalage) {
   setwd(file.path(working_directory, file_name))
   dir.create('StalAge')
   setwd(file.path(getwd(),'/StalAge'))
   
   write.csv(dating_tb_new %>% filter(entity_id == entid & get(du) == 'yes' & date_type != 'Event; hiatus') %>% select(dating_id, corr_age, corr_age_uncert, depth_dating) %>% arrange(., depth_dating), 'ages.csv', row.names = F, col.names = T)
   write.csv(sample_tb_new %>% filter(entity_id == entid) %>% select(sample_id, depth_sample) %>% arrange(., depth_sample), 'depths.csv', row.names = F, col.names = T)
 }
 
 if (linInterp) {
   setwd(file.path(working_directory, file_name))
   dir.create('linInterp')
   setwd(file.path(getwd(),'/linInterp'))
   
   write.csv(sample_tb_new %>% filter(entity_id == entid) %>% select(sample_id, depth_sample) %>% arrange(., depth_sample), 'depths.csv', row.names = F, col.names = T)
   write.csv(dating_tb_new %>% filter(entity_id == entid & get(du) == 'yes' & date_type != 'Event; hiatus') %>% select(dating_id, corr_age, corr_age_uncert, depth_dating, date_type) %>% 
               arrange(., depth_dating),'ages.csv', row.names = F, col.names = T)
 }
 
 if (copRa) {
   setwd(file.path(working_directory, file_name))
   dir.create('copRa')
   setwd(file.path(getwd(),'/copRa'))
   
   write.csv(sample_tb_new %>% filter(entity_id == entid) %>% select(sample_id, depth_sample) %>% arrange(., depth_sample), 'depths.csv', row.names = F, col.names = T)
   write.csv(dating_tb_new %>% filter(entity_id == entid & get(du) == 'yes' & date_type != 'Event; hiatus') %>% select(dating_id, corr_age, corr_age_uncert, depth_dating, date_type) %>% 
               arrange(., depth_dating),'ages.csv', row.names = F, col.names = T)
 }
 
 if (linReg) {
   setwd(file.path(working_directory, file_name))
   dir.create('linReg')
   setwd(file.path(getwd(),'/linReg'))
   
   write.csv(sample_tb_new %>% filter(entity_id == entid & is.na(hiatus) & is.na(gap)) %>% select(sample_id, depth_sample) %>% arrange(., depth_sample), 'depths.csv', row.names = F, col.names = T)
   write.csv(sample_tb_new %>% filter(entity_id == entid) %>% select(sample_id, depth_sample) %>% arrange(., depth_sample), 'id.csv', row.names = F, col.names = T)
   write.csv(dating_tb_new %>% filter(entity_id == entid & get(du) == 'yes' & date_type != 'Event; hiatus') %>% select(dating_id, corr_age, corr_age_uncert, depth_dating, date_type) %>% 
               arrange(., depth_dating),'ages.csv', row.names = F, col.names = T)
 }
 
 return(file_name)
} 

add_hiatus <- function(data, hiatus, stalage = F, bchron =F, lininterp=F) {
  
  age <- unlist(data[,2])
  depth_dating <- unlist(data[,4])
  #uncert <- unlist(data[,3])
  x_out <- unlist(hiatus$depth_sample)
  
  e <- approx(x = depth_dating, y = age, xout = x_out, method = 'linear')
  
  if(lininterp){
    e <- approx(x = depth_dating, y = age, xout = x_out, method = 'linear')
    h <- data.frame(dating_id = hiatus$sample_id,corr_age = e$y, corr_age_uncert = rep(NA, length(hiatus$depth_sample)),depth_dating = e$x, date_type = rep('Hiatus',length(hiatus$depth_sample)))
    new <- rbind(data, h)
    new <- new %>% arrange(., corr_age) %>% mutate(corr_age_uncert = if_else(is.na(corr_age_uncert), ((lead(corr_age_uncert)+lead(corr_age)-corr_age)+(corr_age-(lag(corr_age)-lag(corr_age_uncert))))/2, as.double(corr_age_uncert)))
  }
  
  if(stalage){
    e <- approx(x = depth_dating, y = age, xout = x_out, method = 'linear')
    h <- data.frame(dating_id = hiatus$sample_id,corr_age = e$y, corr_age_uncert = rep(NA, length(hiatus$depth_sample)),depth_dating = e$x)
    new <- rbind(data, h)
    new <- new %>% arrange(., corr_age) %>% mutate(corr_age_uncert = if_else(is.na(corr_age_uncert), ((lead(corr_age_uncert)+lead(corr_age)-corr_age)+(corr_age-(lag(corr_age)-lag(corr_age_uncert))))/2, as.double(corr_age_uncert)))
  }
  
  if(bchron){
    x_out <- x_out/10
    e <- approx(x = depth_dating, y = age, xout = x_out, method = 'linear')
    h <- data.frame(dating_id = hiatus$sample_id,corr_age = e$y, corr_age_uncert = rep(NA, length(hiatus$depth_sample)),depth_dating_new = e$x, thickness_new = rep(NA, length(hiatus$depth_sample)), calib_curve_new = 'normal')
    new <- rbind(data, h)
    new <- new %>% arrange(., corr_age) %>% mutate(corr_age_uncert = if_else(is.na(corr_age_uncert), ((lead(corr_age_uncert)+lead(corr_age)-corr_age)+(corr_age-(lag(corr_age)-lag(corr_age_uncert))))/2, as.double(corr_age_uncert)),
                                                   thickness_new = if_else(is.na(thickness_new), 0.01, as.double(thickness_new)),
                                                   #thickness_new = if_else(is.na(thickness_new), (lead(thickness_new)+lag(thickness_new))/2, as.double(thickness_new)),
                                                   calib_curve_new = if_else(is.na(calib_curve_new), 'normal', calib_curve_new))
                                                  
  }
  
  return(new) 
}

write_files <- function(entid, site_tb, dating_tb, sample_tb, bacon = F, bchron = F, stalage = F, linInterp = F, copRa = F, linReg = F, working_directory) {
  site_tb <- as.data.frame(site_tb) 
  dating_tb <- as.data.frame(dating_tb)
  sample_tb <- as.data.frame(sample_tb)
  
  entity_from_base <- site_tb %>% filter(depth_ref == 'from base') %>% distinct(entity_id)
  sample_from_base <- sample_tb %>% filter(entity_id %in% entity_from_base$entity_id) %>% select(entity_id,depth_sample) %>% group_by(entity_id) %>% summarise(max = max(depth_sample))
  
  dating_from_base <- full_join(dating_tb, sample_from_base, by = 'entity_id') %>% group_by(entity_id) %>% 
    mutate(depth_conv = if_else(entity_id %in% entity_from_base$entity_id, max-depth_dating, NA_real_)) %>% 
    mutate(depth_dating = if_else(!is.na(depth_conv), depth_conv, depth_dating)) %>%
    select(-depth_conv) %>% arrange(., depth_dating, .by_group = T)
  
  sampling_from_base <- full_join(sample_tb, sample_from_base, by = 'entity_id') %>% group_by(entity_id) %>% 
    mutate(depth_conv = if_else(entity_id %in% entity_from_base$entity_id, max-depth_sample, NA_real_)) %>% 
    mutate(depth_sample = if_else(!is.na(depth_conv), depth_conv, depth_sample)) %>%
    select(-depth_conv) %>% arrange(., depth_sample, .by_group = T)
  
  dating_tb_new <- dating_from_base %>% group_by(entity_id) %>%arrange(depth_dating) %>% mutate(depth_dating_new = depth_dating/10, 
                                                                                                corr_age_uncert = (corr_age_uncert_pos + corr_age_uncert_neg)/4,
                                                                                                thickness_new = if_else(is.na(dating_thickness), 0.5, dating_thickness/20)) %>% ungroup() %>% 
    mutate(calib_curve_new = case_when(
      calib_used == 'NULL' ~ 'normal',
      calib_used == 'unknown' ~ 'unknown',
      calib_used == "not calibrated" ~ 'ask again',
      calib_used == "INTCAL13 NH" ~ 'intcal13'
    ))%>% 
    mutate(cc = case_when(
      date_type == 'C14' ~ 1,
      date_type != 'C14' ~ 0
    )) 
  
  sample_tb_new <- sampling_from_base %>% group_by(entity_id) %>% mutate(depth_sample_new = depth_sample/10) %>% ungroup()
  
  
  entity_name <- unique((site_tb %>% filter(entity_id == entid))$entity_name)
  file_name <- paste(entid,"-",entity_name, sep = '')
  setwd(working_directory)
  dir.create(file.path(working_directory,file_name))
  setwd(file.path(working_directory, file_name))
  #write.csv(unique((site_tb %>% filter(entity_id == entid))$notes), "notes.txt", row.names = F, col.names = F)
  write.csv(site_tb %>% filter(entity_id == entid) %>% select(site_id, site_name, latitude, longitude, entity_id, entity_name, depth_ref, speleothem_type, contact, data_DOI_URL, ref_id, citation, publication_DOI, notes), 'info.csv', row.names = F, col.names = T)
  write.csv(sample_tb_new %>% filter(entity_id == entid & hiatus =='H') %>% select(sample_id, depth_sample) %>% arrange(., depth_sample), 'hiatus.csv', row.names = F, col.names = T)
  write.csv(dating_tb_new %>% filter(entity_id == entid & date_used == 'yes' & date_type != 'Event; hiatus') %>% select(depth_dating, corr_age, corr_age_uncert_pos, corr_age_uncert_neg, date_type) %>% arrange(., depth_dating), 'used_dates.csv', row.names = F, col.names = T)
  write.csv(dating_tb_new %>% filter(entity_id == entid  & date_type != 'Event; hiatus') %>% filter(date_used == 'no' | date_used == 'unknown') %>% select(depth_dating, corr_age) %>% arrange(., depth_dating), 'not_used_dates.csv', row.names = F, col.names = T)
  write.csv(sample_tb_new %>% filter(entity_id == entid) %>% select(sample_id, depth_sample, interp_age, d13C_measurement, d18O_measurement) %>% arrange(., depth_sample), "proxy_data.csv", row.names = F, col.names = T)
  write.csv(sample_tb_new %>% filter(entity_id == entid) %>% select(sample_id, depth_sample, interp_age, interp_age_uncert_pos, interp_age_uncert_neg, age_model_type) %>% arrange(., depth_sample), "original_chronology.csv", row.names = F, col.names = T)
  
  if (bacon) {
    setwd(file.path(working_directory, file_name))
    dir.create('Bacon_runs')
    setwd(file.path(getwd(),'/Bacon_runs'))
    dir.create(file_name)
    setwd(file.path(working_directory, file_name,'Bacon_runs', file_name))
    
    hiatus_bacon <- sample_tb_new %>% filter(entity_id == entid & hiatus =='H') %>% arrange(depth_sample) %>% mutate(depth_sample_bacon = depth_sample/10)  %>% select(sample_id, depth_sample_bacon)
    id <- sample_tb_new %>% filter(entity_id == entid & is.na(hiatus) & is.na(gap)) %>% select(sample_id, depth_sample_new) %>% arrange(., depth_sample_new)
    
    write.csv(id$sample_id, 'sample_id.csv', row.names = F, col.names = F)
    write.table(id$depth_sample_new, paste(entid,'-',entity_name,'_depths.txt', sep = ''), row.names = F, col.names = F)
    write.csv(dating_tb_new %>% filter(entity_id == entid & date_used == 'yes' & date_type != 'Event; hiatus') %>% select(dating_id, corr_age, corr_age_uncert, depth_dating_new, cc) %>% arrange(., depth_dating_new), paste(entid,'-',entity_name,'.csv', sep = ''), row.names = F, col.names = T)
    write.csv(hiatus_bacon,'hiatus_bacon.csv', row.names = F, col.names = T)
  }
  
  if (bchron){
    setwd(file.path(working_directory, file_name))
    dir.create('Bchron')
    setwd(file.path(getwd(),'/Bchron'))
    
    write.csv(dating_tb_new %>% filter(entity_id == entid & date_used == 'yes' & date_type != 'Event; hiatus') %>% 
                select(dating_id, corr_age, corr_age_uncert, depth_dating_new, thickness_new, calib_curve_new) %>% arrange(., depth_dating_new), 'ages.csv', row.names = F, col.names = T)
    write.csv(sample_tb_new %>% filter(entity_id == entid) %>% select(sample_id, depth_sample_new) %>% arrange(., depth_sample_new), 'depths.csv', row.names = F, col.names = T)
  }
  
  if (stalage) {
    setwd(file.path(working_directory, file_name))
    dir.create('StalAge')
    setwd(file.path(getwd(),'/StalAge'))
    
    write.csv(dating_tb_new %>% filter(entity_id == entid & date_used == 'yes' & date_type != 'Event; hiatus') %>% select(dating_id, corr_age, corr_age_uncert, depth_dating) %>% arrange(., depth_dating), 'ages.csv', row.names = F, col.names = T)
    write.csv(sample_tb_new %>% filter(entity_id == entid) %>% select(sample_id, depth_sample) %>% arrange(., depth_sample), 'depths.csv', row.names = F, col.names = T)
  }
  
  if (linInterp) {
    setwd(file.path(working_directory, file_name))
    dir.create('linInterp')
    setwd(file.path(getwd(),'/linInterp'))
    
    write.csv(sample_tb_new %>% filter(entity_id == entid) %>% select(sample_id, depth_sample) %>% arrange(., depth_sample), 'depths.csv', row.names = F, col.names = T)
    write.csv(dating_tb_new %>% filter(entity_id == entid & date_used == 'yes' & date_type != 'Event; hiatus') %>% select(dating_id, corr_age, corr_age_uncert, depth_dating, date_type) %>% 
                arrange(., depth_dating),'ages.csv', row.names = F, col.names = T)
  }
  
  if (copRa) {
    setwd(file.path(working_directory, file_name))
    dir.create('copRa')
    setwd(file.path(getwd(),'/copRa'))
    
    write.csv(sample_tb_new %>% filter(entity_id == entid) %>% select(sample_id, depth_sample) %>% arrange(., depth_sample), 'depths.csv', row.names = F, col.names = T)
    write.csv(dating_tb_new %>% filter(entity_id == entid & date_used == 'yes' & date_type != 'Event; hiatus') %>% select(dating_id, corr_age, corr_age_uncert, depth_dating, date_type) %>% 
                arrange(., depth_dating),'ages.csv', row.names = F, col.names = T)
  }
  
  if (linReg) {
    setwd(file.path(working_directory, file_name))
    dir.create('linReg')
    setwd(file.path(getwd(),'/linReg'))
    
    write.csv(sample_tb_new %>% filter(entity_id == entid & is.na(hiatus) & is.na(gap)) %>% select(sample_id, depth_sample) %>% arrange(., depth_sample), 'depths.csv', row.names = F, col.names = T)
    write.csv(sample_tb_new %>% filter(entity_id == entid) %>% select(sample_id, depth_sample) %>% arrange(., depth_sample), 'id.csv', row.names = F, col.names = T)
    write.csv(dating_tb_new %>% filter(entity_id == entid & date_used == 'yes' & date_type != 'Event; hiatus') %>% select(dating_id, corr_age, corr_age_uncert, depth_dating, date_type) %>% 
                arrange(., depth_dating),'ages.csv', row.names = F, col.names = T)
  }
  
  return(file_name)
} 

add_hiatus <- function(data, hiatus, stalage = F, bchron =F, lininterp=F) {
  
  age <- unlist(data[,2])
  depth_dating <- unlist(data[,4])
  #uncert <- unlist(data[,3])
  x_out <- unlist(hiatus$depth_sample)
  
  e <- approx(x = depth_dating, y = age, xout = x_out, method = 'linear')
  
  if(lininterp){
    e <- approx(x = depth_dating, y = age, xout = x_out, method = 'linear')
    h <- data.frame(dating_id = hiatus$sample_id,corr_age = e$y, corr_age_uncert = rep(NA, length(hiatus$depth_sample)),depth_dating = e$x, date_type = rep('Hiatus',length(hiatus$depth_sample)))
    new <- rbind(data, h)
    new <- new %>% arrange(., corr_age) %>% mutate(corr_age_uncert = if_else(is.na(corr_age_uncert), ((lead(corr_age_uncert)+lead(corr_age)-corr_age)+(corr_age-(lag(corr_age)-lag(corr_age_uncert))))/2, as.double(corr_age_uncert)))
  }
  
  if(stalage){
    e <- approx(x = depth_dating, y = age, xout = x_out, method = 'linear')
    h <- data.frame(dating_id = hiatus$sample_id,corr_age = e$y, corr_age_uncert = rep(NA, length(hiatus$depth_sample)),depth_dating = e$x)
    new <- rbind(data, h)
    new <- new %>% arrange(., corr_age) %>% mutate(corr_age_uncert = if_else(is.na(corr_age_uncert), ((lead(corr_age_uncert)+lead(corr_age)-corr_age)+(corr_age-(lag(corr_age)-lag(corr_age_uncert))))/2, as.double(corr_age_uncert)))
  }
  
  if(bchron){
    x_out <- x_out/10
    e <- approx(x = depth_dating, y = age, xout = x_out, method = 'linear')
    h <- data.frame(dating_id = hiatus$sample_id,corr_age = e$y, corr_age_uncert = rep(NA, length(hiatus$depth_sample)),depth_dating_new = e$x, thickness_new = rep(NA, length(hiatus$depth_sample)), calib_curve_new = 'normal')
    new <- rbind(data, h)
    new <- new %>% arrange(., corr_age) %>% mutate(corr_age_uncert = if_else(is.na(corr_age_uncert), ((lead(corr_age_uncert)+lead(corr_age)-corr_age)+(corr_age-(lag(corr_age)-lag(corr_age_uncert))))/2, as.double(corr_age_uncert)),
                                                   thickness_new = if_else(is.na(thickness_new), 0.01, as.double(thickness_new)),
                                                   #thickness_new = if_else(is.na(thickness_new), (lead(thickness_new)+lag(thickness_new))/2, as.double(thickness_new)),
                                                   calib_curve_new = if_else(is.na(calib_curve_new), 'normal', calib_curve_new))
    
  }
  
  return(new) 
}


get_bacon_median_quantile <- function(depth_eval, hiatus, bacon_mcmc) {
  #bacon_mcmc <- sapply(depth_eval, Bacon.Age.d)
  bacon_age <- apply(bacon_mcmc,2,median)
  bacon_quantile <- apply(bacon_mcmc, 2, function(x){quantile(x, probs = c(0.05,0.95), na.rm = T)})
  
  data <- cbind(depth_eval, bacon_age, bacon_age_uncert_pos = bacon_quantile[2,]-bacon_age, 
                bacon_age_uncert_neg = bacon_age - bacon_quantile[1,])
  h <- data.frame(depth_eval = hiatus$depth_sample_bacon, bacon_age = replicate(dim(hiatus)[1],NA), 
                  bacon_age_uncert_pos = replicate(dim(hiatus)[1],NA),
                  bacon_age_uncert_neg = replicate(dim(hiatus)[1],NA))
  data <- rbind(data, h)
  data <- data[order(data[,1]),] 
  
  return(data)
}

get_bacon_median_quantile_new <- function(depth_eval, hiatus, bacon_mcmc) {
  #bacon_mcmc <- sapply(depth_eval, Bacon.Age.d)
  bacon_age <- apply(bacon_mcmc,1,median)
  bacon_quantile <- apply(bacon_mcmc, 1, function(x){quantile(x, probs = c(0.025,0.975), na.rm = T)})
  
  data <- cbind(depth_eval, bacon_age, bacon_age_uncert_pos = bacon_quantile[2,]-bacon_age, 
                bacon_age_uncert_neg = bacon_age - bacon_quantile[1,])
  h <- data.frame(depth_eval = hiatus$depth_sample_bacon, bacon_age = replicate(dim(hiatus)[1],NA), 
                  bacon_age_uncert_pos = replicate(dim(hiatus)[1],NA),
                  bacon_age_uncert_neg = replicate(dim(hiatus)[1],NA))
  data <- rbind(data, h)
  data <- data[order(data[,1]),] 
  
  return(data)
}

get_lin_interp_new <- function(data, depth_eval, method, hiatus) {
  
  library(plyr)
  library(rlist)
  library(Hmisc)
  
  #if (!empty(data.frame(hiatus))) {
  #  data <- add_hiatus(data, hiatus)
  #}
  
  age <- unlist(data[,2])
  depth_dating <- unlist(data[,1])
  x_out <- unlist(depth_eval)
  
  e <- approxExtrap(x = depth_dating, y = age, xout = x_out)
  linInterp <- as_tibble(data.frame(depth_eval = unlist(e$x), lin_interp_age = unlist(e$y)))
  
  if (!empty(data.frame(hiatus))) {
    linInterp <- linInterp %>% rowwise() %>% mutate(lin_interp_age = if_else(depth_eval %in% hiatus, NA_real_, lin_interp_age))
  }
  
  return(linInterp)
}

get_copRa <- function(data, depth_eval, method, hiatus) {
  
  library(plyr)
  library(rlist)
  library(Hmisc)
  library(pracma)
  
  #if (!empty(data.frame(hiatus))) {
  #  data <- add_hiatus(data, hiatus)
  #}
  
  age <- unlist(data[,2])
  depth_dating <- unlist(data[,1])
  x_out <- unlist(depth_eval)
  x_out_1 <- x_out[which(x_out <= max(depth_dating) & x_out >= min(depth_dating))]
  x_out_2 <- x_out[which(x_out > max(depth_dating) | x_out < min(depth_dating))]
  
  e <- pchip(xi = depth_dating, yi = age, x = x_out_1)
  copRa <- as_tibble(data.frame(depth_eval = x_out_1, copRa_age = e))
  
  if(!empty(data.frame(x_out_2))){
    e_2 <- approxExtrap(x = depth_dating, y = age, xout = x_out_2)
    copRa_2 <- as_tibble(data.frame(depth_eval = unlist(e_2$x), copRa_age = unlist(e_2$y)))
    copRa <- rbind(copRa, copRa_2)
    copRa <- copRa %>% arrange(., depth_eval)
  }
  
  if (!empty(data.frame(hiatus))) {
    copRa <- copRa %>% rowwise() %>% mutate(copRa_age = if_else(depth_eval %in% hiatus, NA_real_, copRa_age))
  }
  
  return(copRa)
}

get_copRa_falsch <- function(data, depth_eval, method, hiatus) {
  
  library(plyr)
  library(rlist)
  library(Hmisc)
  
  #if (!empty(data.frame(hiatus))) {
  #  data <- add_hiatus(data, hiatus)
  #}
  
  age <- unlist(data[,2])
  depth_dating <- unlist(data[,1])
  x_out <- unlist(depth_eval)
  
  e <- spline(x = depth_dating, y = age, xout = x_out_1, method = method)
  copRa <- as_tibble(data.frame(depth_eval = unlist(e$x), copRa_age = unlist(e$y)))
  
  if (!empty(data.frame(hiatus))) {
    copRa <- copRa %>% rowwise() %>% mutate(copRa_age = if_else(depth_eval %in% hiatus, NA_real_, copRa_age))
  }
  
  return(copRa)
}

get_lin_interp <- function(data, depth_eval, method, hiatuses) {
  
  library(plyr)
  library(rlist)
  
  age <- data[,2]
  depth_dating <- data[,1]
  
  if (empty(data.frame(hiatuses))) {
    x_out <- depth_eval
    x_e <- depth_dating
    if (length(x_e)<2) {
      #lI <- cbind(x_out)
      lI <- cbind(depth_dating, age)
    } else {
      if (empty(data.frame(x_out))) {
        e <- approx(x = depth_dating, y = age, method = method, rule = 2)
      } else {
        e <- approx(x = depth_dating, y = age, xout = depth_eval, method = method, rule = 2)
      }
      lI <- cbind(unlist(e[[1]]), unlist(e[[2]]))
    }
    
  } else {
    
    j <- length(hiatuses)
    x_e <- depth_dating[depth_dating < hiatuses[1]]
    x_out <- depth_eval[depth_eval <= x_e[length(x_e)]]
    if(length(x_e)==1){
      x_e_new <- c(depth_dating[1], x_e, min(depth_dating[depth_dating > hiatuses[1]])) # changed
      x_out <- depth_eval[depth_eval < x_e_new[length(x_e_new)] & depth_eval > x_e_new[1]]
    }
    
    if (length(x_e)<2) {
      lI <- cbind(x_out, rep(NA, length(x_out)))
      #lI <- c(x_e, age[depth_dating < hiatuses[1]])
    } else {
      
      if (empty(data.frame(x_out))) {
        e <- approx(x = x_e, y = age[depth_dating < hiatuses[1]], method = method, rule = 2)
      } else {
        e <- approx(x = x_e, y = age[depth_dating < hiatuses[1]], xout = x_out, method = method, rule = 2)
      }
      lI <- cbind(unlist(e[[1]]), unlist(e[[2]]))
      
    }
    
    if (j>2) {
      #print(j)
      for (i in seq(2,j)){
        #print(i)
        x_e <- depth_dating[depth_dating < hiatuses[i] & depth_dating > hiatuses[i-1]]
        x_out <- depth_eval[depth_eval < x_e[length(x_e)] & depth_eval > x_e[1]]
        if(length(x_e)==1){
          x_e_new <- c(max(depth_dating[depth_dating < hiatuses[i-1]]), x_e, min(depth_dating[depth_dating > hiatuses[i]]))
          x_out <- depth_eval[depth_eval < x_e_new[length(x_e_new)] & depth_eval > x_e_new[1]]
        }
        #print('>2')
        
        if (length(x_e)<2) {
          lI <- rbind(lI,cbind(x_out, rep(NA, length(x_out))))
          #lI <- rbind(lI,cbind(x_e, age[depth_dating < hiatuses[i]]))
        } else {
          if (empty(data.frame(x_out))) {
            e <- approx(x = x_e, y = age[depth_dating < hiatuses[i] & depth_dating > hiatuses[i-1]], method = method, rule = 2)
          } else {
            e <- approx(x = x_e, y = age[depth_dating < hiatuses[i] & depth_dating > hiatuses[i-1]], xout = x_out, method = method, rule = 2)
          }
          
          lI <- rbind(lI,cbind(unlist(e[[1]]), unlist(e[[2]])))
          
        }
      }
    }
    
    x_e <- depth_dating[depth_dating > hiatuses[j]]
    x_out <- depth_eval[depth_eval >= x_e[1]]
    if(length(x_e)==1){
      x_e_new <- c(max(depth_dating[depth_dating < hiatuses[i-1]]), x_e, min(depth_dating[depth_dating > hiatuses[i]]))
      x_out <- depth_eval[depth_eval < x_e_new[length(x_e_new)] & depth_eval > x_e_new[1]]
    }
    
    if (length(x_e)<2) {
      lI <- rbind(lI,cbind(x_out, rep(NA, length(x_out))))
      #lI <- rbind(lI,cbind(x_e, age[depth_dating > hiatuses[j]]))
    } else {
      if (empty(data.frame(x_out))) {
        e <- approx(x = x_e, y = age[depth_dating > hiatuses[j]], method = method, rule = 2)
      } else {
        e <- approx(x = x_e, y = age[depth_dating > hiatuses[j]], xout = x_out, method = method, rule = 2)
      }
      
      lI <- rbind(lI,cbind(unlist(e[[1]]), unlist(e[[2]])))
    }
    
  }
  
  data <- data.frame(depth_eval = lI[,1], lin_interp_age = lI[,2])
  h <- data.frame(depth_eval = hiatuses, lin_interp_age = replicate(length(hiatuses),NA))
  data <- rbind(data, h)
  d <- data.frame(depth_eval = depth_eval)
  data <- merge(d,data, by.y = 'depth_eval', all.x = T, all.y = T)
  ### new ----
  data <- as.data.frame(lapply(data, unlist))
  ### --------
  data <- data[order(data[,1]),] 
  
  return(data)
}


linear_regression <- function(data, hiatus){ # data = c("depth","age")
  
  # initialize
  j <- length(hiatus)
  m <- list(list())
  
  # calculate slope and intercept for each section between hiati
  
  idx <- data[,1] < hiatus[1]
  m[[1]] <- lm(data[idx,2]~data[idx,1])
  if(is.null(dim(data[idx,]))) {m[[1]] <- NA}
  
  if (j>=2) {
    for (i in seq(from = 2, to = j)){
      idx <- (data[,1] < hiatus[i]) & (data[,1] > hiatus[i-1])
      m[[i]] <- lm(data[idx,2]~data[idx,1])
      if(is.null(dim(data[idx,]))) {m[[i]] <- NA}
    }
    
    idx <- data[,1] > hiatus[length(hiatus)]
    m[[j+1]] <- lm(data[idx,2]~data[idx,1])
    if(is.null(dim(data[idx,]))) {m[[j+1]] <- NA}
    
  } else {
    idx <- data[,1] > hiatus[length(hiatus)]
    m[[j+1]] <- lm(data[idx,2]~data[idx,1])
    if(is.null(dim(data[idx,]))) {m[[j+1]] <- NA}
  }
  
  # retrun data
  return(m)
}

# depths = c("depth in cm")
linear_regression_ages <- function(m, depth_eval, hiatus) {
  
  # initialize
  d <- length(unlist(depth_eval))
  j <- length(m)
  lin_reg_age <- replicate(d, 0)
  
  for (i in seq(1,j)) {
    if(!is.na(m[[i]])) {
      for (k in c(1,2)) {
        if (is.na(m[[i]]$coefficients[[k]])) {
          m[[i]]$coefficients[[k]] = 0
        }
      }
    }
  }
  
  # add column with ages to depth table 
  depth <- cbind(depth_eval, lin_reg_age)
  
  idx <- depth[,1] < hiatus[1]
  if(is.na(m[[1]])){
    depth[idx,2] <- NA
  } else {
    depth[idx,2]<- m[[1]]$coefficients[[1]] + depth[idx,1]*m[[1]]$coefficients[[2]]
  }
  
  
  if (j>2) {
    for (i in seq(2,length(hiatus))){
      idx <- (depth[,1] < hiatus[i]) & (depth[,1] > hiatus[i-1])
      #depth[idx,2]<- m[[i]]$coefficients[[1]] + depth[idx,1]*m[[i]]$coefficients[[2]]
      if(is.na(m[[i]])){
        depth[idx,2] <- NA
      } else {
        depth[idx,2]<- m[[i]]$coefficients[[1]] + depth[idx,1]*m[[i]]$coefficients[[2]]
      }
    }
    
    idx <- depth[,1] > hiatus[length(hiatus)]
    #depth[idx,2]<- m[[length(hiatus)+1]]$coefficients[[1]] + depth[idx,1]*m[[length(hiatus)+1]]$coefficients[[2]]
    if(is.na(m[[1]])){
      depth[idx,2] <- NA
    } else {
      depth[idx,2]<- m[[length(hiatus)+1]]$coefficients[[1]] + depth[idx,1]*m[[length(hiatus)+1]]$coefficients[[2]]
    }
  } else {
    idx <- depth[,1] > hiatus[length(hiatus)]
    #depth[idx,2]<- m[[length(hiatus)+1]]$coefficients[[1]] + depth[idx,1]*m[[length(hiatus)+1]]$coefficients[[2]]
    if(is.na(m[[1]])){
      depth[idx,2] <- NA
    } else {
      depth[idx,2]<- m[[length(hiatus)+1]]$coefficients[[1]] + depth[idx,1]*m[[length(hiatus)+1]]$coefficients[[2]]
    }
  }
  data <- data.frame(depth_eval = depth[,1], lin_reg_age = depth[,2])
  #data <- data[-which(depth_eval == hiatus),]
  h <- data.frame(depth_eval = hiatus, lin_reg_age = replicate(length(hiatus),NA))
  data <- rbind(data, h)
  data <- data[order(data[,1]),] 
  
  return(data)
  
}

lin_reg_no_hiatus <- function(data, depth_eval) {
  m_lR<- lm(data[,2]~data[,1])
  
  age_lR <- replicate(length(unlist(depth_eval)), 0)
  depth <- cbind(depth_eval, age_lR)
  depth[,2]<- m_lR[1]$coefficients[[1]] + depth[,1]*m_lR[1]$coefficients[[2]]
  return(depth)
  
}

plot_AM <- function(age, unknown, used, legend, hiatuses, date_type) {
  j <- ncol(age)
  age[,1] <- age[,1]*10
  unknown[,2] <- unknown[,2]*10
  used[,2] <- used[,2]*10
  matplot(age[,seq(from =2, to = j)], age[,1], type = 'l', lty = seq(2,j), col = seq(2,j), ylim = c(max(age[,1]),0),
          xlab = 'Age [yrs BP]', ylab = 'Depth from top [mm]', lwd = 2)
  points(used[,1], used[,2], pch = 4, col = "orange")
  
  if(length(hiatuses) != 0) {
    abline(h = hiatuses*10, col = 'brown', lty = 6, lwd = 2)
    legend('bottomleft', legend= c(legend,'Hiatus'), lwd = 2,col = c(seq(2,j),'brown'), lty = c(seq(2,j),6), box.lwd = 0,box.col = "white",bg = "white", title = "Additional AM information", cex = 0.65)
  } else {
    legend('bottomleft', legend= legend, col = seq(2,j), lwd = 2,lty = seq(2,j), box.lwd = 0,box.col = "white",bg = "white", title = "Additional AM information", cex = 0.65)
  }
  
  if (length(unknown) != 0) {
    points(unknown[,1],unknown[,2], pch = 25, col ="black")
    legend("topright",legend = c(paste(date_type,"- used"),"U/Th ages - not used"), pch = c(4,25),col = c("orange", "black") ,box.lwd = 0,box.col = "white",bg = "white", title = "AM date types", cex = 0.65)
  } else {
    legend("topright",legend = c(paste(date_type,"- used")), pch = c(4),col = c("orange") ,box.lwd = 0,box.col = "white",bg = "white", title = "AM date types", cex = 0.65)
  }
}

reversal_filter <- function(age, filter){
  d <- diff(age)
  r <- which(d<0)
  if(filter && any(diff(r)<=2)){r <- r[-(which(diff(r)<=2)+1)]}  # UEBERGANGSLOESUNG 
  idx_rev <- matrix(cbind(r,r+1), nrow = length(r), ncol = 2)
  n <- 2
  while (any(apply(matrix(d[idx_rev], nrow = length(r), ncol = n), 1, sum, na.rm = TRUE)<0)) {
    m <- matrix(d[idx_rev], nrow = length(r), ncol = n)
    idx_rev <- cbind(idx_rev,rep(NA,length(r)))
    idx <- which(apply(m, 1, sum, na.rm = TRUE)<0)
    
    idx_rev[idx,dim(idx_rev)[2]] <- r[idx]+n
    n <- n+1 
    if(n > 7) {stop('ERROR: too many reversal! Delete outliers.')}
  }
  
  return(idx_rev)
}

get_median_quantiles <- function(upd){
  age_median <- apply(upd, 1, median)
  age_sd <- apply(upd, 1, function(x){quantile(x, probs = c(0.05,0.95), na.rm = T)})
  
  return(cbind(age_median, t(age_sd)))
}

get_median_quantiles_copRa <- function(upd, q1, q2){
  age_median <- apply(upd, 1, median)
  age_sd <- apply(upd, 1, function(x){quantile(x, probs = c(q1,q2), na.rm = T)})
  
  return(cbind(age_median, t(age_sd)))
}

get_depth_uncertainty <- function(median, sd, depth_dating, depth_sample) {
  #data <- data.frame(depth_dating, median)
  lI_median <- get_lin_interp(cbind(depth_dating, median), depth_sample, 'linear', hiatus)
  lI_sd <- get_lin_interp(cbind(depth_dating, median+sd), depth_sample, 'linear', hiatus)
  #lI_sd_neg <- get_lin_interp(cbind(depth_dating, median-sd), depth_sample, 'linear', hiatus)
  
  lin_interp_mc <- merge(lI_median, lI_sd, by = 'depth_eval')
  lin_interp_mc[,3] <- lin_interp_mc[,3]-lin_interp_mc[,2]
  lI_mc <<- data.frame(depth_eval = lin_interp_mc[,1], median_age = lin_interp_mc[,2], sd_age = lin_interp_mc[,3])
}


mc_linReg <- function(N, hiatus, depth_dating, age_ensemble, depth_sample) {
  #sample <- rep(NA, N)
  for (i in 1:N){
    if (i == 1) {
      if (empty(data.frame(hiatus))){
        age_mc <- lin_reg_no_hiatus(cbind(depth_dating, age_ensemble[i,]), depth_sample)
      }else{
        m <- linear_regression(cbind(depth_dating, age_ensemble[i,]), hiatus) # hiatus = hiatus[[1]]
        age_mc <- linear_regression_ages(m, depth_sample, hiatus)
      }
      sample_ensemble <- age_mc
    } else {
      if (empty(data.frame(hiatus))){
        age_mc <- lin_reg_no_hiatus(cbind(depth_dating, age_ensemble[i,]), depth_sample)
      } else {
        m <- linear_regression(cbind(depth_dating, age_ensemble[i,]), hiatus) 
        age_mc <- linear_regression_ages(m, depth_sample, hiatus)
      }
      sample_ensemble <- cbind(sample_ensemble,age_mc[,2])
    }
  }
  return(sample_ensemble)
}

mc_linInt <- function(N, hiatus, depth_dating, age_ensemble, depth_sample){
  for(j in 1:N){
    if(j == 1){
      age_mc <- get_lin_interp_new(cbind(depth_dating, age_ensemble[j,]), depth_sample, 'linear',hiatus)
      sample_ensemble <- age_mc
      #print(dim(age))
    } else {
      age_mc <- get_lin_interp_new(cbind(depth_dating, age_ensemble[j,]), depth_sample, 'linear',hiatus)
      #print(dim(age))
      sample_ensemble <- cbind(sample_ensemble,age_mc[,2])
    }
  }
  return(sample_ensemble)
  
}

mc_copRa <- function(N, hiatus, depth_dating, age_ensemble, depth_sample){
  for(j in 1:N){
    if(j == 1){
      age_mc <- get_copRa(cbind(depth_dating, age_ensemble[j,]), depth_sample, 'natural',hiatus)
      sample_ensemble <- age_mc
      #print(dim(age))
    } else {
      age_mc <- get_copRa(cbind(depth_dating, age_ensemble[j,]), depth_sample, 'natural',hiatus)
      #print(dim(age))
      sample_ensemble <- cbind(sample_ensemble,age_mc[,2])
    }
  }
  return(sample_ensemble)
  
}


mc_ages <- function(age, age_error, N, c14 = F, filter = T) { # N number of MC simulations
  #set.seed(2019-02-04)
  library(clam)
  
  number <- 0
  d <- 0
  age_ensemble_linReg_final <- NA
  age_ensemble_linInt_final <- NA
  
  idx_rev <- reversal_filter(age, filter)
  print(c('Reversals:',idx_rev))
  
  while(d < N){
    #print('hello')
    k<- N-d
    print(k)
    
    number <- number+1
    print(n)
    
    if(c14){
      age_ensemble <- apply(cbind(age,age_error), 1, function(x)
        sample(calibrate(x[1], x[2], postbomb = 4, cc = 3, graph = FALSE, yrsteps = 0.05)$calib[,1],k, replace =F, 
               prob = calibrate(x[1], x[2], postbomb = 4, cc = 3, graph = FALSE, yrsteps = 0.05)$calib[,2])) 
      age_ensemble_linReg <- age_ensemble
      
    } else {
      
      age_ensemble <- apply(cbind(age,age_error), 1, function(x) rnorm(n = k, mean = x[1], sd = x[2])) # calculate N deviates for each age 
      age_ensemble_linReg <- age_ensemble
      
    }
    
    #age_ensemble_linReg <- age_ensemble
    
    #if(dim(age_ensemble)[1]==dim(age_ensemble_linReg)[1]){print('yes')} else {print('no')}
    
    # test for reversals 
    if (length(idx_rev)>0) {
      if (k==1){
        for (n in 1:dim(idx_rev)[1]){
          age_ensemble[sample(idx_rev[n,], length(idx_rev[n, !is.na(idx_rev[n,])])-1, replace = FALSE)] <- NA
        }
        
      } else {
        for (i in seq(from = 1, to = k)){
          for (n in 1:dim(idx_rev)[1]){
            age_ensemble[i,sample(idx_rev[n,!is.na(idx_rev[n,])], length(idx_rev[n, !is.na(idx_rev[n,])])-1, replace = FALSE)] <- NA
          }
          
        }
      }
    }
    
    #print(age_ensemble)
    
    if(k==1){
      age_ensemble_diff <- diff(age_ensemble)
      run <- T
      
      if (any(age_ensemble_diff < 0 & !is.na(age_ensemble_diff))){
        run <- F
        print('k=1')
      }
      
    } else {
      age_ensemble_diff <- apply(age_ensemble, 1, diff) # each column contains the derivatives for one MC run -> N columns
      run <- T
      if (k ==2) {print('k=2')}
      for (i in seq(1,k)) {
        if(any(age_ensemble_diff[,i] < 0 & !is.na(age_ensemble_diff[,i]))){
          if (k ==2) {print('yes')}
          if (is.null(dim(age_ensemble)[1])) {
            run <- F
            print('dim = 0')
          } else {
            if (k ==2) {print('geloscht')}
            age_ensemble <- age_ensemble[-i,]
            age_ensemble_linReg <- age_ensemble_linReg[-i,]
          }
          
        }
        
      }
      
    }
    
    if (run) {
      if (k == N){
        age_ensemble_linInt_final <- age_ensemble
        age_ensemble_linReg_final <- age_ensemble_linReg
        #if(dim(age_ensemble)[1]==dim(age_ensemble_linReg)[1]){print('yes')}else{print('no')}
      } else {
        age_ensemble_linInt_final <- rbind(age_ensemble_linInt_final, age_ensemble)
        age_ensemble_linReg_final <- rbind(age_ensemble_linReg_final, age_ensemble_linReg)
        #if(dim(age_ensemble)[1]==dim(age_ensemble_linReg)[1]){print('yes')}else{print('no')}
      }
    }
    
    run <- T
    d <- dim(age_ensemble_linInt_final)[1]
    
    
    if(number>10 && d > 1500){return(list(age_ensemble_linReg_final, age_ensemble_linInt_final))
      break
      } else if(number > 10 && d < 1500) {stop('ERROR: too many iterations, check data!')}
    print(c('Dim lI:', dim(age_ensemble_linInt_final)[1]))
    print(c('Dim lR:', dim(age_ensemble_linReg_final)[1]))
  }
  return(list(age_ensemble_linReg_final, age_ensemble_linInt_final)) # ensemble
}
