library(plyr)
library(tidyverse)
wd <- 'wo auch immer deine SISAL csv files liegen'
prefix <- '' # ab v1c heißen die csv files anderes


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
  sample_tb <- join_all(list(sample,hiatus, gap, original_chronology, sisal_chronology, d13C, d18O), by = 'sample_id', type = 'left', match = 'all') %>%
    mutate_at(vars(entity_id, sample_id, sample_thickness, depth_sample, interp_age, interp_age_uncert_pos, interp_age_uncert_neg, COPRA_age,
                   COPRA_age_uncert_pos, COPRA_age_uncert_neg, linear_age, linear_age_uncert_pos, linear_age_uncert_neg, d13C_measurement,
                   d13C_precision, d18O_measurement, d18O_precision), as.numeric)

  return(list(site_tb, dating_tb, sample_tb))
}

data <- load_data(prefix, wd)

site_tb <- as.data.frame(data[1])
dating_tb <- as.data.frame(data[2])
sample_tb <- as.data.frame(data[3])
