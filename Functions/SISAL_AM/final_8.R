library(plyr)
library(tidyverse)

setwd("/home/ariana/R/SISAL_AM")
source(file = 'run_functions.R')
source(file = 'SISAL_AM_add_functions.R')
source(file = 'SISAL_AM_plot_functions.R')
source(file = "StalAge_1_0_mine.R")
source(file = 'New_methods.R')
source(file = 'lin_interp_reversal.R')

#setwd("/home/ariana/Documents/v2")
run<- tibble(entity_id = c(231,232,233,258,363,639,691))
runFile <- tibble(entity_id = run$entity_id, Bacon = rep(T, length(run$entity_id)), Bchron = rep(T, length(run$entity_id)), 
                  StalAge = rep(T, length(run$entity_id)), linInterp = rep(T, length(run$entity_id)), copRa = rep(T, length(run$entity_id)),
                  linReg = rep(T, length(run$entity_id)), 
                  working_directory = rep("/stacycode/ariana/SISAL_final_8", length(run$entity_id)), 
                  wd = rep('/home/ariana/SISAL Data/interimv2', length(run$entity_id)),
                  du = c('date_used_COPRA','date_used_COPRA','date_used_COPRA','date_used_linear','date_used_COPRA','date_used_COPRA','date_used_Bchron'))


load_SISAL <- function(prefix, runFile) {
  
  j <- 1
  
  sapply(1: dim(runFile)[1],
         function(j, x){
           y <- x[j,]
           data <- load_data(prefix,y$wd)
           run_SISAL_chrono(y$entity_id, data, y$working_directory, y$Bacon, y$Bchron, y$StalAge, y$linInterp, y$copRa, y$linReg,j,y$du)
           graphics.off()
           gc()
         },
         x = runFile
  )
  
  
  #data <- load_data(wd)
  #run_SISAL_chrono(entid, data, working_directory, runFile[,2:5])
}

run_SISAL_chrono <- function(entid, data,working_directory, bacon, bchron, stalage, linInterp, copRa, linReg,j,du) {
  file_name <- write_files(entid, data[1], data[2], data[3], bacon, bchron, stalage, linInterp, copRa, linReg, working_directory,du)
  
  err <- NULL
  tryCatch({
    i = 7
    runLinReg(working_directory, file_name)
    print('LinReg done')
    graphics.off()},
    error = function(e){err <<- append(err, paste("ERROR in linReg:",conditionMessage(e), "\n"))
    #error = function(e){print(paste("ERROR in linReg:",conditionMessage(e), "\n"))
    runFile[j,i] <<- F
    linReg <<- F})
  
  tryCatch({
    i = 5
    runLinInterp(working_directory, file_name)
    print('LinInterp done')
    graphics.off()},
    error = function(e){err <<- append(err, paste("ERROR in linInterp:",conditionMessage(e), "\n"))
    #error = function(e){print(paste("ERROR in linInterp:",conditionMessage(e), "\n"))
    runFile[j,i] <<- F
    linInterp <<- F})
  
  tryCatch({
    i = 6
    runcopRa(working_directory, file_name, q1 = 0.05, q2 = 0.95)
    print('copRa done')
    graphics.off()},
    error = function(e){err <<- append(err, paste("ERROR in copRa:",conditionMessage(e), "\n"))
    #error = function(e){print(paste("ERROR in copRa:",conditionMessage(e), "\n"))
    runFile[j,i] <<- F
    copRa <<- F})
  
  tryCatch({
    i = 2
    runBACON(working_directory, file_name)
    print('Bacon done')
    graphics.off()},
    error = function(e){err <<- append(err, paste("ERROR in Bacon:",conditionMessage(e), "\n"))
    #error = function(e){print(paste("ERROR in Bacon:",conditionMessage(e), "\n"))
    runFile[j,i] <<- F
    bacon <<- F})
  
  tryCatch({
    i = 3
    runBchron_new(working_directory, file_name)
    print('Bchron done')
    graphics.off()},
    error = function(e){err <<- append(err, paste("ERROR in Bchron:",conditionMessage(e), "\n"))
    #error = function(e){print(paste("ERROR in Bchron new:",conditionMessage(e), "\n"))
    runFile[j,i] <<- F
    bchron <<- F})
  
  tryCatch({
    i = 4
    runStalAge_new(working_directory, file_name)
    print('StalAge done')
    graphics.off()},
    error = function(e){err <<- append(err, paste("ERROR in StalAge:",conditionMessage(e), "\n"))
    #error = function(e){print(paste("ERROR in StalAge:",conditionMessage(e), "\n"))
    runFile[j,i] <<- F
    stalage <<- F})
  
  setwd(file.path(working_directory, file_name))
  print('wd changed')
  
  tryCatch({
    plotAMC(working_directory, file_name, copra = copRa, b = bacon, bc = bchron, stalage = stalage, linreg = linReg, lininterp = linInterp)
    print('plot done')
    graphics.off()},
    error = function(e){err <<- append(err, paste("ERROR in plot:",conditionMessage(e), "\n"))})
  #error = function(e){print(paste("ERROR in plot:",conditionMessage(e), "\n"))})
  
  write.table(err, 'errors.txt', row.names = F, col.names = F)
  
  setwd(file.path(working_directory))
  write.csv(runFile, 'runFile5.csv', row.names = F, col.names = T)
  
  print(paste(file_name, 'done'))
}

load_SISAL('', runFile) # v1b
setwd("~/Documents/v2")
write.csv(runFile, 'runFile_final.csv', row.names = F, col.names = T)

working_directory <- '/stacycode/ariana/SISAL_final_8/'
file_name <- '691-WS-5d'

tryCatch({
  i = 7
  runLinReg(working_directory, file_name)
  print('LinReg done')
  graphics.off()},
  error = function(e){err <<- append(err, paste("ERROR in linReg:",conditionMessage(e), "\n"))
  #error = function(e){print(paste("ERROR in linReg:",conditionMessage(e), "\n"))
  runFile[j,i] <<- F
  linReg <<- F})

tryCatch({
  i = 5
  runLinInterp(working_directory, file_name)
  print('LinInterp done')
  graphics.off()},
  error = function(e){err <<- append(err, paste("ERROR in linInterp:",conditionMessage(e), "\n"))
  #error = function(e){print(paste("ERROR in linInterp:",conditionMessage(e), "\n"))
  runFile[j,i] <<- F
  linInterp <<- F})

tryCatch({
  i = 6
  runcopRa(working_directory, file_name, q1 = 0.05, q2 = 0.95)
  print('copRa done')
  graphics.off()},
  error = function(e){err <<- append(err, paste("ERROR in copRa:",conditionMessage(e), "\n"))
  #error = function(e){print(paste("ERROR in copRa:",conditionMessage(e), "\n"))
  runFile[j,i] <<- F
  copRa <<- F})

tryCatch({
  i = 2
  runBACON(working_directory, file_name)
  print('Bacon done')
  graphics.off()},
  error = function(e){err <<- append(err, paste("ERROR in Bacon:",conditionMessage(e), "\n"))
  #error = function(e){print(paste("ERROR in Bacon:",conditionMessage(e), "\n"))
  runFile[j,i] <<- F
  bacon <<- F})

tryCatch({
  i = 3
  runBchron_new(working_directory, file_name)
  print('Bchron done')
  graphics.off()},
  error = function(e){err <<- append(err, paste("ERROR in Bchron:",conditionMessage(e), "\n"))
  #error = function(e){print(paste("ERROR in Bchron new:",conditionMessage(e), "\n"))
  runFile[j,i] <<- F
  bchron <<- F})

tryCatch({
  i = 4
  runStalAge_new(working_directory, file_name)
  print('StalAge done')
  graphics.off()},
  error = function(e){err <<- append(err, paste("ERROR in StalAge:",conditionMessage(e), "\n"))
  #error = function(e){print(paste("ERROR in StalAge:",conditionMessage(e), "\n"))
  runFile[j,i] <<- F
  stalage <<- F})

sisalv2interim <- read.SISAL.files('~/SISAL Data/interimv2/',prefix = '')

r <- merge_SISAL_chrono_3(runFile, se,sisalv2interim$entity, sisalv2interim$sample)
final8 <- r[[2]] %>% filter(entity_id %in% run$entity_id)
rnew <- r[[1]]


final8new <- final8 %>% mutate(lin_reg_age = if_else(entity_id %in% c(231,233,363,691,639,258),NA_real_, lin_reg_age),
                               lin_reg_age_uncert_pos = if_else(entity_id %in% c(231,233,363,691,639,258),NA_real_,lin_reg_age_uncert_pos),
                               lin_reg_age_uncert_neg = if_else(entity_id %in% c(231,233,363,691,639,258),NA_real_,lin_reg_age_uncert_neg),
                               lin_interp_age = if_else(entity_id %in% c(233,238,363,691,258),NA_real_,lin_interp_age),
                               lin_interp_age_uncert_pos = if_else(entity_id %in% c(233,238,363,691,258),NA_real_,lin_interp_age_uncert_pos),
                               lin_interp_age_uncert_neg = if_else(entity_id %in% c(233,238,363,691,258),NA_real_,lin_interp_age_uncert_neg),
                               copRa_age = if_else(entity_id %in% c(233,238,363,691,258),NA_real_,copRa_age),
                               copRa_age_uncert_pos = if_else(entity_id %in% c(233,238,363,691,258),NA_real_,copRa_age_uncert_pos),
                               copRa_age_uncert_neg = if_else(entity_id %in% c(233,238,363,691,258),NA_real_,copRa_age_uncert_neg),
                               StalAge_age = if_else(entity_id %in% c(231, 232,233,238,363,691,639,258),NA_real_,StalAge_age),
                               StalAge_age_uncert_pos = if_else(entity_id %in% c(231,232,233,238,363,691,639,258),NA_real_,StalAge_age_uncert_pos),
                               StalAge_age_uncert_neg = if_else(entity_id %in% c(231,232,233,238,363,691,639,258),NA_real_,StalAge_age_uncert_neg),
                               bacon_age = if_else(entity_id %in% c(233,238,363,691,258),NA_real_,bacon_age),
                               bacon_age_uncert_pos = if_else(entity_id %in% c(233,238,363,691,258),NA_real_,bacon_age_uncert_pos),
                               bacon_age_uncert_neg = if_else(entity_id %in% c(233,238,363,691,258),NA_real_,bacon_age_uncert_neg),
                               bchron_age = if_else(entity_id %in% c(233,238,363,691,258),NA_real_,bchron_age),
                               bchron_age_uncert_pos = if_else(entity_id %in% c(233,238,363,691,258),NA_real_,bchron_age_uncert_pos),
                               bchron_age_uncert_neg = if_else(entity_id %in% c(233,238,363,691,258),NA_real_,bchron_age_uncert_neg)) %>%
  mutate(lin_reg_age = if_else(entity_id == 231 & depth_sample < 44.75 & depth_sample > 10,NA_real_, lin_reg_age),
         lin_reg_age_uncert_pos = if_else(entity_id == 231 & depth_sample < 44.75 & depth_sample > 10,NA_real_,lin_reg_age_uncert_pos),
         lin_reg_age_uncert_neg = if_else(entity_id == 231 & depth_sample < 44.75 & depth_sample > 10,NA_real_,lin_reg_age_uncert_neg),
         lin_interp_age = if_else(entity_id == 231 & depth_sample < 44.75 & depth_sample > 10,NA_real_,lin_interp_age),
         lin_interp_age_uncert_pos = if_else(entity_id == 231 & depth_sample < 44.75 & depth_sample > 10,NA_real_,lin_interp_age_uncert_pos),
         lin_interp_age_uncert_neg = if_else(entity_id == 231 & depth_sample < 44.75 & depth_sample > 10,NA_real_,lin_interp_age_uncert_neg),
         copRa_age = if_else(entity_id == 231 & depth_sample < 44.75 & depth_sample > 10,NA_real_,copRa_age),
         copRa_age_uncert_pos = if_else(entity_id == 231 & depth_sample < 44.75 & depth_sample > 10,NA_real_,copRa_age_uncert_pos),
         copRa_age_uncert_neg = if_else(entity_id == 231 & depth_sample < 44.75 & depth_sample > 10,NA_real_,copRa_age_uncert_neg),
         StalAge_age = if_else(entity_id == 231 & depth_sample < 44.75 & depth_sample > 10,NA_real_,StalAge_age),
         StalAge_age_uncert_pos = if_else(entity_id == 231 & depth_sample < 44.75 & depth_sample > 10,NA_real_,StalAge_age_uncert_pos),
         StalAge_age_uncert_neg = if_else(entity_id == 231 & depth_sample < 44.75 & depth_sample > 10,NA_real_,StalAge_age_uncert_neg),
         bacon_age = if_else(entity_id == 231 & depth_sample < 44.75 & depth_sample > 10,NA_real_,bacon_age),
         bacon_age_uncert_pos = if_else(entity_id == 231 & depth_sample < 44.75 & depth_sample > 10,NA_real_,bacon_age_uncert_pos),
         bacon_age_uncert_neg = if_else(entity_id == 231 & depth_sample < 44.75 & depth_sample > 10,NA_real_,bacon_age_uncert_neg),
         bchron_age = if_else(entity_id == 231 & depth_sample < 44.75 & depth_sample > 10,NA_real_,bchron_age),
         bchron_age_uncert_pos = if_else(entity_id == 231 & depth_sample < 44.75 & depth_sample > 10,NA_real_,bchron_age_uncert_pos),
         bchron_age_uncert_neg = if_else(entity_id == 231 & depth_sample < 44.75 & depth_sample > 10,NA_real_,bchron_age_uncert_neg)) %>% 
  mutate(lin_reg_age = if_else(entity_id == 232 & depth_sample < 103.5 & depth_sample > 31,NA_real_, lin_reg_age),
         lin_reg_age_uncert_pos = if_else(entity_id == 232 & depth_sample < 103.5 & depth_sample > 31,NA_real_,lin_reg_age_uncert_pos),
         lin_reg_age_uncert_neg = if_else(entity_id == 232 & depth_sample < 103.5 & depth_sample > 31,NA_real_,lin_reg_age_uncert_neg),
         lin_interp_age = if_else(entity_id == 232 & depth_sample < 103.5 & depth_sample > 31,NA_real_,lin_interp_age),
         lin_interp_age_uncert_pos = if_else(entity_id == 232 & depth_sample < 103.5 & depth_sample > 31,NA_real_,lin_interp_age_uncert_pos),
         lin_interp_age_uncert_neg = if_else(entity_id == 232 & depth_sample < 103.5 & depth_sample > 31,NA_real_,lin_interp_age_uncert_neg),
         copRa_age = if_else(entity_id == 232 & depth_sample < 103.5 & depth_sample > 31,NA_real_,copRa_age),
         copRa_age_uncert_pos = if_else(entity_id == 232 & depth_sample < 103.5 & depth_sample > 31,NA_real_,copRa_age_uncert_pos),
         copRa_age_uncert_neg = if_else(entity_id == 232 & depth_sample < 103.5 & depth_sample > 31,NA_real_,copRa_age_uncert_neg),
         StalAge_age = if_else(entity_id == 232 & depth_sample < 103.5 & depth_sample > 31,NA_real_,StalAge_age),
         StalAge_age_uncert_pos = if_else(entity_id == 232 & depth_sample < 103.5 & depth_sample > 31,NA_real_,StalAge_age_uncert_pos),
         StalAge_age_uncert_neg = if_else(entity_id == 232 & depth_sample < 103.5 & depth_sample > 31,NA_real_,StalAge_age_uncert_neg),
         bacon_age = if_else(entity_id == 232 & depth_sample < 103.5 & depth_sample > 31,NA_real_,bacon_age),
         bacon_age_uncert_pos = if_else(entity_id == 232 & depth_sample < 103.5 & depth_sample > 31,NA_real_,bacon_age_uncert_pos),
         bacon_age_uncert_neg = if_else(entity_id == 232 & depth_sample < 103.5 & depth_sample > 31,NA_real_,bacon_age_uncert_neg),
         bchron_age = if_else(entity_id == 232 & depth_sample < 103.5 & depth_sample > 31,NA_real_,bchron_age),
         bchron_age_uncert_pos = if_else(entity_id == 232 & depth_sample < 103.5 & depth_sample > 31,NA_real_,bchron_age_uncert_pos),
         bchron_age_uncert_neg = if_else(entity_id == 232 & depth_sample < 103.5 & depth_sample > 31,NA_real_,bchron_age_uncert_neg)) 

final8newnew <- left_join(final8new,n, by='sample_id') %>% 
  mutate(lin_reg_age = if_else(!is.na(linear_regress_age) & is.na(lin_reg_age),linear_regress_age, lin_reg_age),
         lin_reg_age_uncert_pos = if_else(!is.na(linear_regress_age_uncert_pos) & is.na(lin_reg_age_uncert_pos),linear_regress_age_uncert_pos, lin_reg_age_uncert_pos),
         lin_reg_age_uncert_neg = if_else(!is.na(linear_regress_age_uncert_neg) & is.na(lin_reg_age_uncert_neg),linear_regress_age_uncert_neg, lin_reg_age_uncert_neg),
         lin_interp_age = if_else(!is.na(linear_age) & is.na(lin_interp_age),linear_age, lin_interp_age),
         lin_interp_age_uncert_pos = if_else(!is.na(linear_age_uncert_pos) & is.na(lin_interp_age_uncert_pos),linear_age_uncert_pos, lin_interp_age_uncert_pos),
         lin_interp_age_uncert_neg = if_else(!is.na(linear_age_uncert_neg) & is.na(lin_interp_age_uncert_neg),linear_age_uncert_neg, lin_interp_age_uncert_neg),
         copRa_age = if_else(!is.na(COPRA_age) & is.na(copRa_age),COPRA_age, copRa_age),
         copRa_age_uncert_pos = if_else(!is.na(COPRA_age_uncert_pos) & is.na(copRa_age_uncert_pos),COPRA_age_uncert_pos, copRa_age_uncert_pos),
         copRa_age_uncert_neg = if_else(!is.na(COPRA_age_uncert_neg) & is.na(copRa_age_uncert_neg),COPRA_age_uncert_neg, copRa_age_uncert_neg),
         bacon_age = if_else(!is.na(Bacon_age) & is.na(bacon_age),Bacon_age, bacon_age),
         bacon_age_uncert_pos = if_else(!is.na(Bacon_age_uncert_pos) & is.na(bacon_age_uncert_pos),Bacon_age_uncert_pos, bacon_age_uncert_pos),
         bacon_age_uncert_neg = if_else(!is.na(Bacon_age_uncert_neg) & is.na(bacon_age_uncert_neg),Bacon_age_uncert_neg, bacon_age_uncert_neg),
         bchron_age = if_else(!is.na(Bchron_age) & is.na(bchron_age),Bchron_age, bchron_age),
         bchron_age_uncert_pos = if_else(!is.na(Bchron_age_uncert_pos) & is.na(bchron_age_uncert_pos),Bchron_age_uncert_pos, bchron_age_uncert_pos),
         bchron_age_uncert_neg = if_else(!is.na(Bchron_age_uncert_neg) & is.na(bchron_age_uncert_neg),Bchron_age_uncert_neg, bchron_age_uncert_neg)) %>% 
  select(entity_id, sample_id, depth_sample, lin_reg_age, lin_reg_age_uncert_pos, lin_reg_age_uncert_neg, lin_interp_age, lin_interp_age_uncert_pos, lin_interp_age_uncert_neg,
         copRa_age, copRa_age_uncert_pos, copRa_age_uncert_neg, StalAge_age, StalAge_age_uncert_pos, StalAge_age_uncert_neg,
         bacon_age, bacon_age_uncert_pos, bacon_age_uncert_neg, bchron_age, bchron_age_uncert_pos, bchron_age_uncert_neg) %>% 
  mutate(OxCal_age = rep(NA_real_, length(final8new$sample_id)),
         OxCal_age_uncert_pos = rep(NA_real_, length(final8new$sample_id)),
         OxCal_age_uncert_neg = rep(NA_real_, length(final8new$sample_id))) %>%
  filter(!(sample_id %in% sisalv2interim$hiatus$sample_id))
  

s632 <- sisalv2interim$sample %>% filter(entity_id == 632) %>% filter(!(sample_id %in% sisalv2interim$hiatus$sample_id))
f632 <- sisalv2interim$sisal_chronology %>% filter(sample_id %in% s632$sample_id) %>%
  dplyr::rename(lin_reg_age = linear_regress_age, lin_reg_age_uncert_pos = linear_regress_age_uncert_pos, lin_reg_age_uncert_neg = linear_regress_age_uncert_neg) %>% 
  mutate(entity_id = rep(632,length(s632$sample_id)), depth_sample = s632$depth_sample) %>% select(entity_id, sample_id, depth_sample,lin_reg_age, lin_reg_age_uncert_pos, lin_reg_age_uncert_neg)
final8newnewnew <- bind_rows(final8newnew,f632)

s <- sisalv2interim$dating %>% filter(entity_id %in% c(231,232,233,258,363,632,639,691))

final8date_file <- tibble(entity_id = s$entity_id, dating_id = s$dating_id) %>%
  mutate(date_used_lin_reg = c((s %>% filter(entity_id %in% c(231,232,233)))$date_used_COPRA,(s %>% filter(entity_id == 258))$date_used_linear, (s %>% filter(entity_id == 363))$date_used_COPRA,
                               (s %>% filter(entity_id==632))$date_used_linear_regress,
                               (s %>% filter(entity_id == 639))$date_used_COPRA, (s %>% filter(entity_id == 691))$date_used_Bchron),
         date_used_lin_interp = c((s %>% filter(entity_id %in% c(231,232,233)))$date_used_COPRA,(s %>% filter(entity_id == 258))$date_used_linear, (s %>% filter(entity_id == 363))$date_used_COPRA,
                                  (s %>% filter(entity_id==632))$date_used_linear_regress,
                                  (s %>% filter(entity_id == 639))$date_used_COPRA, (s %>% filter(entity_id == 691))$date_used_Bchron), 
         date_used_copRa = c((s %>% filter(entity_id %in% c(231,232,233)))$date_used_COPRA,(s %>% filter(entity_id == 258))$date_used_linear, (s %>% filter(entity_id == 363))$date_used_COPRA,
                             (s %>% filter(entity_id==632))$date_used_linear_regress,
                             (s %>% filter(entity_id == 639))$date_used_COPRA, (s %>% filter(entity_id == 691))$date_used_Bchron), 
         date_used_StalAge = c((s %>% filter(entity_id %in% c(231,232,233)))$date_used_COPRA,(s %>% filter(entity_id == 258))$date_used_linear, (s %>% filter(entity_id == 363))$date_used_COPRA,
                               (s %>% filter(entity_id==632))$date_used_linear_regress,
                               (s %>% filter(entity_id == 639))$date_used_COPRA, (s %>% filter(entity_id == 691))$date_used_Bchron), 
         date_used_bacon = c((s %>% filter(entity_id %in% c(231,232,233)))$date_used_COPRA,(s %>% filter(entity_id == 258))$date_used_linear, (s %>% filter(entity_id == 363))$date_used_COPRA,
                             (s %>% filter(entity_id==632))$date_used_linear_regress,
                             (s %>% filter(entity_id == 639))$date_used_COPRA, (s %>% filter(entity_id == 691))$date_used_Bchron), 
         date_used_bchron = c((s %>% filter(entity_id %in% c(231,232,233)))$date_used_COPRA,(s %>% filter(entity_id == 258))$date_used_linear, (s %>% filter(entity_id == 363))$date_used_COPRA,
                              (s %>% filter(entity_id==632))$date_used_linear_regress,
                              (s %>% filter(entity_id == 639))$date_used_COPRA, (s %>% filter(entity_id == 691))$date_used_Bchron)) %>%
  mutate(date_used_lin_reg = if_else(entity_id %in% c(231,233,363,691,639,258),'NULL', date_used_lin_reg),
         date_used_lin_interp = if_else(entity_id %in% c(233,632,363,691,231),'NULL',date_used_lin_interp),
         date_used_copRa = if_else(entity_id %in% c(258,632,691),'NULL',date_used_copRa),
         date_used_StalAge = if_else(entity_id %in% c(231,232,233,258,363,632,639,691),'NULL',date_used_StalAge),
         date_used_bacon = if_else(entity_id %in% c(233,258,363,632),'NULL',date_used_bacon),
         date_used_bchron = if_else(entity_id %in% c(233,258,363,632),'NULL',date_used_bchron),
         date_used_OxCal = rep('NULL',94))




