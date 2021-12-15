# packages to install
#rJava, xlsx, plyr, tidyverse, rbacon, bchron, clam, tidyr, tibble, Hmisc, plotrix
#library(rJava)
#library(xlsx)
library(plyr)
library(tidyverse)

setwd("/home/ariana/R/SISAL_AM")
source(file = 'run_functions.R')
source(file = 'SISAL_AM_add_functions.R')
source(file = 'SISAL_AM_plot_functions.R')
source(file = "StalAge_1_0_mine.R")
source(file = 'New_methods.R')
source(file = 'lin_interp_reversal.R')

setwd("~/Documents/rbacon test version 2.3.9.1")
#run<- read.csv('run.csv', header = T)
#run<- read.csv('run2.csv', header = T)
#run <- run %>% filter(entity_id %in% c(2,4,8,12,13,17,26))
#run <- run %>% filter(entity_id %in% c(237,319,320, 305,312))
run <- data.frame(entity_id = c(237,319,320,305,312))
run <- run1 %>% arrange(., entity_id)
runFile <- tibble(entity_id = run$entity_id, Bacon = rep(T, length(run$entity_id)), Bchron = rep(T, length(run$entity_id)), 
                  StalAge = rep(T, length(run$entity_id)), linInterp = rep(T, length(run$entity_id)), copRa = rep(T, length(run$entity_id)),
                  linReg = rep(T, length(run$entity_id)), 
                  working_directory = rep("/home/ariana/Documents/full Run 2", length(run$entity_id)), 
                  wd = rep('/home/ariana/SISAL Data/sisalv1b_csv', length(run$entity_id)))
runFile <- tibble(entity_id = run$entity_id, Bacon = rep(T, length(run$entity_id)), Bchron = rep(F, length(run$entity_id)), 
                  StalAge = rep(F, length(run$entity_id)), linInterp = rep(F, length(run$entity_id)),linReg = rep(F, length(run$entity_id)), 
                  working_directory = rep("/home/ariana/Documents/rbacon test version 2.3.9.1", length(run$entity_id)), 
                  wd = rep('/home/ariana/SISAL Data/sisalv1b_csv', length(run$entity_id)))
runFile <- tibble(entity_id = run$entity_id, Bacon = rep(T, length(run$entity_id)), Bchron = rep(T, length(run$entity_id)), 
                  StalAge = rep(T, length(run$entity_id)), linInterp = rep(T, length(run$entity_id)), copRa = rep(T, length(run$entity_id)), 
                  linReg = rep(T, length(run$entity_id)), 
                  working_directory = rep("/home/ariana/Documents/SISAL v1c", length(run$entity_id)), 
                  wd = rep("~/SISAL Data/v1c_20190712", length(run$entity_id)))
runFile <- tibble(entity_id = run$entity_id, Bacon = rep(T, length(run$entity_id)), Bchron = rep(F, length(run$entity_id)), 
                  StalAge = rep(F, length(run$entity_id)), linInterp = rep(F, length(run$entity_id)), copRa = rep(F, length(run$entity_id)), 
                  linReg = rep(F, length(run$entity_id)), 
                  working_directory = rep("/home/ariana/Documents/test final", length(run$entity_id)), 
                  wd = rep("~/SISAL Data/sisalv1b_csv", length(run$entity_id)))

#working_directory <- "/home/ariana/Documents/test new hiatus"
wd <- '/home/ariana/SISAL Data/sisalv1b_csv'
wd <- "~/SISAL Data/v1c_20190712"
wd <- "~/SISAL Data/interim_v2_csv_20190812"
prefix <- 'v1c_SISAL_'
working_directory <- "/home/ariana/Documents/SISAL v1c"
file_name <- '11-KS06-A-H'


load_SISAL <- function(prefix, runFile) {
  
  j <- 1
  
  sapply(1: dim(runFile)[1],
        function(j, x){
          y <- x[j,]
          data <- load_data(prefix,y$wd)
          run_SISAL_chrono(y$entity_id, data, y$working_directory, y$Bacon, y$Bchron, y$StalAge, y$linInterp, y$copRa, y$linReg,j)
          graphics.off()
          gc()
        },
        x = runFile
  )
  
  
  #data <- load_data(wd)
  #run_SISAL_chrono(entid, data, working_directory, runFile[,2:5])
}

run_SISAL_chrono <- function(entid, data,working_directory, bacon, bchron, stalage, linInterp, copRa, linReg,j) {
  file_name <- write_files(entid, data[1], data[2], data[3], bacon, bchron, stalage, linInterp, copRa, linReg, working_directory)
  
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
  write.csv(runFile, 'runFile.csv', row.names = F, col.names = T)
  
  print(paste(file_name, 'done'))
}

load_SISAL('', runFile) # v1b
load_SISAL('v1c_SISAL_', runFile) # v1c
setwd(file.path(working_directory))
write.csv(runFile, 'runFile9.csv', row.names = F, col.names = T)



# 10: MC not working
# 11, 14, 16, 19, 20: Bacon not working
# 48: no depth_sample -> nothing working
# 72: aufgehängt
# 48, 207, 261, 281, 322, 327, 351, 417: keine depth_sample
setwd('/home/ariana/SISAL Data/sisalv1b_csv')
sample <- read.csv('sample.csv', header = T, stringsAsFactors = F)#, colClasses = c('numeric', 'numeric', 'character', 'numeric','character', 'character'))
entity <- read.csv('entity.csv', header = T, stringsAsFactors = F)
dating <- read.csv('dating.csv', header = T, stringsAsFactors = F) 
dating <- dating %>% mutate_at(vars(dating_id, entity_id, depth_dating, dating_thickness, X14C_correction, corr_age, corr_age_uncert_pos, corr_age_uncert_neg), as.numeric)
se <- left_join(sample, entity, by = 'entity_id')
SISAL_chronology_new_bacon <- data.frame(entity_id = se$entity_id, sample_id = se$sample_id, 
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
                               StalAge_age_uncert_neg = rep(NA_real_, length(sample$sample_id)))
setwd("~/Documents/lin_interp_new_scan")
r1 <- read.csv('runFile1.csv', header = T) %>% mutate_at(vars(working_directory), as.character) %>% select(-wd)
setwd("~/Documents/Hiatus & Reversals (non tractable)")
r2 <- read.csv('runFile1.csv', header = T) %>% mutate_at(vars(working_directory), as.character) %>% select( -wd)
s1 <- merge_SISAL_chrono(r1,SISAL_chronology_new, entity)
s2 <- merge_SISAL_chrono(r2,SISAL_chronology_new, entity)
test <- merge_SISAL_chrono(runFile %>% select(-wd), SISAL_chronology_new_bacon, entity)
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

get_growth_rates <- function(SISAL_chronology){
  SISAL_chronology_GR <- SISAL_chronology %>% mutate(bacon_gr = ()/(lead(bacon_age) -bacon_age)
  
  
} 


#depth <- data.frame(depth_sample = depth)
s1_new <- s1 %>% select(sample_id, lin_interp_age, lin_interp_age_uncert_pos, lin_interp_age_uncert_neg )
s2_new <- s2 %>% select(-c(lin_interp_age, lin_interp_age_uncert_pos, lin_interp_age_uncert_neg))
s <- left_join(s1_new, s2_new, by="sample_id") %>% filter(entity_id %in% r1$entity_id)

s_new <- s2
depth <- sample %>% mutate_at(vars(sample_id, entity_id, depth_sample, sample_thickness), as.numeric) %>% filter(sample_id %in% s_new$sample_id) %>% select(sample_id, depth_sample) 
sc_new_new_2 <- left_join(s2, depth, by = "sample_id") %>% select(entity_id, sample_id, depth_sample, everything()) %>% filter(entity_id %in% r2$entity_id)
sc_new_new_1 <- left_join(s1, depth, by = "sample_id") %>% select(entity_id, sample_id, depth_sample, everything()) %>% filter(entity_id %in% r1$entity_id)

eval_SISAL_chrono <- function(SISAL_chronology, depths){
  
  s_new <- SISAL_chronology
  s_new_fil_1 <- s1 %>% filter(entity_id %in% r1$entity_id)
  SISAL_chronology_eval_reversal_1 <- s_new_fil_1 %>% ungroup() %>% group_by(entity_id) %>% 
    mutate(., diff_lr = lead(lin_reg_age) - lin_reg_age,
           diff_li = lead(lin_interp_age) - lin_interp_age,
           diff_bacon = lead(bacon_age) - bacon_age,
           diff_bchron = lead(bchron_age) - bchron_age,
           diff_stalage = lead(StalAge_age) - StalAge_age) %>%
    mutate(rev_lr = if_else( diff_lr< 0, FALSE, TRUE),
           rev_li = if_else( diff_li< 0, FALSE, TRUE),
           rev_bacon = if_else( diff_bacon< 0, FALSE, TRUE),
           rev_bchron = if_else( diff_bchron< 0, FALSE, TRUE),
           rev_StalAge = if_else( diff_stalage< 0, FALSE, TRUE)) %>%
    group_by(entity_id) %>% 
    summarise(lR_reversal = if_else(all(is.na(rev_lr)),-1, if_else(any(!(rev_lr), na.rm = T), 0, 1)),
              lI_reversal = if_else(all(is.na(rev_li)),-1, if_else(any(!(rev_li), na.rm = T), 0, 1)),
              bacon_reversal = if_else(all(is.na(rev_bacon)),-1, if_else(any(!(rev_bacon), na.rm = T), 0, 1)),
              bchron_reversal = if_else(all(is.na(rev_bchron)),-1, if_else(any(!(rev_bchron), na.rm = T), 0, 1)),
              stalage_reversal = if_else(all(is.na(rev_StalAge)),-1, if_else(any(!(rev_StalAge), na.rm = T), 0, 1))) %>% 
    select(entity_id, lI_reversal)
  
  s_new_fil_2 <- s2 %>% filter(entity_id %in% r2$entity_id)
  SISAL_chronology_eval_reversal_2 <- s_new_fil_2 %>% ungroup() %>% group_by(entity_id) %>% 
    mutate(., diff_lr = lead(lin_reg_age) - lin_reg_age,
              diff_li = lead(lin_interp_age) - lin_interp_age,
              diff_bacon = lead(bacon_age) - bacon_age,
              diff_bchron = lead(bchron_age) - bchron_age,
              diff_stalage = lead(StalAge_age) - StalAge_age) %>%
    mutate(rev_lr = if_else( diff_lr< 0, FALSE, TRUE),
           rev_li = if_else( diff_li< 0, FALSE, TRUE),
           rev_bacon = if_else( diff_bacon< 0, FALSE, TRUE),
           rev_bchron = if_else( diff_bchron< 0, FALSE, TRUE),
           rev_StalAge = if_else( diff_stalage< 0, FALSE, TRUE)) %>%
    group_by(entity_id) %>% 
    summarise(lR_reversal = if_else(all(is.na(rev_lr)),-1, if_else(any(!(rev_lr), na.rm = T), 0, 1)),
              lI_reversal = if_else(all(is.na(rev_li)),-1, if_else(any(!(rev_li), na.rm = T), 0, 1)),
              bacon_reversal = if_else(all(is.na(rev_bacon)),-1, if_else(any(!(rev_bacon), na.rm = T), 0, 1)),
              bchron_reversal = if_else(all(is.na(rev_bchron)),-1, if_else(any(!(rev_bchron), na.rm = T), 0, 1)),
              stalage_reversal = if_else(all(is.na(rev_StalAge)),-1, if_else(any(!(rev_StalAge), na.rm = T), 0, 1))) %>% select(-lI_reversal)
    
  SISAL_eval_rev <- left_join(SISAL_chronology_eval_reversal_1, SISAL_chronology_eval_reversal_2, by = 'entity_id')
    
    data2 <- SISAL_eval_rev %>% filter(entity_id %in% r1$entity_id)
    mine.data2 <- as.data.frame(gather(data = data2, key = Class, value = Fit, c(-1)))
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
  
  dt <- dating %>% filter(entity_id %in% s_new$entity_id & date_type != 'Event; hiatus')
  
  d <- data.frame(entity_id = dt$entity_id, 
                  sample_id = dt$dating_id, 
                  depth_sample = dt$depth_dating,
                  lin_reg_age = dt$corr_age,
                  lin_reg_age_uncert_pos = dt$corr_age_uncert_pos,
                  lin_reg_age_uncert_neg = dt$corr_age_uncert_neg,
                  lin_interp_age = dt$corr_age,
                  lin_interp_age_uncert_pos = dt$corr_age_uncert_pos,
                  lin_interp_age_uncert_neg = dt$corr_age_uncert_neg,
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
  d_fil <- d %>% filter(entity_id %in% r1$entity_id)
  
  SISAL_chronology_eval_fit_1 <- bind_rows(sc_new_new_1, d_fil) %>% ungroup() %>% group_by(entity_id) %>% dplyr::arrange(., depth_sample, .by_group =T) %>%
    mutate(rev_lr = if_else(abs(lin_reg_age-lead(lin_reg_age)) < (2*lin_reg_age_uncert_pos + 2*lead(lin_reg_age_uncert_neg)), TRUE, FALSE),
           rev_li = if_else(abs(lin_interp_age-lead(lin_interp_age)) < (2*lin_interp_age_uncert_pos + 2*lead(lin_interp_age_uncert_neg)), TRUE, FALSE),
           rev_bacon = if_else(abs(bacon_age-lead(bacon_age)) < (2*bacon_age_uncert_pos + 2*lead(bacon_age_uncert_neg)), TRUE, FALSE),
           rev_bchron = if_else(abs(bchron_age-lead(bchron_age)) < (2*bchron_age_uncert_pos + 2*lead(bchron_age_uncert_neg)), TRUE, FALSE),
           rev_StalAge = if_else(abs(StalAge_age-lead(StalAge_age)) < (2*StalAge_age_uncert_pos + 2*lead(StalAge_age_uncert_neg)), TRUE, FALSE)) %>%
    summarise(lR_count = sum(!rev_lr, na.rm = T),
              lI_count = sum(!rev_li, na.rm = T),
              bacon_count = sum(!rev_bacon, na.rm = T),
              bchron_count = sum(!rev_bchron, na.rm = T),
              stalage_count = sum(!rev_StalAge, na.rm = T)) %>%
    left_join(., dt_count, by = "entity_id") %>% group_by(entity_id) %>% summarise(lR_fitness = if_else(lR_count <= 0.3*n, 1, 0),
                                                                                   lI_fitness = if_else(lI_count <= 0.3*n, 1, 0),
                                                                                   bacon_fitness = if_else(bacon_count <= 0.3*n, 1, 0),
                                                                                   bchron_fitness = if_else(bchron_count <= 0.3*n, 1, 0),
                                                                                   stalage_fitness = if_else(stalage_count <= 0.3*n, 1, 0)) %>%
    select(entity_id, lI_fitness)
  
  SISAL_chronology_eval_fit_2 <- bind_rows(sc_new_new_2, d_fil) %>% ungroup() %>% group_by(entity_id) %>% dplyr::arrange(., depth_sample, .by_group =T) %>%
    mutate(rev_lr = if_else(abs(lin_reg_age-lead(lin_reg_age)) < (2*lin_reg_age_uncert_pos + 2*lead(lin_reg_age_uncert_neg)), TRUE, FALSE),
           rev_li = if_else(abs(lin_interp_age-lead(lin_interp_age)) < (2*lin_interp_age_uncert_pos + 2*lead(lin_interp_age_uncert_neg)), TRUE, FALSE),
           rev_bacon = if_else(abs(bacon_age-lead(bacon_age)) < (2*bacon_age_uncert_pos + 2*lead(bacon_age_uncert_neg)), TRUE, FALSE),
           rev_bchron = if_else(abs(bchron_age-lead(bchron_age)) < (2*bchron_age_uncert_pos + 2*lead(bchron_age_uncert_neg)), TRUE, FALSE),
           rev_StalAge = if_else(abs(StalAge_age-lead(StalAge_age)) < (2*StalAge_age_uncert_pos + 2*lead(StalAge_age_uncert_neg)), TRUE, FALSE)) %>%
    summarise(lR_count = sum(!rev_lr, na.rm = T),
              lI_count = sum(!rev_li, na.rm = T),
              bacon_count = sum(!rev_bacon, na.rm = T),
              bchron_count = sum(!rev_bchron, na.rm = T),
              stalage_count = sum(!rev_StalAge, na.rm = T)) %>%
    left_join(., dt_count, by = "entity_id") %>% group_by(entity_id) %>% summarise(lR_fitness = if_else(lR_count <= 0.3*n, 1, 0),
                                                                                   lI_fitness = if_else(lI_count <= 0.3*n, 1, 0),
                                                                                   bacon_fitness = if_else(bacon_count <= 0.3*n, 1, 0),
                                                                                   bchron_fitness = if_else(bchron_count <= 0.3*n, 1, 0),
                                                                                   stalage_fitness = if_else(stalage_count <= 0.3*n, 1, 0)) %>% 
    select(-lI_fitness)
  SISAL_eval_fit <- left_join(SISAL_chronology_eval_fit_1,SISAL_chronology_eval_fit_2, by = 'entity_id')
  data <- SISAL_eval_fit %>% filter(entity_id %in% r1$entity_id)
  mine.data <- gather(data = data, key = Class, value = Fit, c(-1))
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
  
  SISAL_dates <-  bind_rows(sc_new_new, d) %>% ungroup() %>% group_by(entity_id) %>% arrange(., depth_sample, .by_group =T) %>% left_join(., d_new, by = "sample_id") %>%
    mutate(depth_new_d = sapply(depth_new, function(x) depth_sample[order(abs(x - depth_sample))][1])) %>%
    mutate(depth_new_d = if_else(!is.na(depth_new), depth_new_d,  NA_real_)) #%>% 
    mutate(date_iqr = if_else(!is.na(depth_new),(lin_reg_age_uncert_pos+lin_reg_age_uncert_neg)/2, NA_real_),
           lR_iqr = if_else(depth_sample %in% depth_new_d,(lin_reg_age_uncert_pos+lin_reg_age_uncert_neg)/2 , NA_real_),
           lI_iqr = if_else(depth_sample %in% depth_new_d, (lin_interp_age_uncert_pos+lin_interp_age_uncert_neg)/2, NA_real_),
           bacon_iqr = if_else(depth_sample %in% depth_new_d, (bacon_age_uncert_pos+bacon_age_uncert_neg)/2, NA_real_),
           bchron_iqr = if_else(depth_sample %in% depth_new_d, (bchron_age_uncert_pos+bchron_age_uncert_neg)/2, NA_real_),
           StalAge_iqr = if_else(depth_sample %in% depth_new_d, (StalAge_age_uncert_pos+StalAge_age_uncert_neg)/2, NA_real_))
    
    sc_new_new_fil <- sc_new_new %>% filter(entity_id %in% runFile$entity_id)
    d_fil <- d %>% filter(entity_id %in% runFile$entity_id)
    SISAL_IQR_check_1 <- bind_rows(sc_new_new_1, d_fil) %>% ungroup() %>% group_by(entity_id) %>% dplyr::arrange(., depth_sample, .by_group =T) %>% left_join(., d_new, by = "sample_id") %>% 
      mutate(depth_new_d = sapply(depth_new, function(x) depth_sample[order(abs(x - depth_sample))][1])) %>% 
      mutate(depth_new_d = if_else(!is.na(depth_new), depth_new_d,  NA_real_)) %>% nest() %>% 
      mutate(filter = purrr::map(.$data, function(x) filter(x, !is.na(depth_new_d))$depth_new_d)) %>% 
      mutate(proxies = purrr::map2(.$data, .$filter, function(x, y) filter(x, depth_sample %in% y & is.na(depth_new) & is.na(depth_new_d)) %>% mutate( 
        lR_iqr = (lin_reg_age_uncert_pos+lin_reg_age_uncert_neg)/2 , 
        lI_iqr =  (lin_interp_age_uncert_pos+lin_interp_age_uncert_neg)/2, 
        bacon_iqr =  (bacon_age_uncert_pos+bacon_age_uncert_neg)/2, 
        bchron_iqr =  (bchron_age_uncert_pos+bchron_age_uncert_neg)/2, 
        StalAge_iqr =  (StalAge_age_uncert_pos+StalAge_age_uncert_neg)/2) %>% select(sample_id, depth_sample, lR_iqr, lI_iqr, bacon_iqr, bchron_iqr, StalAge_iqr))) %>% 
      mutate(filter2 = purrr::map(.$data, function(x) filter(x, !is.na(depth_new))$sample_id)) %>% 
      mutate(dating =  purrr::map2(.$data, .$filter2, function(x, y) filter(x, sample_id %in% y) %>% mutate( 
        date_iqr = (lin_reg_age_uncert_pos+lin_reg_age_uncert_neg)/2) %>% select(sample_id, depth_new_d, date_iqr))) %>%
      mutate(final = purrr::map2(.$proxies, .$dating, function(x,y) left_join(y %>% rename(depth_join = depth_new_d),x %>% rename(depth_join = depth_sample), by= 'depth_join'))) %>%
      unnest(final) %>%
      mutate(lR = if_else((date_iqr + lead(date_iqr))/2<lR_iqr, 1,0),
             lI = if_else((date_iqr + lead(date_iqr))/2<lI_iqr , 1,0),
             bacon = if_else((date_iqr + lead(date_iqr))/2< bacon_iqr , 1,0),
             bchron = if_else((date_iqr + lead(date_iqr))/2<bchron_iqr , 1,0),
             stalage = if_else((date_iqr + lead(date_iqr))/2<StalAge_iqr, 1,0)) %>% group_by(entity_id) %>% 
      summarise(general = n(),
                lR_iqr_1 = if_else(sum(lR, na.rm = T)>= 0.6*general, 1,0),
                lI_iqr_1 = if_else(sum(lI, na.rm = T)>= 0.6*general, 1, 0),
                bacon_iqr_1 = if_else(sum(bacon, na.rm = T)>= 0.6*general,1,0),
                bchron_iqr_1 = if_else(sum(bchron, na.rm = T)>= 0.6*general,1,0),
                stalage_iqr_1 = if_else(sum(stalage, na.rm = T)>= 0.6*general,1,0)) %>%
      select(entity_id, lI_iqr_1)
    
    SISAL_IQR_check_2 <- bind_rows(sc_new_new_2, d_fil) %>% ungroup() %>% group_by(entity_id) %>% dplyr::arrange(., depth_sample, .by_group =T) %>% left_join(., d_new, by = "sample_id") %>% 
      mutate(depth_new_d = sapply(depth_new, function(x) depth_sample[order(abs(x - depth_sample))][1])) %>% 
      mutate(depth_new_d = if_else(!is.na(depth_new), depth_new_d,  NA_real_)) %>% nest() %>% 
      mutate(filter = purrr::map(.$data, function(x) filter(x, !is.na(depth_new_d))$depth_new_d)) %>% 
      mutate(proxies = purrr::map2(.$data, .$filter, function(x, y) filter(x, depth_sample %in% y & is.na(depth_new) & is.na(depth_new_d)) %>% mutate( 
        lR_iqr = (lin_reg_age_uncert_pos+lin_reg_age_uncert_neg)/2 , 
        lI_iqr =  (lin_interp_age_uncert_pos+lin_interp_age_uncert_neg)/2, 
        bacon_iqr =  (bacon_age_uncert_pos+bacon_age_uncert_neg)/2, 
        bchron_iqr =  (bchron_age_uncert_pos+bchron_age_uncert_neg)/2, 
        StalAge_iqr =  (StalAge_age_uncert_pos+StalAge_age_uncert_neg)/2) %>% select(sample_id, depth_sample, lR_iqr, lI_iqr, bacon_iqr, bchron_iqr, StalAge_iqr))) %>% 
      mutate(filter2 = purrr::map(.$data, function(x) filter(x, !is.na(depth_new))$sample_id)) %>% 
      mutate(dating =  purrr::map2(.$data, .$filter2, function(x, y) filter(x, sample_id %in% y) %>% mutate( 
        date_iqr = (lin_reg_age_uncert_pos+lin_reg_age_uncert_neg)/2) %>% select(sample_id, depth_new_d, date_iqr))) %>%
      mutate(final = purrr::map2(.$proxies, .$dating, function(x,y) left_join(y %>% rename(depth_join = depth_new_d),x %>% rename(depth_join = depth_sample), by= 'depth_join'))) %>%
      unnest(final) %>%
      mutate(lR = if_else((date_iqr + lead(date_iqr))/2<lR_iqr, 1,0),
             lI = if_else((date_iqr + lead(date_iqr))/2<lI_iqr , 1,0),
             bacon = if_else((date_iqr + lead(date_iqr))/2< bacon_iqr , 1,0),
             bchron = if_else((date_iqr + lead(date_iqr))/2<bchron_iqr , 1,0),
             stalage = if_else((date_iqr + lead(date_iqr))/2<StalAge_iqr, 1,0)) %>% group_by(entity_id) %>% 
      summarise(general = n(),
                lR_iqr_1 = if_else(sum(lR, na.rm = T)>= 0.6*general, 1,0),
                lI_iqr_1 = if_else(sum(lI, na.rm = T)>= 0.6*general, 1, 0),
                bacon_iqr_1 = if_else(sum(bacon, na.rm = T)>= 0.6*general,1,0),
                bchron_iqr_1 = if_else(sum(bchron, na.rm = T)>= 0.6*general,1,0),
                stalage_iqr_1 = if_else(sum(stalage, na.rm = T)>= 0.6*general,1,0)) %>%
      select(-lI_iqr_1)
    
    SISAL_IQR_check <- left_join(SISAL_IQR_check_1, SISAL_IQR_check_2, by='entity_id')
    
    data3 <- SISAL_IQR_check %>% select(entity_id, lR_iqr_1, lI_iqr_1, bacon_iqr_1, bchron_iqr_1, stalage_iqr_1)
    mine.data3 <- gather(data = data3, key = Class, value = Fit, c(-1))
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
    
    
    SISAL_IQR_1 <- sc_new_new_1 %>% ungroup() %>% group_by(entity_id) %>% arrange(., depth_sample, .by_group =T) %>%
      mutate(lR_iqr = (lin_reg_age_uncert_pos+lin_reg_age_uncert_neg)/2 , 
             lI_iqr =  (lin_interp_age_uncert_pos+lin_interp_age_uncert_neg)/2, 
             bacon_iqr =  (bacon_age_uncert_pos+bacon_age_uncert_neg)/2, 
             bchron_iqr =  (bchron_age_uncert_pos+bchron_age_uncert_neg)/2, 
             StalAge_iqr =  (StalAge_age_uncert_pos+StalAge_age_uncert_neg)/2)
    
    SISAL_young_1 <- SISAL_IQR_1 %>% group_by(entity_id) %>% mutate(general = n()/3) %>%
      slice(.,1:general) %>% group_by(entity_id) %>% summarise(lR_young = mean(lR_iqr, na.rm = T),
                                                               lI_young = mean(lI_iqr, na.rm = T),
                                                               bacon_young = mean(bacon_iqr, na.rm = T),
                                                               bchron_young = mean(bchron_iqr, na.rm = T),
                                                               stalage_young = mean(StalAge_iqr, na.rm = T))
    
    SISAL_old_1 <- SISAL_IQR_1 %>% group_by(entity_id) %>% mutate(general = n()/3) %>%
      slice(.,general*2:n()) %>% group_by(entity_id) %>% summarise(lR_old = mean(lR_iqr, na.rm = T),
                                                               lI_old = mean(lI_iqr, na.rm = T),
                                                               bacon_old = mean(bacon_iqr, na.rm = T),
                                                               bchron_old = mean(bchron_iqr, na.rm = T),
                                                               stalage_old = mean(StalAge_iqr, na.rm = T))
    
    SISAL_IQR_check2_1 <- left_join(SISAL_young_1, SISAL_old_1, by = "entity_id") %>% mutate(lR_iqr_2 = if_else(lR_young < lR_old , 1, 0),
                                                                                       lI_iqr_2 = if_else(lI_young < lI_old, 1,0),
                                                                                       bacon_iqr_2 = if_else(bacon_young < bacon_old,1,0),
                                                                                       bchron_iqr_2 = if_else(bchron_young < bchron_old,1,0),
                                                                                       stalage_iqr_2 = if_else(stalage_young < stalage_old, 1,0)) %>%
      mutate(lR_iqr_2 = if_else(is.na(lR_young) & is.na(lR_old), 0, lR_iqr_2),
             lI_iqr_2 = if_else(is.na(lI_young) & is.na(lI_old), 0, lI_iqr_2),
             bacon_iqr_2 = if_else(is.na(bacon_young) & is.na(bacon_old), 0, bacon_iqr_2),
             bchron_iqr_2 = if_else(is.na(bchron_young) & is.na(bchron_old),0,bchron_iqr_2),
             stalage_iqr_2 = if_else(is.na(stalage_young) & is.na(stalage_old),0,stalage_iqr_2)) %>% 
      select(entity_id, lI_iqr_2)
    
    SISAL_IQR_2 <- sc_new_new_2 %>% ungroup() %>% group_by(entity_id) %>% arrange(., depth_sample, .by_group =T) %>%
      mutate(lR_iqr = (lin_reg_age_uncert_pos+lin_reg_age_uncert_neg)/2 , 
             lI_iqr =  (lin_interp_age_uncert_pos+lin_interp_age_uncert_neg)/2, 
             bacon_iqr =  (bacon_age_uncert_pos+bacon_age_uncert_neg)/2, 
             bchron_iqr =  (bchron_age_uncert_pos+bchron_age_uncert_neg)/2, 
             StalAge_iqr =  (StalAge_age_uncert_pos+StalAge_age_uncert_neg)/2)
    
    SISAL_young_2 <- SISAL_IQR_2 %>% group_by(entity_id) %>% mutate(general = n()/3) %>%
      slice(.,1:general) %>% group_by(entity_id) %>% summarise(lR_young = mean(lR_iqr, na.rm = T),
                                                               lI_young = mean(lI_iqr, na.rm = T),
                                                               bacon_young = mean(bacon_iqr, na.rm = T),
                                                               bchron_young = mean(bchron_iqr, na.rm = T),
                                                               stalage_young = mean(StalAge_iqr, na.rm = T))
    
    SISAL_old_2 <- SISAL_IQR_2 %>% group_by(entity_id) %>% mutate(general = n()/3) %>%
      slice(.,general*2:n()) %>% group_by(entity_id) %>% summarise(lR_old = mean(lR_iqr, na.rm = T),
                                                                   lI_old = mean(lI_iqr, na.rm = T),
                                                                   bacon_old = mean(bacon_iqr, na.rm = T),
                                                                   bchron_old = mean(bchron_iqr, na.rm = T),
                                                                   stalage_old = mean(StalAge_iqr, na.rm = T))
    
    SISAL_IQR_check2_2 <- left_join(SISAL_young_2, SISAL_old_2, by = "entity_id") %>% mutate(lR_iqr_2 = if_else(lR_young < lR_old , 1, 0),
                                                                                       lI_iqr_2 = if_else(lI_young < lI_old, 1,0),
                                                                                       bacon_iqr_2 = if_else(bacon_young < bacon_old,1,0),
                                                                                       bchron_iqr_2 = if_else(bchron_young < bchron_old,1,0),
                                                                                       stalage_iqr_2 = if_else(stalage_young < stalage_old, 1,0)) %>%
      mutate(lR_iqr_2 = if_else(is.na(lR_young) & is.na(lR_old), 0, lR_iqr_2),
             lI_iqr_2 = if_else(is.na(lI_young) & is.na(lI_old), 0, lI_iqr_2),
             bacon_iqr_2 = if_else(is.na(bacon_young) & is.na(bacon_old), 0, bacon_iqr_2),
             bchron_iqr_2 = if_else(is.na(bchron_young) & is.na(bchron_old),0,bchron_iqr_2),
             stalage_iqr_2 = if_else(is.na(stalage_young) & is.na(stalage_old),0,stalage_iqr_2)) %>%
      select(-lI_iqr_2)
    
    SISAL_IQR_check2 <- left_join(SISAL_IQR_check2_1, SISAL_IQR_check2_2, by = 'entity_id')
    
    data4 <- SISAL_IQR_check2 %>% select(entity_id, lR_iqr_2, lI_iqr_2, bacon_iqr_2, bchron_iqr_2, stalage_iqr_2)
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
      scale_fill_manual(guide = guide_legend(title = 'Fit yes-no'),values = c('0'='indianred', '1'='skyblue', 'NA' = 'green'))
    mine.heatmap4
    
    SISAL_eval <- left_join(SISAL_eval_rev, SISAL_eval_fit, by = "entity_id") %>% left_join(., SISAL_IQR_check, by = "entity_id") %>%
      left_join(., SISAL_IQR_check2, by = "entity_id") %>%
      mutate(lR = if_else(lR_reversal == -1, -1, if_else(lR_reversal == 0, 0, 1 + lR_fitness + lR_iqr_1 +lR_iqr_2)),
            lI = if_else(lI_reversal == -1,-1,if_else(lI_reversal == 0, 0, 1 + lI_fitness + lI_iqr_1 + lI_iqr_2)),
            bacon = if_else(bacon_reversal == -1, -1, if_else(bacon_reversal == 0, 0, 1 + bacon_fitness + bacon_iqr_1 + bacon_iqr_2)),
            bchron = if_else(bchron_reversal == -1, -1, if_else(bchron_reversal == 0, 0, 1+ bchron_fitness + bchron_iqr_1 + bchron_iqr_2)),
            stalage = if_else(stalage_reversal == -1,-1, if_else(stalage_reversal == 0, 0, 1 + stalage_fitness + stalage_iqr_1 + stalage_iqr_2)))
    
    data5 <- SISAL_eval %>% select(entity_id, lR, lI, bacon, bchron, stalage)
    names(data5) <- c('entity_id', 'lin Reg.', 'lin. Interp.', 'Bacon', 'Bchron', 'StalAge')
    mine.data5 <- gather(data = data5, key = Class, value = Fit, c(-1))
    mine.data5 <- as.data.frame(mine.data5)
    mine.data5$Fit <- factor(mine.data5$Fit, levels = c('-1', '0', '1', '2', '3', '4', 'NA'))
    
    mine.heatmap5 <- ggplot(data = mine.data5, mapping = aes(x = Class, y = factor(entity_id), fill = Fit)) + 
      geom_tile()  + 
      facet_grid( ~ Class,switch = 'x',scales = 'free_x', space = 'free_x') + 
      theme(plot.title = element_text(hjust = 0.5),
            axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
            axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
            strip.placement = 'outside', 
            strip.background = element_rect(fill = 'white',color = 'black')) + 
      xlab(label = 'Age Model') + 
      ylab(label = 'Entity ID') + 
      scale_fill_manual(guide = guide_legend(title = 'Score'),values = c('-1' = 'black', '0'='indianred', '1'='deepskyblue1', 
                                                                         '2' = 'deepskyblue2', '3' = 'deepskyblue3', '4' = 'deepskyblue4', 'NA' ='yellow'))
    #ggtitle(label = 'Age Model Check') +    
    mine.heatmap5
    
    pdf('age-model-check.pdf', 6,4)
    mine.heatmap5
    dev.off()
    
    pdf('IQR_increase_depth.pdf', 6,4)
    mine.heatmap4
    dev.off()
    
    pdf('increase_bd.pdf', 6,4)
    mine.heatmap3
    dev.off()
    
    pdf('Reversals.pdf', 6,4)
    mine.heatmap2
    dev.off()
    
    pdf('Flexibility.pdf', 6,4)
    mine.heatmap
    dev.off()
}
