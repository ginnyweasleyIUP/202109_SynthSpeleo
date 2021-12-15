#library(plyr)
library(tidyverse)

#dating <- read.csv("~/SISAL Data/sisalv2/dating_old.csv", header = T, stringsAsFactors = F) %>% 
#  mutate_at(vars(corr_age_uncert_pos, corr_age_uncert_neg), as.numeric) %>%
#  mutate(corr_age_uncert_pos = if_else(entity_id %in% c(665,666,667,668), corr_age_uncert_pos*1000, corr_age_uncert_pos),
#         corr_age_uncert_neg = if_else(entity_id %in% c(665,666,667,668), corr_age_uncert_neg*1000, corr_age_uncert_neg))
#write.csv(dating, "~/SISAL Data/sisalv2/dating.csv", row.names = F)



setwd("/stacycode/ariana/R/SISAL_AM")
source(file = 'run_functions.R')
source(file = 'SISAL_AM_add_functions.R')
source(file = 'SISAL_AM_plot_functions.R')
source(file = "StalAge_1_0_mine.R")
source(file = 'New_methods.R')
source(file = 'lin_interp_reversal.R')

#setwd("/home/ariana/Documents/v2_90ci/")
run<- read.csv('/stacywork/ariana/v2_90ci/run2.csv', header = T)
#run <- tibble(entity_id = c(665,666,667,668))
runFile <- tibble(entity_id = run$entity_id, Bacon = rep(F, length(run$entity_id)), Bchron = rep(F, length(run$entity_id)), 
                  StalAge = rep(T, length(run$entity_id)), linInterp = rep(F, length(run$entity_id)), copRa = rep(F, length(run$entity_id)),
                  linReg = rep(F, length(run$entity_id)), 
                  working_directory = rep("/stacywork/ariana/StalAge_95CI/", length(run$entity_id)), 
                  wd = rep('/home/ariana/SISAL Data/sisalv2', length(run$entity_id))) 
#runFile <- runFile %>% filter(entity_id > 665) %>% mutate(Bchron = if_else(entity_id %in% c(115,330,331,378,432,506,510), F, T)) %>% arrange(., entity_id)  # 2: ab 115
#runFile <- runFile %>% mutate(Bchron = if_else(entity_id == 668, F, T))

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
  
  
    if(bchron){
      tryCatch({
        i = 3
        runBchron_new(working_directory, file_name)
        print('Bchron done')
        graphics.off()},
        error = function(e){err <<- append(err, paste("ERROR in Bchron:",conditionMessage(e), "\n"))
        #error = function(e){print(paste("ERROR in Bchron new:",conditionMessage(e), "\n"))
        runFile[j,i] <<- F
        bchron <<- F})
    }
  
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
  write.csv(runFile, 'runFile_v2.5.csv', row.names = F, col.names = T)
  
  print(paste(file_name, 'done'))
}

load_SISAL('', runFile) # v1b
setwd("~/Documents/v2/")
write.csv(runFile, 'runFile_final_v2.csv', row.names = F, col.names = T)