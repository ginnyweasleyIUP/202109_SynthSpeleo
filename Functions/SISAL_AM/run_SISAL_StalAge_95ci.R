library(tidyverse)

install.packages('Hmisc')

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
  
  write.table(err, 'errors.txt', row.names = F, col.names = F)
  
  setwd(file.path(working_directory))
  write.csv(runFile, 'runFile_v2.5.csv', row.names = F, col.names = T)
  
  print(paste(file_name, 'done'))
}

load_SISAL('', runFile) # v1b
setwd("/stacywork/ariana/StalAge_95CI/")
write.csv(runFile, 'runFile_final_v2.csv', row.names = F, col.names = T)