library(plyr)
library(tidyverse)


setwd("/home/ariana/R/SISAL_AM")
source(file = 'run_functions.R')
source(file = 'SISAL_AM_add_functions.R')
source(file = 'SISAL_AM_plot_functions.R')
source(file = "StalAge_1_0_mine.R")
source(file = 'New_methods.R')
source(file = 'lin_interp_reversal.R')

run<- read.csv('/home/ariana/COPRA_vs_copRa/am/run.csv', header = T)

runFile <- tibble(entity_id = run$entity_id, Bacon = rep(F, length(run$entity_id)), Bchron = rep(F, length(run$entity_id)), 
                  StalAge = rep(F, length(run$entity_id)), linInterp = rep(F, length(run$entity_id)), copRa = rep(T, length(run$entity_id)),
                  linReg = rep(F, length(run$entity_id)), 
                  working_directory = rep("/home/ariana/COPRA_vs_copRa/am_new/", length(run$entity_id)), 
                  wd = rep('/home/ariana/SISAL Data/sisalv2', length(run$entity_id))) 

runFile2 <- tibble(entity_id = 175, Bacon = F, Bchron = F, 
                  StalAge = F, linInterp = F, copRa = T,
                  linReg = F, 
                  working_directory = "/home/ariana/COPRA_vs_copRa/am_new/", 
                  wd = '/home/ariana/SISAL Data/sisalv2') 

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
  if(linReg){
    tryCatch({
      i = 7
      runLinReg(working_directory, file_name)
      print('LinReg done')
      graphics.off()},
      error = function(e){err <<- append(err, paste("ERROR in linReg:",conditionMessage(e), "\n"))
      #error = function(e){print(paste("ERROR in linReg:",conditionMessage(e), "\n"))
      runFile[j,i] <<- F
      linReg <<- F})
  }
  
  if(linInterp){
    tryCatch({
      i = 5
      runLinInterp(working_directory, file_name)
      print('LinInterp done')
      graphics.off()},
      error = function(e){err <<- append(err, paste("ERROR in linInterp:",conditionMessage(e), "\n"))
      #error = function(e){print(paste("ERROR in linInterp:",conditionMessage(e), "\n"))
      runFile[j,i] <<- F
      linInterp <<- F})
  }
  
  if(copRa){
    tryCatch({
      i = 6
      runcopRa(working_directory, file_name, q1 = 0.025, q2 = 0.975)
      print('copRa done')
      graphics.off()},
      error = function(e){err <<- append(err, paste("ERROR in copRa:",conditionMessage(e), "\n"))
      #error = function(e){print(paste("ERROR in copRa:",conditionMessage(e), "\n"))
      runFile[j,i] <<- F
      copRa <<- F})
  }
  
  if(bacon){
    tryCatch({
      i = 2
      runBACON(working_directory, file_name)
      print('Bacon done')
      graphics.off()},
      error = function(e){err <<- append(err, paste("ERROR in Bacon:",conditionMessage(e), "\n"))
      #error = function(e){print(paste("ERROR in Bacon:",conditionMessage(e), "\n"))
      runFile[j,i] <<- F
      bacon <<- F})
  }
  
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
  
  if(stalage){
    tryCatch({
      i = 4
      runStalAge_new(working_directory, file_name)
      print('StalAge done')
      graphics.off()},
      error = function(e){err <<- append(err, paste("ERROR in StalAge:",conditionMessage(e), "\n"))
      #error = function(e){print(paste("ERROR in StalAge:",conditionMessage(e), "\n"))
      runFile[j,i] <<- F
      stalage <<- F})
  }
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
  write.csv(runFile, 'runFile_v2.csv', row.names = F, col.names = T)
  
  print(paste(file_name, 'done'))
}

load_SISAL('', runFile) # v1b