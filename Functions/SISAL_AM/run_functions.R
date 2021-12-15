#setwd("~/Documents/Masterarbeit")
#source(file = 'SISAL_AM_add_functions.R')

runLinReg <- function(working_directory, file_name, N = 2000, c14 = F, filter = T){
  
  print('---------------- Read in data -------------')
  setwd(file.path(working_directory, file_name, '/linReg'))
  dating_tb <- read.csv('ages.csv', header = T, stringsAsFactors = F)
  depth_sample <- read.csv('depths.csv', header = T, stringsAsFactors = F, colClasses = c('numeric', 'numeric'))
  id <- read.csv('id.csv', header = T, stringsAsFactors = F, colClasses = c('numeric', 'numeric'))
  
  setwd(file.path(working_directory, file_name))
  hiatus <- read.csv('hiatus.csv', header = T, stringsAsFactors = F,  colClasses = c('numeric', 'numeric'))
  unknown_age <- read.csv('not_used_dates.csv', header = T)
  
  sample <- data.frame(sample_id = id$sample_id, depth_eval = id$depth_sample)
  
  print('---------------Reversal check--------------')
  if(any(dating_tb$date_type == 'C14')) { c14 = T} else {c14 = F}
  #if(any(duplicated(reversal_filter(dating_tb$corr_age, filter = F)))) {filter = T} else {filter = F}
  filter <- T
  print('------- MC Simulations----------')
  mc_runs <- mc_ages_sep_new(linReg = T, age = dating_tb$corr_age, age_error = dating_tb$corr_age_uncert, N = 2000, c14 = c14, working_directory = working_directory, file_name = file_name)
  
  print('--------lin Reg ------------')
  N <- dim(mc_runs)[1]
  lr <- mc_linReg(N, hiatus$depth_sample, dating_tb$depth_dating, mc_runs, depth_sample$depth_sample)
  
  lr <- merge(sample, lr, by = 'depth_eval', all.x = T, all.y = T)
  lr <- select(lr, sample_id, depth_eval, everything())
  
  print('---------------save data Lin Reg --------------')
  setwd(file.path(working_directory, file_name, '/linReg'))
  write.table(as.matrix(mc_runs),'dating_mc_linReg_ensemble.txt', col.names = F, row.names = F)
  write.table(as.matrix(lr),'mc_linReg_ensemble.txt', col.names = F, row.names = F)
  
  print('--------------median and quantiles--------------')
  stats <- get_median_quantiles(lr[,3:N+2])
  age_median <- stats[,1]
  age_sd_low <- stats[,2]
  age_sd_high <- stats[,3]
  lin_reg <- cbind(lr$sample_id,age_median, age_sd_high-age_median, age_median - age_sd_low)
  colnames(lin_reg) <- c('sample_id', 'lin_reg_age', 'lin_reg_age_uncert_pos', 'lin_reg_age_uncert_neg')
  write.csv(lin_reg, 'linReg_chronology.csv',  row.names = F, sep =',')
  
  pdf('final_age_model.pdf', 6,4)
  matplot(x = age_median, y= lr$depth_eval, col = 'black', lty = 1, type = 'l', lwd = 1, 
          ylim = c(max(lr$depth_eval),0), xlab = 'Age [yrs BP]', ylab = 'Depth from top [mm]')
  lines(x=age_sd_high,y=lr$depth_eval, lty = 2, col = 'red')
  lines(x=age_sd_low,y = lr$depth_eval, lty = 2, col = 'red')
  points(x = dating_tb$corr_age, y=dating_tb$depth_dating, lty = 2, col = 'orange', pch = 4)
  arrows(dating_tb$corr_age-dating_tb$corr_age_uncert, dating_tb$depth_dating, dating_tb$corr_age+dating_tb$corr_age_uncert, dating_tb$depth_dating, length=0.05, angle=90, code=3, col = 'orange')
  if (!empty(data.frame(hiatus))) {
    abline(h = hiatus$depth_sample, col = 'grey', lty = 2)
  }
  dev.off()
}

runLinInterp <- function(working_directory, file_name, N = 2000, c14 = F, filter = T){
  
  print('---------------- Read in data -------------')
  setwd(file.path(working_directory, file_name, '/linInterp'))
  dating_tb <- read.csv('ages.csv', header = T, stringsAsFactors = F)
  depth_sample <- read.csv('depths.csv', header = T, stringsAsFactors = F, colClasses = c('numeric', 'numeric'))
  
  setwd(file.path(working_directory, file_name))
  hiatus <- read.csv('hiatus.csv', header = T, stringsAsFactors = F,  colClasses = c('numeric', 'numeric'))
  unknown_age <- read.csv('not_used_dates.csv', header = T)
  
  sample <- data.frame(sample_id = depth_sample$sample_id, depth_eval = depth_sample$depth_sample)
  
  if (!empty(data.frame(hiatus))) {
    dating_tb <- add_hiatus(dating_tb, hiatus, lininterp =T)
    write.csv(dating_tb, 'hiatus_dates_interp.csv', row.names = F, sep = ',')
  }
  
  dating_tb <- scan_lin_interp(dating_tb)
  dating_tb <- scan_fine_lin_interp(dating_tb)
  
  setwd(file.path(working_directory, file_name, '/linInterp'))
  write.csv(dating_tb, 'new_dating_tb.csv', row.names = F)
  
  print('---------------Reversal check--------------')
  if(any(dating_tb$date_type == 'C14')) { c14 = T} else {c14 = F}
  #if(any(duplicated(reversal_filter(dating_tb$corr_age, filter = F)))) {filter = T} else {filter = F}
  #filter <- T
  print('------- MC Simulations----------')
  mc_runs <- mc_ages_sep_new(linInterp = T, age = dating_tb$corr_age, age_error = dating_tb$corr_age_uncert, N = 2000, c14 = c14, working_directory = working_directory, file_name = file_name)
  
  print('------------- lin Interp -----------')
  N <- dim(mc_runs)[1]
  li <- mc_linInt(N, hiatus$depth_sample, dating_tb$depth_dating, mc_runs, depth_sample$depth_sample)
  li <- merge(sample, li, by = 'depth_eval', all.x = T, all.y = T)
  li <- select(li, sample_id, depth_eval, everything())
  
  print('---------------save data Lin Interp --------------')
  setwd(file.path(working_directory, file_name, '/linInterp'))
  write.table(as.matrix(mc_runs[[2]]),'dating_mc_linInt_ensemble.txt', col.names = F, row.names = F)
  write.table(as.matrix(li),'mc_linInt_ensemble.txt', col.names = F, row.names = F)
  
  print('--------------median and quantiles--------------')
  stats <- get_median_quantiles(li[,3:N+2])
  age_median <- stats[,1]
  age_sd_low <- stats[,2]
  age_sd_high <- stats[,3]
  lin_interp <- cbind(li$sample_id, age_median, age_sd_high-age_median, age_median - age_sd_low)
  colnames(lin_interp) <- c('sample_id', 'lin_interp_age', 'lin_interp_age_uncert_pos', 'lin_interp_age_uncert_neg')
  write.csv(lin_interp, 'linInt_chronology.csv', row.names = F, sep = ',')
  
  pdf('final_age_model.pdf', 6,4)
  matplot(x = age_median, y= li$depth_eval, col = 'black', lty = 1, type = 'l', lwd = 1, 
          ylim = c(max(li$depth_eval),0), xlab = 'Age [yrs BP]', ylab = 'Depth from top [mm]')
  lines(x=age_sd_high,y=li$depth_eval, lty = 2, col = 'red')
  lines(x=age_sd_low,y = li$depth_eval, lty = 2, col = 'red')
  points(x = dating_tb$corr_age, y=dating_tb$depth_dating, lty = 2, col = 'orange', pch = 4)
  arrows(dating_tb$corr_age-dating_tb$corr_age_uncert, dating_tb$depth_dating, dating_tb$corr_age+dating_tb$corr_age_uncert, dating_tb$depth_dating, length=0.05, angle=90, code=3, col = 'orange')
  if (!empty(data.frame(hiatus))) {
    abline(h = hiatus$depth_sample, col = 'grey', lty = 2)
  }
  dev.off()
}

runStalAge <- function(working_directory, file_name){
  
  #Einlesen der Daten
  print('#-----------READ IN DATA-------------#')
  setwd(file.path(working_directory, file_name))
  hiatus <- read.csv('hiatus.csv', header = T, stringsAsFactors = F,  colClasses = c('numeric', 'numeric'))
  
  setwd(file.path(working_directory, file_name, '/StalAge'))
  Daten_orig <- read.csv("ages.csv", header = T)
  Daten_orig <- Daten_orig %>% mutate(corr_age_uncert = corr_age_uncert*2) # 2-sigma
  
  if (!empty(data.frame(hiatus))) {
    Daten_orig <- add_hiatus(Daten_orig, hiatus)
  }
  
  names(Daten_orig) <- c('dating_id', 'age','error', 'depth')
  
  depths <- read.csv("depths.csv", header = T)
  names(depths) <- c('sample_id', 'dft')
  
  
  
  #StalAge Modellierung
  print('#----------SCAN DATA FOR REVERSALS ----------------#')
  Daten <- scan(Daten_orig$depth, Daten_orig$age, Daten_orig$error)
  Daten <- scan_fine(Daten)
  
  print('#----------FIT DATA ----------------#')
  fit <- age_model(Daten, depths$dft, Daten_orig)
  names(fit) <- c('depth_sample', 'StalAge_age','StalAge_age_uncert_pos', 'StalAge_age_uncert_neg')
  
  print('#----------WRITE DATA TO CSV FILES----------------#')
  sample_id <- depths$sample_id
  depth_eval <- depths$dft
  StalAge_age <- fit$StalAge_age
  StalAge_age_uncert_pos <- fit$StalAge_age_uncert_pos - StalAge_age
  StalAge_age_uncert_neg <- StalAge_age - fit$StalAge_age_uncert_neg
  
  d <- cbind(sample_id, depth_eval, StalAge_age, StalAge_age_uncert_pos, StalAge_age_uncert_neg)
  h <- data.frame(sample_id = hiatus$sample_id, depth_eval = hiatus$depth_sample, StalAge_age = replicate(dim(hiatus)[1],NA), 
                  StalAge_age_uncert_pos = replicate(dim(hiatus)[1],NA),
                  StalAge_age_uncert_neg = replicate(dim(hiatus)[1],NA))
  d <- rbind(d, h)
  d <- d[order(d[,2]),] 
  
  
  #chronology <- cbind(sample_id,StalAge_age, StalAge_age_uncert_pos - StalAge_age, StalAge_age - StalAge_age_uncert_neg)
  #colnames(chronology) <- c('sample_id','StalAge_age', 'StalAge_age_uncert_pos', 'StalAge_age_uncert_neg')
  
  setwd(file.path(working_directory, file_name, '/StalAge'))
  write.csv(d, file="StalAge_results.csv", 
            col.names = T, row.names = F)
  write.table(d[,c(1,3,4,5)], file="StalAge_chronology.csv", col.names = T,sep = ',' ,row.names = F)
  write.csv(Daten, file = 'StalAge_date_used.csv', col.names = T, row.names = F)
  
  pdf('final_age_model.pdf', 6,4)
  matplot(x = StalAge_age, y= depth_eval, col = 'black', lty = 1, type = 'l', lwd = 1, 
          ylim = c(max(depth_eval),0), xlab = 'Age [kyrs BP]', ylab = 'Depth from top [mm]')
  lines(x=StalAge_age_uncert_pos,y=depth_eval, lty = 2, col = 'red')
  lines(x=StalAge_age_uncert_neg,y = depth_eval, lty = 2, col = 'red')
  
  #legend("topright",pt.cex = 1.5, legend = c("C14 - used","U/Th ages - not used"), pch = c(4,25),col = c("orange", "black") ,bty = 'n', title = "AM date types", cex = 0.8)
  #legend('bottomleft', legend= c('Original','Hiatus'), col = c(2,'grey'), lwd = 2, lty = c(2,1), bty = 'n', title = "Additional AM information", cex = 0.8)
  
  #axis(side = 1, lwd = 2, lwd.ticks = 2, cex.axis = 1.1)
  #axis(side = 2,  lwd = 2, lwd.ticks = 2, cex.axis = 1.1, labels = seq(0,3000, by = 1000), at= seq(0,3000, by = 1000))
  #mtext(side = 1, text = 'Age [kyrs BP]', line = 3, cex = 1.5)
  #mtext(side = 2, text = 'Depth from top [mm]', line = 2, cex  =1.5, padj = -0.4)
  
  dev.off()
 
}

runBACON <- function(working_directory, file_name, postbomb = 0, cc = 0){
  library(rbacon)
  library(forecast)
  
  setwd(file.path(working_directory, file_name, '/Bacon_runs/', file_name))
  depth_eval <- matrix(read.table(paste(file_name,'_depths.txt', sep = ''), col.names = ''))[[1]]
  sample_id <-  read.csv('sample_id.csv', header = T, stringsAsFactors = F,  colClasses = c('numeric'))
  hiatus <- read.csv('hiatus_bacon.csv', header = T, stringsAsFactors = F,  colClasses = c('numeric', 'numeric'))
  core <- read.csv(paste(file_name,'.csv', sep = ''), header = T, stringsAsFactors = F,  colClasses = c('numeric', 'numeric', 'numeric', 'numeric'))
  
  
  accMean<- sapply(c(1, 2, 5), function(x) x * 10^(-1:2))
  ballpacc <- lm(core[,2] * 1.1 ~ core[, 4])$coefficients[2]
  ballpacc <- abs(accMean - ballpacc)
  ballpacc <- ballpacc[ballpacc > 0]
  accMean <- accMean[order(ballpacc)[1]]
  
  k <- seq(floor(min(depth_eval, na.rm = T)), ceiling(max(depth_eval, na.rm = T)), by = 5)
  if (k < 10) {
    thickness <- pretty(5 * (k/10), 10)
    thickness <- min(thickness[thickness > 0])
  } else if (k > 20) {
    thickness <- max(pretty(5 * (k/20)))
  }
  
  j <- 2000
  tho <- c()
  
  
  setwd(file.path(working_directory, file_name))
  unknown_age <- read.csv('not_used_dates.csv', header = T)
  print('#------------ run bacon ---------------#')
  if (dim(hiatus)[1] == 0) {
    tryCatch({Bacon(core = file_name, depths.file = TRUE, thick = thickness, acc.mean = accMean, postbomb = postbomb, cc = cc, suggest = F, ask = F, ssize = j, th0 = tho)},
            error = function(e){write.table(x=paste("ERROR in Bacon:",conditionMessage(e)),file = 'bacon_error.txt')})
  } else {
    tryCatch({Bacon( core = file_name, depths.file = TRUE, thick = thickness, acc.mean = accMean ,postbomb = postbomb, hiatus.depths = hiatus$depth_sample_bacon, cc =cc, suggest = F, ask = F, ssize = j, th0 = tho)},
            error = function(e){write.table(x=paste("ERROR in Bacon:",conditionMessage(e)),file = 'bacon_error.txt')})
  }
  
  print('--------------save data -----------------')
  bacon_mcmc <- sapply(depth_eval, Bacon.Age.d)
  bacon_age <- get_bacon_median_quantile(depth_eval, hiatus, bacon_mcmc)
  bacon_mcmc <- rbind(depth_eval, bacon_mcmc)
  bacon_mcmc <- t(bacon_mcmc)
  bacon_mcmc <- cbind(sample_id, bacon_mcmc)
  
  h <- cbind(hiatus, matrix(NA,nrow = dim(hiatus)[1], ncol = dim(bacon_mcmc)[2]-2))
  names(h) <- names(bacon_mcmc)
  
  bacon_mcmc <- rbind(bacon_mcmc, h)
  bacon_mcmc <- bacon_mcmc[order(bacon_mcmc[,2]),]
  
  sample_id <- bacon_mcmc[,1]
  
  setwd(file.path(working_directory, file_name,'/Bacon_runs'))
  write.table(bacon_mcmc, 'mc_bacon_ensemble.txt', col.names = F, row.names = F)
  write.csv(cbind(sample_id, bacon_age[,2:4]),'bacon_chronology.csv' ,col.names = T, row.names = F)
  
  pdf('final_age_model.pdf', 6,4)
  matplot(x = bacon_age[,2], y= bacon_age[,1]*10, col = 'black', lty = 1, type = 'l', lwd = 1,
          ylim = c(max(depth_eval*10),0), xlab = 'Age [yrs BP]', ylab = 'Depth from top [mm]')
  lines(x=bacon_age[,3]+bacon_age[,2],y=bacon_age[,1]*10, lty = 2, col = 'red')
  lines(x=bacon_age[,2]-bacon_age[,4],y = bacon_age[,1]*10, lty = 2, col = 'red')
  points(x = core$corr_age, y=core$depth_dating_new*10, lty = 2, col = 'orange', pch = 4)
  arrows(core$corr_age-core$corr_age_uncert, core$depth_dating_new*10, core$corr_age+core$corr_age_uncert, core$depth_dating_new*10, length=0.05, angle=90, code=3, col = 'orange')
  if (!empty(data.frame(hiatus))) {
    abline(h = hiatus$depth_sample_bacon*10, col = 'grey', lty = 2)
  }
  dev.off()
  
}
runBchron <- function(working_directory,file_name){ # sample ids missing
  library(Bchron)
  
  setwd(file.path(working_directory, file_name, '/Bchron'))
  dating_tb <- read.csv('ages.csv', header = T, stringsAsFactors = F)
  depth_sample <- read.csv('depths.csv', header = T, stringsAsFactors = F, colClasses = c('numeric', 'numeric'))
  
  depth_eval <- depth_sample$depth_sample
  
  setwd(file.path(working_directory, file_name))
  hiatus <- read.csv('hiatus.csv', header = T, stringsAsFactors = F,  colClasses = c('numeric', 'numeric'))
  
  run <- Bchronology(ages=dating_tb$corr_age,
                     ageSds=dating_tb$corr_age_uncert,
                     positions = dating_tb$depth_dating_new,
                     positionThicknesses = dating_tb$thickness_new,
                     calCurves = dating_tb$calib_curve_new,
                     ids = dating_tb$dating_id,
                     predictPositions = depth_eval,
                     jitterPositions = T)
  
  mcmc <- run$thetaPredict
  bchron_age <- apply(mcmc,2,median)
  bchron_quantile <- apply(mcmc, 2, function(x){quantile(x, probs = c(0.025,0.975), na.rm = T)})
  
  d <- cbind(depth_eval, bchron_age, 
             bchron_age_uncert_pos = bchron_quantile[2,]-bchron_age, 
             bchron_age_uncert_neg = bchron_age - bchron_quantile[1,])
  h <- data.frame(depth_eval = hiatus$depth_sample, bchron_age = replicate(dim(hiatus)[1],NA), 
                  bchron_age_uncert_pos = replicate(dim(hiatus)[1],NA),
                  bchron_age_uncert_neg = replicate(dim(hiatus)[1],NA))
  d <- rbind(d, h)
  d <- d[order(d[,1]),] 
  
  mcmc <- t(mcmc)
  mcmc <- cbind(depth_sample, mcmc)
  h <- cbind(hiatus, matrix(NA,nrow = dim(hiatus)[1], ncol = dim(mcmc)[2]-2))
  names(h) <- names(mcmc)
  mcmc <- rbind(mcmc, h)
  mcmc <- mcmc[order(mcmc[,2]),]
  
  sample_id <- mcmc[,1]
  
  setwd(file.path(working_directory, file_name, '/Bchron'))
  write.table(mcmc, 'bchron_ensemble.txt', col.names = F, row.names = F)
  write.csv(cbind(sample_id, d[,2:4]),'bchron_chronology.csv' ,col.names = T, row.names = F)
}
