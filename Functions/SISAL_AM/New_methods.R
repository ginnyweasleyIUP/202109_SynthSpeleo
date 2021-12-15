convert_from_base <- function(original_tb, dating_tb){
  depth_sample <- original_tb$depth_sample
  depth_dating <- dating_tb$depth_dating
  
  depth_sample_conv <- data.frame(depth_new = max(depth_sample) - depth_sample)
  depth_dating_conv <- data.frame(depth_new = max(depth_sample) - depth_dating)
  
  original_tb_new <- bind_cols(original_tb, depth_sample_conv) %>% select(-depth_sample)
  dating_tb_new <- bind_cols(original_tb, depth_dating_conv) %>% select(-depth_dating)
  
  return(list(original_tb_new, dating_tb_new))
}

plot_AM_anonymous <- function(m, folder, dating_tb, hiatus, age_median, depth_eval, age_sd_high, ag_sd_low, unknown) {
  matplot(x = age_median, y= depth_eval, col = 'black', lty = 1, type = 'l', lwd = 1, 
          ylim = c(max(depth_eval),0), xlab = 'Age [yrs BP]', ylab = 'Depth from top [mm]')
  lines(x=age_sd_high,y=depth_eval, lty = 2, col = 'red')
  lines(x=age_sd_low,y = depth_eval, lty = 2, col = 'red')
  points(x = dating_tb$corr_age, y=dating_tb$depth_dating, lty = 2, col = 'orange', pch = 4)
  arrows(dating_tb$corr_age-dating_tb$corr_age_uncert, dating_tb$depth_dating, dating_tb$corr_age+dating_tb$corr_age_uncert, dating_tb$depth_dating, length=0.05, angle=90, code=3, col = 'orange')
  
  if (!empty(data.frame(hiatus))) {
    abline(h = hiatus$depth_sample, col = 'grey', lty = 2)
  }
  
  if (length(unknown) != 0) {
    points(unknown[,1],unknown[,3], pch = 25, col ="red", cex = 2)
  }
  
  
  
}

mc_ages_sep <- function(linReg = F, linInterp = F, age, age_error, N, c14 = F, filter = T) { # N number of MC simulations
  #set.seed(2019-02-04)
  library(clam)
  
  number <- 0
  d <- 0
  age_ensemble_final <- NA
  
  idx_rev <- reversal_filter(age, filter)
  print(c('Reversals:',idx_rev))
  
  if(linReg) {
    if(c14){
      age_ensemble <- apply(cbind(age,age_error), 1, function(x)
        sample(calibrate(x[1], x[2], postbomb = 4, cc = 3, graph = FALSE, yrsteps = 0.05)$calib[,1],N, replace =F, 
               prob = calibrate(x[1], x[2], postbomb = 4, cc = 3, graph = FALSE, yrsteps = 0.05)$calib[,2])) 
      age_ensemble_final <- age_ensemble
      
    } else {
      
      age_ensemble <- apply(cbind(age,age_error), 1, function(x) rnorm(n = N, mean = x[1], sd = x[2])) # calculate N deviates for each age 
      age_ensemble_final <- age_ensemble
      
    }
    return(age_ensemble_final)
  }
  
  
  if(linInterp){
    while(d < N){
      #print('hello')
      k<- N-d
      #print(k)
      
      number <- number+1
      #print(n)
      
      if(c14){
        age_ensemble <- apply(cbind(age,age_error), 1, function(x)
          sample(calibrate(x[1], x[2], postbomb = 4, cc = 3, graph = FALSE, yrsteps = 0.05)$calib[,1],k, replace =F, 
                 prob = calibrate(x[1], x[2], postbomb = 4, cc = 3, graph = FALSE, yrsteps = 0.05)$calib[,2])) 
       # age_ensemble_final <- age_ensemble
        
      } else {
        
        age_ensemble <- apply(cbind(age,age_error), 1, function(x) rnorm(n = k, mean = x[1], sd = x[2])) # calculate N deviates for each age 
        #age_ensemble_linReg <- age_ensemble
        
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
          #print('k=1')
        }
        
      } else {
        age_ensemble_diff <- apply(age_ensemble, 1, diff) # each column contains the derivatives for one MC run -> N columns
        run <- T
        #if (k ==2) {print('k=2')}
        for (i in seq(1,k)) {
          if(any(age_ensemble_diff[,i] < 0 & !is.na(age_ensemble_diff[,i]))){
            #if (k ==2) {print('yes')}
            if (is.null(dim(age_ensemble)[1])) {
              run <- F
              #print('dim = 0')
            } else {
              #if (k ==2) {print('geloscht')}
              age_ensemble <- age_ensemble[-i,]
              #age_ensemble_linReg <- age_ensemble_linReg[-i,]
            }
            
          }
          
        }
        
      }
      
      if (run) {
        if (k == N){
          age_ensemble_final <- age_ensemble
          #age_ensemble_linReg_final <- age_ensemble_linReg
          #if(dim(age_ensemble)[1]==dim(age_ensemble_linReg)[1]){print('yes')}else{print('no')}
        } else {
          age_ensemble_final <- rbind(age_ensemble_final, age_ensemble)
          #age_ensemble_linReg_final <- rbind(age_ensemble_linReg_final, age_ensemble_linReg)
          #if(dim(age_ensemble)[1]==dim(age_ensemble_linReg)[1]){print('yes')}else{print('no')}
        }
      }
      
      run <- T
      d <- dim(age_ensemble_final)[1]
      
      
      if(number>10 && d > 1500){return(age_ensemble_final)
        break
      } else if(number > 10 && d < 1500) {stop('ERROR: too many iterations, check data!')}
      #print(c('Dim lI:', dim(age_ensemble_linInt_final)[1]))
      #print(c('Dim lR:', dim(age_ensemble_linReg_final)[1]))
    }
    return(age_ensemble_final) # ensemble
  }
  
}


runBchron_new <- function(working_directory,file_name){ # sample ids missing
  library(Bchron)
  
  setwd(file.path(working_directory, file_name, '/Bchron'))
  dating_tb <- read.csv('ages.csv', header = T, stringsAsFactors = F)
  depth_sample <- read.csv('depths.csv', header = T, stringsAsFactors = F, colClasses = c('numeric', 'numeric'))
  
  depth_eval <- depth_sample$depth_sample # in cm !!!!!!
  
  setwd(file.path(working_directory, file_name))
  hiatus <- read.csv('hiatus.csv', header = T, stringsAsFactors = F,  colClasses = c('numeric', 'numeric'))
  unknown_age <- read.csv('not_used_dates.csv', header = T)
  
  #dating_tb <- dating_tb %>% add_row(.,dating_id =1,corr_age =3000,corr_age_uncert= 60,depth_dating_new= 0,thickness_new =0.1,calib_curve_new='normal', .before = 1)
  if (!empty(data.frame(hiatus))) {
    dating_tb <- add_hiatus(dating_tb, hiatus, bchron = T)
    write.csv(dating_tb, 'hiatus_dates_bchron.csv', row.names = F, sep = ',')
  }
  
  R.utils::withTimeout(
  run <- Bchronology(ages=dating_tb$corr_age,
                     ageSds=dating_tb$corr_age_uncert,
                     positions = dating_tb$depth_dating_new,
                     positionThicknesses = dating_tb$thickness_new,
                     calCurves = dating_tb$calib_curve_new,
                     ids = dating_tb$dating_id,
                     predictPositions = depth_eval),#,
                    # Janica removed jitterPositions 13.11.21 at 02:58
                     #jitterPositions = T),
  timeout = 900,onTimeout = 'error')
  
  mcmc <- run$thetaPredict
  bchron_age <- apply(mcmc,2,median)
  bchron_quantile <- apply(mcmc, 2, function(x){quantile(x, probs = c(0.025,0.975), na.rm = T)})
  
  d <- cbind(depth_eval, bchron_age, 
             bchron_age_uncert_pos = bchron_quantile[2,]-bchron_age, 
             bchron_age_uncert_neg = bchron_age - bchron_quantile[1,])
  d <- data.frame(d)
  
  if (!empty(data.frame(hiatus))) {
    hiatus_new <- hiatus %>% mutate(depth_sample = depth_sample/10)
    d <- d %>% mutate(bchron_age = if_else(depth_eval %in% hiatus_new$depth_sample, NA_real_, bchron_age),
                      bchron_age_uncert_pos = if_else(depth_eval %in% hiatus_new$depth_sample, NA_real_, bchron_age_uncert_pos),
                      bchron_age_uncert_neg = if_else(depth_eval %in% hiatus_new$depth_sample, NA_real_, bchron_age_uncert_neg))
  }
  
  setwd(file.path(working_directory, file_name, '/Bchron'))
  write.table(mcmc, 'bchron_ensemble.txt', col.names = F, row.names = F)
  b_chrono <- data.frame(sample_id = depth_sample$sample_id, 
                         bchron_age = d[,2],
                         bchron_age_uncert_pos = d[,3],
                         bchron_age_uncert_neg = d[,4])
  write.csv(b_chrono,'bchron_chronology.csv' ,col.names = T, row.names = F)
  
  x_min <- min(b_chrono$bchron_age, na.rm = T)
  x_max <- max(b_chrono$bchron_age, na.rm = T)
  
  pdf('final_age_model.pdf', 6,4)
  matplot(x = b_chrono$bchron_age, y= d$depth_eval*10, col = 'black', lty = 1, type = 'l', lwd = 1, 
          ylim = c(max(d$depth_eval*10),0), xlab = 'Age [yrs BP]', ylab = 'Depth from top [mm]', xlim = c(3000, x_max))
  lines(x=b_chrono$bchron_age - b_chrono$bchron_age_uncert_neg,y=d$depth_eval*10, lty = 2, col = 'red')
  lines(x=b_chrono$bchron_age_uncert_pos + b_chrono$bchron_age,y = d$depth_eval*10, lty = 2, col = 'red')
  points(x = dating_tb$corr_age, y=dating_tb$depth_dating_new*10, lty = 2, col = 'orange', pch = 4)
  arrows(dating_tb$corr_age-dating_tb$corr_age_uncert, dating_tb$depth_dating_new*10, dating_tb$corr_age+dating_tb$corr_age_uncert, dating_tb$depth_dating_new*10, length=0.05, angle=90, code=3, col = 'orange')
  if (!empty(data.frame(hiatus))) {
    abline(h = hiatus$depth_sample, col = 'grey', lty = 2)
  }
  dev.off()
}

runStalAge_new <- function(working_directory, file_name){
  
  #Einlesen der Daten
  print('#-----------READ IN DATA-------------#')
  setwd(file.path(working_directory, file_name))
  hiatus <- read.csv('hiatus.csv', header = T, stringsAsFactors = F,  colClasses = c('numeric', 'numeric'))
  unknown_age <- read.csv('not_used_dates.csv', header = T)
  
  setwd(file.path(working_directory, file_name, '/StalAge'))
  Daten_orig <- read.csv("ages.csv", header = T)
  Daten_orig <- Daten_orig %>% mutate(corr_age_uncert = 2*corr_age_uncert) %>% select(dating_id, corr_age, corr_age_uncert, depth_dating) 
  
  if (!plyr::empty(data.frame(hiatus))) {
    Daten_orig <- add_hiatus(Daten_orig, hiatus, stalage=T)
    write.csv(Daten_orig, 'hiatus_dates_StalAge.csv', row.names = F, sep = ',')
  }
  
  names(Daten_orig) <- c('dating_id', 'age','error', 'depth')
  
  depths <- read.csv("depths.csv", header = T)
  names(depths) <- c('sample_id', 'dft')
  
  
  
  #StalAge Modellierung
  print('#----------SCAN DATA FOR REVERSALS ----------------#')
  Daten <- scan(Daten_orig$depth, Daten_orig$age, Daten_orig$error)
  Daten <- scan_fine(Daten)
  
  print('#----------FIT DATA ----------------#')
  fit <- age_model(Daten, depths$dft, Daten_orig, q1 = 0.025, q2 = 0.975)
  names(fit) <- c('depth_sample', 'StalAge_age','StalAge_age_uncert_pos', 'StalAge_age_uncert_neg')
  
  print('#----------WRITE DATA TO CSV FILES----------------#')
  sample_id <- depths$sample_id
  depth_eval <- depths$dft
  StalAge_age <- fit$StalAge_age
  StalAge_age_uncert_pos <- fit$StalAge_age_uncert_pos - StalAge_age
  StalAge_age_uncert_neg <- StalAge_age - fit$StalAge_age_uncert_neg
  
  d <- cbind(sample_id, depth_eval, StalAge_age, StalAge_age_uncert_pos, StalAge_age_uncert_neg)
  d <- data.frame(d)
  
  if(!plyr::empty(data.frame(hiatus))){
    d <- d %>% rowwise() %>% mutate(StalAge_age = if_else(depth_eval %in% hiatus$depth_sample, NA_real_, StalAge_age),
    StalAge_age_uncert_pos = if_else(depth_eval %in% hiatus$depth_sample, NA_real_, StalAge_age_uncert_pos),
                                    StalAge_age_uncert_neg = if_else(depth_eval %in% hiatus$depth_sample, NA_real_, StalAge_age_uncert_neg))
  }
  
  
  #chronology <- cbind(sample_id,StalAge_age, StalAge_age_uncert_pos - StalAge_age, StalAge_age - StalAge_age_uncert_neg)
  #colnames(chronology) <- c('sample_id','StalAge_age', 'StalAge_age_uncert_pos', 'StalAge_age_uncert_neg')
  
  setwd(file.path(working_directory, file_name, '/StalAge'))
  write.csv(d, file="StalAge_results.csv", 
            col.names = T, row.names = F)
  write.table(d[,c(1,3,4,5)], file="StalAge_chronology.csv", col.names = T,sep = ',' ,row.names = F)
  write.csv(Daten, file = 'StalAge_date_used.csv', col.names = T, row.names = F)
  
  pdf('final_age_model.pdf', 6,4)
  matplot(x = StalAge_age, y= depth_eval, col = 'black', lty = 1, type = 'l', lwd = 1, 
          ylim = c(max(depth_eval),0), xlab = 'Age [yrs BP]', ylab = 'Depth from top [mm]')
  lines(x=StalAge_age_uncert_pos+StalAge_age,y=depth_eval, lty = 2, col = 'red')
  lines(x=StalAge_age - StalAge_age_uncert_neg,y = depth_eval, lty = 2, col = 'red')
  points(x = Daten_orig$age, y=Daten_orig$depth, lty = 2, col = 'orange', pch = 4)
  arrows(Daten_orig$age-Daten_orig$error, Daten_orig$depth, Daten_orig$age+Daten_orig$error, Daten_orig$depth, length=0.05, angle=90, code=3, col = 'orange')
  if (!empty(data.frame(hiatus))) {
    abline(h = hiatus$depth_sample, col = 'grey', lty = 2)
  }
  dev.off()
  
}

runcopRa <- function(working_directory, file_name, N = 2000, c14 = F, filter = T, q1, q2){
  
  print('---------------- Read in data -------------')
  setwd(file.path(working_directory, file_name, '/copRa'))
  dating_tb <- read.csv('ages.csv', header = T, stringsAsFactors = F)
  depth_sample <- read.csv('depths.csv', header = T, stringsAsFactors = F, colClasses = c('numeric', 'numeric'))
  
  setwd(file.path(working_directory, file_name))
  hiatus <- read.csv('hiatus.csv', header = T, stringsAsFactors = F,  colClasses = c('numeric', 'numeric'))
  unknown_age <- read.csv('not_used_dates.csv', header = T)
  
  sample <- data.frame(sample_id = depth_sample$sample_id, depth_eval = depth_sample$depth_sample)
  
  if (!empty(data.frame(hiatus))) {
    dating_tb <- add_hiatus(dating_tb, hiatus, lininterp =T)
    write.csv(dating_tb, 'hiatus_dates_copRa.csv', row.names = F, sep = ',')
  }
  
  dating_tb <- scan_lin_interp(dating_tb)
  dating_tb <- scan_fine_lin_interp(dating_tb)
  
  setwd(file.path(working_directory, file_name, '/copRa'))
  write.csv(dating_tb, 'new_dating_tb.csv', row.names = F)
  
  print('---------------Reversal check--------------')
  if(any(dating_tb$date_type == 'C14')) { c14 = T} else {c14 = F}
  #if(any(duplicated(reversal_filter(dating_tb$corr_age, filter = F)))) {filter = T} else {filter = F}
  #filter <- T
  print('------- MC Simulations----------')
  mc_runs <- mc_ages_sep_new(linInterp = T, age = dating_tb$corr_age, age_error = dating_tb$corr_age_uncert, N = 2000, c14 = c14, working_directory = working_directory, file_name = file_name)
  
  print('------------- copRa -----------')
  N <- dim(mc_runs)[1]
  copRa <- mc_copRa(N, hiatus$depth_sample, dating_tb$depth_dating, mc_runs, depth_sample$depth_sample)
  copRa <- merge(sample, copRa, by = 'depth_eval', all.x = T, all.y = T)
  copRa <- select(copRa, sample_id, depth_eval, everything())
  
  print('---------------save data copRa --------------')
  setwd(file.path(working_directory, file_name, '/copRa'))
  write.table(as.matrix(mc_runs[[2]]),'dating_mc_copRa_ensemble.txt', col.names = F, row.names = F)
  write.table(as.matrix(copRa),'mc_copRa_ensemble.txt', col.names = F, row.names = F)
  
  print('--------------median and quantiles--------------')
  stats <- get_median_quantiles_copRa(copRa[,3:N+2], q1,q2)
  age_median <- stats[,1]
  age_sd_low <- stats[,2]
  age_sd_high <- stats[,3]
  COPrA <- cbind(copRa$sample_id, age_median, age_sd_high-age_median, age_median - age_sd_low)
  colnames(COPrA) <- c('sample_id', 'copRa_age', 'copRa_uncert_pos', 'copRa_uncert_neg')
  write.csv(COPrA, 'copRa_chronology.csv', row.names = F, sep = ',')
  
  pdf('final_age_model.pdf', 6,4)
  matplot(x = age_median, y= copRa$depth_eval, col = 'black', lty = 1, type = 'l', lwd = 1, 
          ylim = c(max(copRa$depth_eval),0), xlab = 'Age [yrs BP]', ylab = 'Depth from top [mm]')
  lines(x=age_sd_high,y=copRa$depth_eval, lty = 2, col = 'red')
  lines(x=age_sd_low,y = copRa$depth_eval, lty = 2, col = 'red')
  points(x = dating_tb$corr_age, y=dating_tb$depth_dating, lty = 2, col = 'orange', pch = 4)
  arrows(dating_tb$corr_age-dating_tb$corr_age_uncert, dating_tb$depth_dating, dating_tb$corr_age+dating_tb$corr_age_uncert, dating_tb$depth_dating, length=0.05, angle=90, code=3, col = 'orange')
  if (!empty(data.frame(hiatus))) {
    abline(h = hiatus$depth_sample, col = 'grey', lty = 2)
  }
  dev.off()
}

runBACON_new <- function(working_directory, file_name, postbomb = 0, cc = 0 ){
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
  
  j <- 10000
  tho <- c()
  
  
  setwd(file.path(working_directory, file_name))
  print('#------------ run bacon ---------------#')
  if (dim(hiatus)[1] == 0) {
    tryCatch({Bacon(core = file_name, depths.file = TRUE, thick = thickness, acc.mean = accMean, postbomb = postbomb, cc = cc, suggest = F, ask = F, ssize = j, th0 = tho)},
             error = function(e){write.table(x=paste("ERROR in Bacon:",conditionMessage(e)),file = 'bacon_error.txt')})
  } else {
    tryCatch({Bacon( core = file_name, depths.file = TRUE, thick = thickness, acc.mean = accMean ,postbomb = postbomb, hiatus.depths = hiatus$depth_sample_bacon, cc =cc, suggest = F, ask = F, ssize = j, th0 = tho)},
             error = function(e){write.table(x=paste("ERROR in Bacon:",conditionMessage(e)),file = 'bacon_error.txt')})
  }
  
  setwd(file.path(working_directory, file_name, '/Bacon_runs/',file_name))
  output <- read.table(paste(info$prefix, ".out", sep = ""))
  log_like <- output[,length(output)]
  ac <- Acf(log_like, lag.max = 1000)
  n <- which(ac$acf[,1,1]<0.3)[1] -1
  
  if(n>9){
    setwd(file.path(working_directory, file_name))
    print('#------------ run bacon ---------------#')
    if (dim(hiatus)[1] == 0) {
      tryCatch({Bacon(core = file_name, depths.file = TRUE, thick = thickness, acc.mean = accMean, postbomb = postbomb, cc = cc, suggest = F, ask = F, ssize = j, th0 = tho)},
               error = function(e){write.table(x=paste("ERROR in Bacon:",conditionMessage(e)),file = 'bacon_error.txt')})
    } else {
      tryCatch({Bacon( core = file_name, depths.file = TRUE, thick = thickness, acc.mean = accMean ,postbomb = postbomb, hiatus.depths = hiatus$depth_sample_bacon, cc =cc, suggest = F, ask = F, ssize = j, th0 = tho)},
               error = function(e){write.table(x=paste("ERROR in Bacon:",conditionMessage(e)),file = 'bacon_error.txt')})
    }
    setwd(file.path(working_directory, file_name, '/Bacon_runs/',file_name))
    output <- read.table(paste(info$prefix, ".out", sep = ""))
    log_like <- output[,length(output)]
    setwd(file.path(working_directory, file_name, '/Bacon_runs/'))
    pdf(paste('Acf_1.pdf', sep = ''), 4,6)
    ac <- Acf(log_like, lag.max = 10)
    dev.off()
    n <- which(ac$acf[,1,1]<0.3)[1] -1
    thinner(1-1/(n-1))
    
  } else {
    
    setwd(file.path(working_directory, file_name, '/Bacon_runs/'))
    pdf(paste('Acf_1.pdf', sep = ''), 4,6)
    ac <- Acf(log_like, lag.max = 10)
    dev.off()
    
    thinner(1-1/(n-1))
    
  }
  setwd(file.path(working_directory, file_name, '/Bacon_runs/',file_name))
  output <- read.table(paste(info$prefix, ".out", sep = ""))
  log_like <- output[,length(output)]
  setwd(file.path(working_directory, file_name, '/Bacon_runs/'))
  pdf(paste('Acf_2.pdf', sep = ''), 4,6)
  ac <- Acf(log_like, lag.max = 10)
  dev.off()
  
  n <- which(ac$acf[,1,1]<0.3)[1] -1
  print(n)
  
  print('--------------save data -----------------')
  bacon_mcmc <- sapply(depth_eval, Bacon.Age.d)
  bacon_age <- get_bacon_median_quantile_new(depth_eval, hiatus, bacon_mcmc)
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
