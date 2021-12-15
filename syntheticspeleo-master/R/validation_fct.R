run_SISAL_chrono <- function(speleo, working_directory, file_name, hiatus) {
  saveSynSpeleo(speleo, working_directory, file_name, hiatus)

  err <- NULL
  tryCatch({
    i = 7
    runLinReg(working_directory, file_name)
    print('LinReg done')},
    error = function(e){err <<- append(err, paste("ERROR in linReg:",conditionMessage(e), "\n"))
    #error = function(e){print(paste("ERROR in linReg:",conditionMessage(e), "\n"))
    #runFile[j,i] <<- F
    linReg <<- F})

  tryCatch({
    i = 5
    runLinInterp(working_directory, file_name)
    print('LinInterp done')},
    error = function(e){err <<- append(err, paste("ERROR in linInterp:",conditionMessage(e), "\n"))
    ##  #error = function(e){print(paste("ERROR in linInterp:",conditionMessage(e), "\n"))
    # r#unFile[j,i] <<- F
    linInterp <<- F})

  tryCatch({
    i = 6
    runcopRa(working_directory, file_name, q1 = 0.05, q2 = 0.95)
    print('copRa done')},
    error = function(e){err <<- append(err, paste("ERROR in copRa:",conditionMessage(e), "\n"))
    #error = function(e){print(paste("ERROR in copRa:",conditionMessage(e), "\n"))
    #runFile[j,i] <<- F
    copRa <<- F})

  tryCatch({
    i = 2
    runBACON(working_directory, file_name)
    print('Bacon done')},
    error = function(e){err <<- append(err, paste("ERROR in Bacon:",conditionMessage(e), "\n"))
    #error = function(e){print(paste("ERROR in Bacon:",conditionMessage(e), "\n"))
    #runFile[j,i] <<- F
    bacon <<- F})

  tryCatch({
    i = 3
    runBchron_new(working_directory, file_name)
    print('Bchron done')},
    error = function(e){err <<- append(err, paste("ERROR in Bchron:",conditionMessage(e), "\n"))
    #error = function(e){print(paste("ERROR in Bchron new:",conditionMessage(e), "\n"))
    #runFile[j,i] <<- F
    bchron <<- F})

  tryCatch({
    i = 4
    runStalAge_new(working_directory, file_name)
    print('StalAge done')},
    error = function(e){err <<- append(err, paste("ERROR in StalAge:",conditionMessage(e), "\n"))
    #error = function(e){print(paste("ERROR in StalAge:",conditionMessage(e), "\n"))
    #runFile[j,i] <<- F
    stalage <<- F})

  setwd(file.path(working_directory, file_name))
  print('wd changed')

  write.table(err, 'errors.txt', row.names = F, col.names = F)

  setwd(file.path(working_directory))

  print(paste(file_name, 'done'))
  gc()
}


saveSynSpeleo <- function(speleo, working_directory, file_name, hiatus){
  setwd(working_directory)
  dating_tb <- speleo$Dtable[,1:4]
  colnames(dating_tb) <- c('dating_id', 'depth_dating','corr_age', 'corr_age_uncert')
  dating_tb <- dating_tb %>% mutate(depth_dating_new = depth_dating/10)

  sample_tb <- speleo$proxy$proxymat[,1:5]
  colnames(sample_tb) <- c('sample_id', 'depth_sample', 'interp_age', 'd13C_measurement', 'd18O_measurement')
  sample_tb <- sample_tb %>% mutate(depth_sample_new = depth_sample/10)

  dir.create(file.path(working_directory,file_name))
  setwd(file.path(working_directory, file_name))
  write.csv(hiatus, 'hiatus.csv', row.names = F, col.names = T)
  write.csv(sample_tb, "proxy_data.csv", row.names = F, col.names = T)
  write.csv(sample_tb, "original_chronology.csv", row.names = F, col.names = T)
  write.csv(hiatus, "not_used_dates.csv", row.names = F, col.names = T)


  setwd(file.path(working_directory, file_name))
  dir.create('Bacon_runs')
  setwd(file.path(getwd(),'/Bacon_runs'))
  dir.create(file_name)
  setwd(file.path(working_directory, file_name,'Bacon_runs', file_name))

  hiatus_bacon <-hiatus %>% mutate(depth_sample_bacon = depth_sample/10)  %>% select(sample_id, depth_sample_bacon)

  write.csv(sample_tb$sample_id, 'sample_id.csv', row.names = F, col.names = F)
  write.table(sample_tb$depth_sample_new, paste(file_name,'_depths.txt', sep = ''), row.names = F, col.names = F)
  write.csv(dating_tb %>% select(dating_id, corr_age, corr_age_uncert, depth_dating_new) %>% mutate(cc = 0) %>% arrange(., depth_dating_new), paste(file_name,'.csv', sep = ''), row.names = F, col.names = T)
  write.csv(hiatus_bacon,'hiatus_bacon.csv', row.names = F, col.names = T)


  setwd(file.path(working_directory, file_name))
  dir.create('Bchron')
  setwd(file.path(getwd(),'/Bchron'))

  write.csv(dating_tb %>% select(dating_id, corr_age, corr_age_uncert, depth_dating_new) %>% mutate(thickness_new = 0.5, calib_curve_new = 'normal') %>% arrange(., depth_dating_new), 'ages.csv', row.names = F, col.names = T)
  write.csv(sample_tb %>% select(sample_id, depth_sample_new) %>% arrange(., depth_sample_new), 'depths.csv', row.names = F, col.names = T)

  setwd(file.path(working_directory, file_name))
  dir.create('StalAge')
  setwd(file.path(getwd(),'/StalAge'))

  write.csv(dating_tb %>% select(dating_id, corr_age, corr_age_uncert, depth_dating) %>% arrange(., depth_dating), 'ages.csv', row.names = F, col.names = T)
  write.csv(sample_tb %>% select(sample_id, depth_sample) %>% arrange(., depth_sample), 'depths.csv', row.names = F, col.names = T)

  setwd(file.path(working_directory, file_name))
  dir.create('linInterp')
  setwd(file.path(getwd(),'/linInterp'))

  write.csv(dating_tb %>% select(dating_id, corr_age, corr_age_uncert, depth_dating) %>% mutate(date_type = 'syn. U/Th') %>% arrange(., depth_dating), 'ages.csv', row.names = F, col.names = T)
  write.csv(sample_tb %>% select(sample_id, depth_sample) %>% arrange(., depth_sample),'depths.csv', row.names = F, col.names = T)

  setwd(file.path(working_directory, file_name))
  dir.create('copRa')
  setwd(file.path(getwd(),'/copRa'))

  write.csv(dating_tb %>% select(dating_id, corr_age, corr_age_uncert, depth_dating)%>% mutate(date_type = 'syn. U/Th') %>% arrange(., depth_dating), 'ages.csv', row.names = F, col.names = T)
  write.csv(sample_tb %>% select(sample_id, depth_sample) %>% arrange(., depth_sample),'depths.csv', row.names = F, col.names = T)

  setwd(file.path(working_directory, file_name))
  dir.create('linReg')
  setwd(file.path(getwd(),'/linReg'))

  write.csv(bind_rows(sample_tb %>% select(sample_id, depth_sample),hiatus) %>% arrange(., depth_sample), 'id.csv', row.names = F, col.names = T)
  write.csv(dating_tb %>% select(dating_id, corr_age, corr_age_uncert, depth_dating) %>% mutate(date_type = 'syn. U/Th') %>% arrange(., depth_dating), 'ages.csv', row.names = F, col.names = T)
  write.csv(sample_tb %>% select(sample_id, depth_sample) %>% arrange(., depth_sample),'depths.csv', row.names = F, col.names = T)
}

merge_synSpeleo <- function(runFile, dating_tb, hiatus_tb, working_directory){
  j <- 0
  for(i in runFile$entity_id) {
    print(i)
    j <- j+1
    y <- runFile %>% filter(entity_id == i)
    file_name <- i
    #file_name <- paste(i,"-",entity_name, sep = '')
    file_name_wd <- file_name
    file_name <- list()

    dates <- dating %>% filter(entity_id == i) #%>% select(dating_id, depth_dating, date_used, date_type, corr_age, corr_age_uncert_pos, corr_age_uncert_neg)

    setwd(file.path(working_directory, file_name_wd))
    hiatus_tb <- read.csv("hiatus.csv", header = T, colClasses = c('numeric', 'numeric'))
    proxy_tb <- read.csv('proxy_data.csv', header = T, stringsAsFactors = F)

    file_name$hiatus <- hiatus_tb
    file_name$proxy <- proxy_tb
    file_name$dating <- dates

    setwd(file.path(working_directory, file_name_wd, '/linReg'))
    if(file.exists('mc_linReg_ensemble.txt')){
      AM <- read.csv('linReg_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% distinct(sample_id, .keep_all = T)
      #AM_new <- correct_am(i,AM, rm_hiatus, rm_dt, dating_from_base, h)
      ens <- read.table('mc_linReg_ensemble.txt', header = F, colClasses = 'numeric', row.names = NULL)
      names <- colnames(ens)
      ens_new <- ens %>% as_tibble() %>% rename(sample_id = names[1]) %>% distinct(sample_id, .keep_all = T) %>% as.data.frame()
      colnames(ens_new) <- NULL
      #gr <- growth.rates(ens_new)
      #colnames(gr) <- c('lin_reg_gr','lin_reg_gr_uncert_pos', 'lin_reg_gr_uncert_neg')
      file_name$linReg$AM <- AM
      file_name$linReg$ens <- ens_new
      #file_name$linReg$gr <- gr
    } else {
      file_name$linReg$AM <- list()
      file_name$linReg$ens <- list()
      #file_name$linReg$gr <- list()
    }
    print('2')
    setwd(file.path(working_directory, file_name_wd, '/linInterp'))
    if(file.exists('mc_linInt_ensemble.txt')){
      AM <- read.csv('linInt_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% distinct(sample_id, .keep_all = T)
      ens <- read.table('mc_linInt_ensemble.txt', header = F, colClasses = 'numeric', row.names = NULL)
      names <- colnames(ens)
      ens_new <- ens %>% as_tibble() %>% rename(sample_id = names[1]) %>% distinct(sample_id, .keep_all = T) %>% as.data.frame()
      colnames(ens_new) <- NULL
      #gr <- growth.rates(ens_new)
      #colnames(gr) <- c('lin_interp_gr','lin_interp_gr_uncert_pos', 'lin_interp_gr_uncert_neg')
      file_name$linInterp$AM <- AM
      file_name$linInterp$ens <- ens_new
      #file_name$linInterp$gr <- gr
    } else {
      file_name$linInterp$AM <- list()
      file_name$linInterp$ens <- list()
      #file_name$linInterp$gr <- list()
    }
    print('3')
    setwd(file.path(working_directory, file_name_wd, '/copRa'))
    if(file.exists('mc_copRa_ensemble.txt')){
      AM <- read.csv('copRa_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% distinct(sample_id, .keep_all = T)
      ens <- read.table('mc_copRa_ensemble.txt', header = F, colClasses = 'numeric', row.names = NULL)
      names <- colnames(ens)
      ens_new <- ens %>% as_tibble() %>% rename(sample_id = names[1]) %>% distinct(sample_id, .keep_all = T) %>% as.data.frame()
      colnames(ens_new) <- NULL
      #colnames(gr) <- c('copRa_gr','copRa_gr_uncert_pos','copRa_gr_uncert_neg')
      #gr <- growth.rates(ens_new)
      #colnames(gr) <- c('copRa_gr','copRa_gr_uncert_pos', 'copRa_gr_uncert_neg')
      file_name$copRa$AM <- AM
      file_name$copRa$ens <- ens_new
      #file_name$copRa$gr <- gr
    } else {
      file_name$copRa$AM <- list()
      file_name$copRa$ens <- list()
      #file_name$copRa$gr <- list()
    }
    print('4')
    setwd(file.path(working_directory, file_name_wd, '/StalAge'))
    if(file.exists('StalAge_chronology.csv')){
      AM <- read.csv('StalAge_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric'))
      file_name$StalAge$AM <- AM
      file_name$StalAge$ens <- list()
      #file_name$StalAge$gr <- list()
    } else {
      file_name$StalAge$AM <- list()
      file_name$StalAge$ens <- list()
      #file_name$StalAge$gr <- list()
    }
    print('5')
    setwd(file.path(working_directory, file_name_wd, '/Bchron'))
    if(file.exists('bchron_ensemble.txt')){
      AM <- read.csv('bchron_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric'))
      ens <- t(read.table('bchron_ensemble.txt', header = F, colClasses = 'numeric', row.names = NULL))
      ens <- as.data.frame(cbind(AM$sample_id, proxy_tb$depth_sample, ens))
      #gr <- growth.rates(ens)
      #colnames(gr) <- c('bchron_gr','bchron_gr_uncert_pos', 'bchron_gr_uncert_neg')
      colnames(ens) <- NULL
      file_name$Bchron$AM <- AM
      file_name$Bchron$ens <- ens
      #file_name$Bchron$gr <- gr
    } else {
      file_name$Bchron$AM <- list()
      file_name$Bchron$ens <- list()
      #file_name$Bchron$gr <- list()
    }

    print('6')
    setwd(file.path(working_directory, file_name_wd, '/Bacon_runs'))
    if(file.exists('mc_bacon_ensemble.txt')){
      AM <- read.csv('bacon_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric'))
      ens <- read.table('mc_bacon_ensemble.txt', header = F, colClasses = 'numeric', row.names = NULL)
      #gr <- growth.rates(ens)
      #colnames(gr) <- c('bacon_gr','bacon_gr_uncert_pos', 'bacon_gr_uncert_neg')
      colnames(ens) <- NULL

      file_name$Bacon$AM <- AM
      file_name$Bacon$ens <- ens
      #file_name$Bacon$gr <- gr
    } else {
      file_name$Bacon$AM <- list()
      file_name$Bacon$ens <- list()
      #file_name$Bacon$gr <- list()
    }

    save(file_name, file = file.path("~/synSpeleo/",paste(file_name_wd,".RData", sep = '')))
  }
}


correct_am_synSpeleo <- function(exp.chrono, file_path, file_path_new, eID) {

  for(i in eID$entity_id) {
    print(i)
    entity_name <- i
    file_name <- get(load(file.path(file_path, paste(entity_name, '.RData', sep = ''))))
    exp.chrono.file <- exp.chrono %>% ungroup() %>% filter(entity_id == i)
    file_name$linReg$AM <- exp.chrono.file %>% select(sample_id, lin_reg_age, lin_reg_age_uncert_pos, lin_reg_age_uncert_neg)
    file_name$linInterp$AM <- exp.chrono.file %>% select(sample_id, lin_interp_age, lin_interp_age_uncert_pos, lin_interp_age_uncert_neg)
    file_name$copRa$AM <- exp.chrono.file %>% select(sample_id, copRa_age, copRa_age_uncert_pos, copRa_age_uncert_neg)
    file_name$StalAge$AM <- exp.chrono.file %>% select(sample_id, StalAge_age, StalAge_age_uncert_pos, StalAge_age_uncert_neg)
    file_name$Bacon$AM <- exp.chrono.file %>% select(sample_id, bacon_age, bacon_age_uncert_pos, bacon_age_uncert_neg)
    file_name$Bchron$AM <- exp.chrono.file %>% select(sample_id, bchron_age, bchron_age_uncert_pos, bchron_age_uncert_neg)

    save(file_name, file = file.path(file_path_new, paste(entity_name, '.RData', sep = '')))
  }
}


correct_ens_synSpeleo <- function(eID, file_path, hiatus_tb){

  for(i in eID$entity_id) {
    print(i)
    entity_name <- i
    file_name <- get(load(file.path(file_path, paste(entity_name, '.RData', sep = ''))))
    #exp.chrono.file <- exp %>% ungroup() %>% filter(entity_id == i)

    h <- hiatus_tb %>% filter(entity_id == entity_name)
    print(h)

    if(!plyr::empty(data.frame(file_name$linReg$ens))){
      if(plyr::empty(data.frame(file_name$linReg$AM))){
        file_name$linReg$ens <- list()
      } else{
        file_name$linReg$ens[which(is.na(file_name$linReg$AM$lin_reg_age)),] <- NA_real_
      }

    }

    if(!plyr::empty(data.frame(file_name$linInterp$ens))){
      if(plyr::empty(data.frame(file_name$linInterp$AM))){
        file_name$linInterp$ens <- list()
      } else{
        ens <- file_name$linInterp$ens
        dim <- dim(ens)
        hiatus <- cbind(h$dating_id,h$sampling.depths, matrix(NA, nrow=length(h$dating_id),ncol=dim[2]-2))
        ens_new <- rbind(as.matrix(ens),hiatus)
        ens_new <- ens_new[order(ens_new[,2]),]
        file_name$linInterp$ens <- ens_new
        file_name$linInterp$ens[which(is.na(file_name$linInterp$AM$lin_interp_age)),] <- NA_real_
      }
    }

    if(!plyr::empty(data.frame(file_name$copRa$ens))){
      if(plyr::empty(data.frame(file_name$copRa$AM))){
        file_name$copRa$ens <- list()
      } else{
        ens <- file_name$copRa$ens
        dim <- dim(ens)
        hiatus <- cbind(h$dating_id,h$sampling.depths, matrix(NA, nrow=length(h$dating_id),ncol=dim[2]-2))
        ens_new <- rbind(as.matrix(ens),hiatus)
        ens_new <- ens_new[order(ens_new[,2]),]
        file_name$copRa$ens <- ens_new
        file_name$copRa$ens[which(is.na(file_name$copRa$AM$copRa_age)),] <- NA_real_
      }
    }

    if(!plyr::empty(data.frame(file_name$Bacon$ens))){
      if(plyr::empty(data.frame(file_name$Bacon$AM))){
        file_name$Bacon$ens <- list()
      } else{
        file_name$Bacon$ens[which(is.na(file_name$Bacon$AM$bacon_age)),] <- NA_real_
      }
    }

    if(!plyr::empty(data.frame(file_name$Bchron$ens))){
      if(plyr::empty(data.frame(file_name$Bchron$AM))){
        file_name$Bchron$ens <- list()
      } else{
        ens <- file_name$Bchron$ens
        dim <- dim(ens)
        hiatus <- cbind(h$dating_id,h$sampling.depths, matrix(NA, nrow=length(h$dating_id),ncol=dim[2]-2))
        ens_new <- rbind(as.matrix(ens),hiatus)
        ens_new <- ens_new[order(ens_new[,2]),]
        file_name$Bchron$ens <- ens_new
        file_name$Bchron$ens[which(is.na(file_name$Bchron$AM$bchron_age)),] <- NA_real_
      }
    }

    save(file_name, file = file.path(file_path, paste(entity_name, '.RData', sep = '')))
  }
}

merge_synSpeleo_chrono <- function(runFile, working_directory, se){

  SISAL_chronology_new <- data.frame(entity_id = se$entity_id,
                                     sample_id = se$sample_id,
                                     depth_sample = se$depth_sample,
                                     lin_reg_age = rep(NA_real_, length(se$sample_id)),
                                     lin_reg_age_uncert_pos = rep(NA_real_, length(se$sample_id)),
                                     lin_reg_age_uncert_neg = rep(NA_real_, length(se$sample_id)),
                                     lin_interp_age = rep(NA_real_, length(se$sample_id)),
                                     lin_interp_age_uncert_pos = rep(NA_real_, length(se$sample_id)),
                                     lin_interp_age_uncert_neg = rep(NA_real_, length(se$sample_id)),
                                     copRa_age = rep(NA_real_, length(se$sample_id)),
                                     copRa_age_uncert_pos = rep(NA_real_, length(se$sample_id)),
                                     copRa_age_uncert_neg = rep(NA_real_, length(se$sample_id)),
                                     StalAge_age = rep(NA_real_, length(se$sample_id)),
                                     StalAge_age_uncert_pos = rep(NA_real_, length(se$sample_id)),
                                     StalAge_age_uncert_neg = rep(NA_real_, length(se$sample_id)),
                                     bacon_age = rep(NA_real_, length(se$sample_id)),
                                     bacon_age_uncert_pos = rep(NA_real_, length(se$sample_id)),
                                     bacon_age_uncert_neg = rep(NA_real_, length(se$sample_id)),
                                     bchron_age = rep(NA_real_, length(se$sample_id)),
                                     bchron_age_uncert_pos = rep(NA_real_, length(se$sample_id)),
                                     bchron_age_uncert_neg = rep(NA_real_, length(se$sample_id))
  )

  for(i in runFile$entity_id) {
    print(i)
    entity_name <- i
    file_name <- i
    #y <- runFile$add

    setwd(file.path(working_directory, file_name, '/linReg'))
    if(file.exists('linReg_chronology.csv')){
      lR_chrono <- read.csv('linReg_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% mutate(sample_id = sample_id + y)  %>% distinct(sample_id, .keep_all = T)
      if(i == 'class2') {lR_chrono <- lR_chrono %>% filter(!(sample_id %in% c(350,450,451))) %>% distinct(sample_id, .keep_all = T)}
      names(lR_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
      SISAL_chronology_new <- full_join(lR_chrono, SISAL_chronology_new, by = 'sample_id') %>%
        mutate(lin_reg_age = if_else(sample_id %in% lR_chrono$sample_id, age, lin_reg_age),
               lin_reg_age_uncert_pos = if_else(sample_id %in% lR_chrono$sample_id, uncert_pos, lin_reg_age_uncert_pos),
               lin_reg_age_uncert_neg = if_else(sample_id %in% lR_chrono$sample_id, uncert_neg, lin_reg_age_uncert_neg)) %>%
        select(-age, -uncert_pos, -uncert_neg)
    }

    setwd(file.path(working_directory, file_name, '/linInterp'))
    if(file.exists('linInt_chronology.csv')){
      lI_chrono <- read.csv('linInt_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% mutate(sample_id = sample_id + y) %>% distinct(sample_id, .keep_all = T)
      names(lI_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
      SISAL_chronology_new <- full_join(lI_chrono, SISAL_chronology_new, by = 'sample_id')  %>%
        mutate(lin_interp_age = if_else(sample_id %in% lI_chrono$sample_id, age, lin_interp_age),
               lin_interp_age_uncert_pos = if_else(sample_id %in% lI_chrono$sample_id, uncert_pos, lin_interp_age_uncert_pos),
               lin_interp_age_uncert_neg = if_else(sample_id %in% lI_chrono$sample_id, uncert_neg, lin_interp_age_uncert_neg)) %>%
        select(-age, -uncert_pos, -uncert_neg)
    }

    setwd(file.path(working_directory, file_name, '/copRa'))
    if(file.exists('copRa_chronology.csv')){
      copRa_chrono <- read.csv('copRa_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% mutate(sample_id = sample_id + y) %>% distinct(sample_id, .keep_all = T)
      names(copRa_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
      SISAL_chronology_new <- full_join(copRa_chrono, SISAL_chronology_new, by = 'sample_id')  %>%
        mutate(copRa_age = if_else(sample_id %in% copRa_chrono$sample_id, age, copRa_age),
               copRa_age_uncert_pos = if_else(sample_id %in% copRa_chrono$sample_id, uncert_pos, copRa_age_uncert_pos),
               copRa_age_uncert_neg = if_else(sample_id %in% copRa_chrono$sample_id, uncert_neg, copRa_age_uncert_neg)) %>%
        select(-age, -uncert_pos, -uncert_neg)
    }


    setwd(file.path(working_directory, file_name, '/Bchron'))
    if(file.exists('bchron_chronology.csv')){
      bchron_chrono <- read.csv('bchron_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% mutate(sample_id = sample_id + y)
      names(bchron_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
      SISAL_chronology_new <- full_join(bchron_chrono, SISAL_chronology_new, by = 'sample_id') %>%
        mutate(bchron_age = if_else(sample_id %in% bchron_chrono$sample_id, age, bchron_age),
               bchron_age_uncert_pos = if_else(sample_id %in% bchron_chrono$sample_id, uncert_pos, bchron_age_uncert_pos),
               bchron_age_uncert_neg = if_else(sample_id %in% bchron_chrono$sample_id, uncert_neg, bchron_age_uncert_neg)) %>%
        select(-age, -uncert_pos, -uncert_neg)
    }


    setwd(file.path(working_directory, file_name, '/StalAge'))
    if(file.exists('StalAge_chronology.csv')){
      stalage_chrono <- read.csv('StalAge_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% mutate(sample_id = sample_id + y)
      names(stalage_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
      SISAL_chronology_new <- full_join(stalage_chrono, SISAL_chronology_new, by = 'sample_id')  %>%
        mutate(StalAge_age = if_else(sample_id %in% stalage_chrono$sample_id, age, StalAge_age),
               StalAge_age_uncert_pos = if_else(sample_id %in% stalage_chrono$sample_id, uncert_pos, StalAge_age_uncert_pos),
               StalAge_age_uncert_neg = if_else(sample_id %in% stalage_chrono$sample_id, uncert_neg, StalAge_age_uncert_neg)) %>%
        select(-age, -uncert_pos, -uncert_neg)
    }

    setwd(file.path(working_directory, file_name, '/Bacon_runs'))
    if(file.exists('bacon_chronology.csv')){
      bacon_chrono <- read.csv('bacon_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% mutate(sample_id = sample_id + y)
      if(i == 'class2') {bacon_chrono <- bacon_chrono %>% filter(!(sample_id %in% c(350,450,451))) %>% distinct(sample_id, .keep_all = T)}
      names(bacon_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
      SISAL_chronology_new <- full_join(bacon_chrono, SISAL_chronology_new, by = 'sample_id') %>%
        mutate(bacon_age = if_else(sample_id %in% bacon_chrono$sample_id, age, bacon_age),
               bacon_age_uncert_pos = if_else(sample_id %in% bacon_chrono$sample_id, uncert_pos, bacon_age_uncert_pos),
               bacon_age_uncert_neg = if_else(sample_id %in% bacon_chrono$sample_id, uncert_neg, bacon_age_uncert_neg)) %>%
        select(-age, -uncert_pos, -uncert_neg)
    }
  }
  return(SISAL_chronology_new)
}
