growth.rates <- function(ens){
  dim <- dim(ens)
  
  gr <- apply(ens[,3:dim[2]], 2,function(x){
    dx <- diff(x)
    if(any(dx == 0)){
      k <- which(dx == 0)
      dx[k] <- runif(length(k))
    }
    return(diff(ens[,2])/dx)})
  
  gr_stat <-get_median_quantiles(upd = gr)
  return(data.frame(median.gr = gr_stat[,1], gr_uncert_pos = (gr_stat[,3]-gr_stat[,1]), gr_uncert_neg = (gr_stat[,1]-gr_stat[,2])))
}

merge_SISAL_ensemble_final <- function(runFile, dating){
  j <- 0
  entity <- read.csv('/home/ariana/SISAL Data/sisalv2/entity.csv', header = T, stringsAsFactors = F)
  h <- read.csv('/home/ariana/SISAL Data/sisalv2/hiatus.csv', header = T, stringsAsFactors = F)
  
  for(i in runFile$entity_id) {
    print(i)
    j <- j+1
    y <- runFile %>% filter(entity_id == i) 
    entity_name <- (entity %>% filter(entity_id == i))$entity_name
    file_name <- paste(i,"-",entity_name, sep = '')
    file_name_wd <- file_name
    file_name <- list()
    
    dates <- dating %>% filter(entity_id == i) %>% select(dating_id, depth_dating, date_used, date_type, corr_age, corr_age_uncert_pos, corr_age_uncert_neg)
    
    setwd(file.path(y$working_directory, file_name_wd))
    hiatus_tb <- read.csv("hiatus.csv", header = T, colClasses = c('numeric', 'numeric'))
    info_tb <- read.csv("info.csv", header = T, stringsAsFactors = F)
    orig_tb <- read.csv('original_chronology.csv', header = T, stringsAsFactors = F)
    proxy_tb <- read.csv('proxy_data.csv', header = T, stringsAsFactors = F)
    
    file_name$hiatus <- hiatus_tb
    file_name$info <- info_tb
    file_name$origAM <- orig_tb
    file_name$proxy <- proxy_tb
    file_name$dating <- dates
    
    setwd(file.path(y$working_directory, file_name_wd, '/linReg'))
    if(file.exists('mc_linReg_ensemble.txt')){
      AM <- read.csv('linReg_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% distinct(sample_id, .keep_all = T)
      #AM_new <- correct_am(i,AM, rm_hiatus, rm_dt, dating_from_base, h)
      ens <- read.table('mc_linReg_ensemble.txt', header = F, colClasses = 'numeric', row.names = NULL)
      names <- colnames(ens)
      ens_new <- ens %>% as_tibble() %>% rename(sample_id = names[1]) %>% distinct(sample_id, .keep_all = T) %>% as.data.frame()
      colnames(ens_new) <- NULL
      gr <- growth.rates(ens_new)
      colnames(gr) <- c('lin_reg_gr','lin_reg_gr_uncert_pos', 'lin_reg_gr_uncert_neg')
      file_name$linReg$AM <- AM
      file_name$linReg$ens <- ens_new
      file_name$linReg$gr <- gr
    } else {
      file_name$linReg$AM <- list()
      file_name$linReg$ens <- list()
      file_name$linReg$gr <- list()
      runFile[j,7] <- F
    }
    print('2')
    setwd(file.path(y$working_directory, file_name_wd, '/linInterp'))
    if(file.exists('mc_linInt_ensemble.txt')){
      AM <- read.csv('linInt_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% distinct(sample_id, .keep_all = T)
      ens <- read.table('mc_linInt_ensemble.txt', header = F, colClasses = 'numeric', row.names = NULL) 
      names <- colnames(ens)
      ens_new <- ens %>% as_tibble() %>% rename(sample_id = names[1]) %>% distinct(sample_id, .keep_all = T) %>% as.data.frame()
      colnames(ens_new) <- NULL
      gr <- growth.rates(ens_new)
      colnames(gr) <- c('lin_interp_gr','lin_interp_gr_uncert_pos', 'lin_interp_gr_uncert_neg')
      file_name$linInterp$AM <- AM
      file_name$linInterp$ens <- ens_new
      file_name$linInterp$gr <- gr
    } else {
      file_name$linInterp$AM <- list()
      file_name$linInterp$ens <- list()
      file_name$linInterp$gr <- list()
      runFile[j,5] <- F
    }
    print('3')
    setwd(file.path(y$working_directory, file_name_wd, '/copRa'))
    if(file.exists('mc_copRa_ensemble.txt')){
      AM <- read.csv('copRa_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% distinct(sample_id, .keep_all = T)
      ens <- read.table('mc_copRa_ensemble.txt', header = F, colClasses = 'numeric', row.names = NULL) 
      names <- colnames(ens)
      ens_new <- ens %>% as_tibble() %>% rename(sample_id = names[1]) %>% distinct(sample_id, .keep_all = T) %>% as.data.frame()
      colnames(ens_new) <- NULL
      colnames(gr) <- c('copRa_gr','copRa_gr_uncert_pos','copRa_gr_uncert_neg')
      gr <- growth.rates(ens_new)
      colnames(gr) <- c('copRa_gr','copRa_gr_uncert_pos', 'copRa_gr_uncert_neg')
      file_name$copRa$AM <- AM
      file_name$copRa$ens <- ens_new
      file_name$copRa$gr <- gr
    } else {
      file_name$copRa$AM <- list()
      file_name$copRa$ens <- list()
      file_name$copRa$gr <- list()
      runFile[j,6] <- F
    }
    print('4')
    setwd(file.path(y$working_directory, file_name_wd, '/StalAge'))
    if(file.exists('StalAge_chronology.csv')){
      AM <- read.csv('StalAge_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric'))
      file_name$StalAge$AM <- AM
      file_name$StalAge$ens <- list()
      file_name$StalAge$gr <- list()
    } else {
      file_name$StalAge$AM <- list()
      file_name$StalAge$ens <- list()
      file_name$StalAge$gr <- list()
      runFile[j, 4] <- F
    }
    print('5')
    if(y$Bchron){
      setwd(file.path(y$working_directory, file_name_wd, '/Bchron'))
      if(file.exists('bchron_ensemble.txt')){
        AM <- read.csv('bchron_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric'))
        ens <- t(read.table('bchron_ensemble.txt', header = F, colClasses = 'numeric', row.names = NULL))
        ens <- as.data.frame(cbind(AM$sample_id, orig_tb$depth_sample, ens))
        gr <- growth.rates(ens)
        colnames(gr) <- c('bchron_gr','bchron_gr_uncert_pos', 'bchron_gr_uncert_neg')
        colnames(ens) <- NULL
        file_name$Bchron$AM <- AM
        file_name$Bchron$ens <- ens
        file_name$Bchron$gr <- gr
      } else {
        file_name$Bchron$AM <- list()
        file_name$Bchron$ens <- list()
        file_name$Bchron$gr <- list()
        runFile[j,3] <- F
      }
    } else {
      file_name$Bchron$AM <- list()
      file_name$Bchron$ens <- list()
      file_name$Bchron$gr <- list()
      runFile[j,3] <- F
    }
    print('6')
    setwd(file.path(y$working_directory, file_name_wd, '/Bacon_runs'))
    if(file.exists('mc_bacon_ensemble.txt')){
      AM <- read.csv('bacon_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) 
      ens <- read.table('mc_bacon_ensemble.txt', header = F, colClasses = 'numeric', row.names = NULL)
      gr <- growth.rates(ens)
      colnames(gr) <- c('bacon_gr','bacon_gr_uncert_pos', 'bacon_gr_uncert_neg')
      colnames(ens) <- NULL
      
      file_name$Bacon$AM <- AM
      file_name$Bacon$ens <- ens
      file_name$Bacon$gr <- gr
    } else {
      file_name$Bacon$AM <- list()
      file_name$Bacon$ens <- list()
      file_name$Bacon$gr <- list()
      runFile[j,2] <- F
    }
    
    save(file_name, file = file.path("/home/ariana/SISAL/v2_upd_2/",paste(file_name_wd,".RData", sep = '')))
  }
  return(runFile)
}

plot_am <- function(chrono,.entity_id, dating., entity_name){
  not.age <- dating. %>% filter(entity_id == .entity_id) %>% filter(date_used == 'no' | date_used == 'unknown') %>%
    mutate(corr_age_uncert = 0.5*(corr_age_uncert_pos + corr_age_uncert_neg)) %>%
    mutate(.pch = if_else(date_type == "MC-ICP-MS U/Th", 2,
                          if_else(date_type == "TIMS", 8,
                                  if_else(date_type == "ICP-MS U/Th Other", 4,
                                          if_else(date_type == "Alpha U/Th", 5,
                                                  if_else(date_type ==  "U/Th unspecified",7, 
                                                          if_else(date_type == 'Event; start of laminations' |
                                                                    date_type == 'Event; end of laminations', 3, 
                                                                  if_else(date_type == 'C14', 12, NA_real_)))))))) 
  
  age <- dating. %>% filter(entity_id == .entity_id) %>% filter( date_used == 'yes') %>%
    mutate(corr_age_uncert = 0.5*(corr_age_uncert_pos + corr_age_uncert_neg)) %>% 
    mutate(.pch = if_else(date_type == "MC-ICP-MS U/Th", 2,
                          if_else(date_type == "TIMS", 8,
                                  if_else(date_type == "ICP-MS U/Th Other", 4,
                                          if_else(date_type == "Alpha U/Th", 5,
                                                  if_else(date_type ==  "U/Th unspecified",7, NA_real_))))))
  
  #new.age <- read.csv(file.path(file_path, paste(.entity_id,'-',entity_name, sep = ''),'/copRa/new_dating_tb.csv'), header = T) %>% 
  #  filter(corr_age_uncert != age$corr_age_uncert) %>%
  #  mutate(.pch = if_else(date_type == "MC-ICP-MS U/Th", 2,
  #                        if_else(date_type == "TIMS", 8,
  #                                if_else(date_type == "ICP-MS U/Th Other", 4,
  #                                        if_else(date_type == "Alpha U/Th", 5,
  #                                                if_else(date_type ==  "U/Th unspecified",7, NA_real_))))))
  
  color = c('black','royalblue3','forestgreen','skyblue','springgreen2','purple', 'sienna1', 'brown2')
  
  n <- rep(0,8)
  if(!plyr::empty(data.frame(chrono$origAM$interp_age))){n[1]<- 1}
  if(!plyr::empty(data.frame(chrono$copRa$AM$copRa_age))){n[2]<- 2}
  if(!plyr::empty(data.frame(chrono$StalAge$AM$StalAge_age))){n[3]<- 3}
  if(!plyr::empty(data.frame(chrono$linInterp$AM$lin_interp_age))){n[4]<- 4}
  if(!plyr::empty(data.frame(chrono$linReg$AM$lin_reg_age))){n[5]<- 5}
  if(!plyr::empty(data.frame(chrono$Bacon$AM$bacon_age))){n[6]<- 6}
  if(!plyr::empty(data.frame(chrono$Bchron$AM$bchron_age))){n[7]<- 7}
  if(!plyr::empty(data.frame(chrono$OxCal$AM$OxCal_age))){n[8]<- 8}
  
  n <- n[which(n !=0)]
  
  matplot(x = cbind(chrono$origAM$interp_age,chrono$copRa$AM$copRa_age, chrono$StalAge$AM$StalAge_age, chrono$linInterp$AM$lin_interp_age, chrono$linReg$AM$lin_reg_age, chrono$Bacon$AM$bacon_age, chrono$Bchron$AM$bchron_age, chrono$OxCal$AM$OxCal_age), 
          y= chrono$origAM$depth_sample, 
          col = color[n], lty = 1, type = 'l', lwd = 1.5, 
          xlim = range(not.age$corr_age, not.age$corr_age +not.age$corr_age_uncert, not.age$corr_age -not.age$corr_age_uncert, 
                       age$corr_age, age$corr_age +age$corr_age_uncert, age$corr_age -age$corr_age_uncert, 
                       chrono$origAM$interp_age, 
                       chrono$copRa$AM$copRa_age, 
                       chrono$StalAge$AM$StalAge_age,
                       chrono$linInterp$AM$lin_interp_age,
                       chrono$linReg$AM$lin_reg_age,
                       chrono$Bacon$AM$bacon_age,
                       chrono$Bchron$AM$bchron_age,
                       chrono$OxCal$AM$OxCal_age,
                       na.rm=TRUE),
          ylim = c(max(range(chrono$origAM$depth_sample, not.age$depth_dating, age$depth_dating, na.rm = T), na.rm = T),min(range(chrono$origAM$depth_sample, not.age$depth_dating, age$depth_dating, na.rm = T), na.rm = T)), xlab = '', ylab = '')
  mtext(side = 1, line = 2, text = 'Median age [yr BP]')
  mtext(side = 2, line = 2, text = 'Depth from top [mm]')
  
  if(!plyr::empty(data.frame(not.age))){
    points(x = not.age$corr_age, y=not.age$depth_dating, lty = 2, col = 'red', pch = not.age$.pch)
    arrows(not.age$corr_age-not.age$corr_age_uncert, not.age$depth_dating, not.age$corr_age+not.age$corr_age_uncert, not.age$depth_dating, length=0.05, angle=90, code=3, col = 'red')
  }
  
  points(x = age$corr_age, y=age$depth_dating, lty = 2, col = 'black', pch = age$.pch)
  arrows(age$corr_age-age$corr_age_uncert, age$depth_dating, age$corr_age+age$corr_age_uncert, age$depth_dating, length=0.05, angle=90, code=3, col = 'black')
  
  #if(!plyr::empty(data.frame(new.age))){
  #points(x = new.age$corr_age, y=new.age$depth_dating, lty = 2, col = 'cor', pch = not.age$.pch)
  #  arrows(new.age$corr_age-new.age$corr_age_uncert, new.age$depth_dating, new.age$corr_age+new.age$corr_age_uncert, new.age$depth_dating, length=0.05, angle=90, code=3, col = 'cornflowerblue')
  #}
  
  if (any(age$date_type=="Event; actively forming")){
    abline(h=age$depth_dating[which(age$date_type=="Event; actively forming")],lty=3,col="orange")
  }
  if (!plyr::empty(data.frame(chrono$hiatus))) {
    abline(h =  chrono$hiatus$depth_sample, col = 'grey', lty = 3)
  }
  
  #orig_am <- chrono$origAM %>% filter(age_model_type != 'NA') %>% distinct(age_model_type)
  
  datetype <- (age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(date_type))$date_type
  pch.datetype <- (age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(.pch))$.pch
  
  if(plyr::empty(data.frame(not.age))){
    legend("topright",legend = c(paste(datetype, '- used')),
           pch = c(pch.datetype),cex=1,
           col=c(rep('black', length(datetype))),ncol=1,title="Date Type",bty="n")
  }
  
  
  
  if(!plyr::empty(data.frame(not.age))){
    not.datetype <- (not.age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(date_type))$date_type
    pch.not.datetype <- (not.age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(.pch))$.pch
    legend("topright",legend = c(paste(datetype, '- used'), paste(not.datetype, '- not used/unkown')),
           pch = c(pch.datetype, pch.not.datetype),cex=1,
           col=c(rep('black', length(datetype)), rep('red', length(not.datetype))),ncol=1,title="Date Type",bty="n")
  }
  
  mtext(line = 0, text = 'Age-depth model')
  
  return(range(not.age$corr_age, not.age$corr_age +not.age$corr_age_uncert, not.age$corr_age -not.age$corr_age_uncert, 
               age$corr_age, age$corr_age +age$corr_age_uncert, age$corr_age -age$corr_age_uncert, 
               chrono$origAM$interp_age, 
               chrono$copRa$AM$copRa_age, 
               chrono$StalAge$AM$StalAge_age,
               chrono$linInterp$AM$lin_interp_age,
               chrono$linReg$AM$lin_reg_age,
               chrono$Bacon$AM$bacon_age,
               chrono$Bchron$AM$bchron_age,
               chrono$OxCal$AM$OxCal_age,
               na.rm=TRUE))
}

plot_gr <- function(chrono,.entity_id, dating., entity_name, x.lim){
  color = c('royalblue3','skyblue','springgreen2','purple', 'sienna1')
  #color = c('royalblue3','skyblue','springgreen2','purple', 'sienna1')
  
  n <- rep(0,6)
  if(!plyr::empty(data.frame(chrono$copRa$AM$copRa_age))){n[1]<- 1}
  if(!plyr::empty(data.frame(chrono$linInterp$AM$lin_interp_age))){n[2]<- 2}
  if(!plyr::empty(data.frame(chrono$linReg$AM$lin_reg_age))){n[3]<- 3}
  if(!plyr::empty(data.frame(chrono$Bacon$AM$bacon_age))){n[4]<- 4}
  if(!plyr::empty(data.frame(chrono$Bchron$AM$bchron_age))){n[5]<- 5}
  
  n <- n[which(n !=0)]
  
  matplot(x = cbind(chrono$copRa$AM$copRa_age[-1], chrono$linInterp$AM$lin_interp_age[-1], chrono$linReg$AM$lin_reg_age[-1], chrono$Bacon$AM$bacon_age[-1], chrono$Bchron$AM$bchron_age[-1]), 
          #x = cbind( chrono$linInterp$AM$lin_interp_age[-1], chrono$linReg$AM$lin_reg_age[-1], chrono$Bacon$AM$bacon_age[-1], chrono$Bchron$AM$bchron_age[-1]), 
          #y = cbind(chrono$copRa$gr$median.gr, chrono$StalAge$gr$median.gr, chrono$linInterp$gr$median.gr, chrono$linReg$gr$median.gr, chrono$Bacon$gr$median.gr, chrono$Bchron$gr$median.gr),
          y = cbind(chrono$copRa$gr$copRa_gr, chrono$linInterp$gr$lin_interp_gr, chrono$linReg$gr$lin_reg_gr, chrono$Bacon$gr$bacon_gr*10, chrono$Bchron$gr$bchron_gr), 
          #y = cbind(chrono$linInterp$gr$lin_interp_gr, chrono$linReg$gr$lin_reg_gr, chrono$Bacon$gr$bacon_gr, chrono$Bchron$gr$bchron_gr), 
          col = color[n], lty = 1, type = 'l', lwd = 1.5, 
          #xlim = range(chrono$copRa$AM$copRa_age + chrono$copRa$AM$copRa_age_uncert_pos, chrono$copRa$AM$copRa_age - chrono$copRa$AM$copRa_age_uncert_neg,
          #             #chrono$StalAge$AM$StalAge_age + chrono$StalAge$AM$StalAge_age_uncert_pos, chrono$StalAge$AM$StalAge_age - chrono$StalAge$AM$StalAge_age_uncert_neg,
          #             chrono$linInterp$AM$lin_interp_age + chrono$linInterp$AM$lin_interp_age_uncert_pos, chrono$linInterp$AM$lin_interp_age - chrono$linInterp$AM$lin_interp_age_uncert_neg,
          #             chrono$linReg$AM$lin_reg_age + chrono$linReg$AM$lin_reg_age_uncert_pos, chrono$linReg$AM$lin_reg_age - chrono$linReg$AM$lin_reg_age_uncert_neg,
          ##             chrono$Bacon$AM$bacon_age + chrono$Bacon$AM$bacon_age_uncert_pos, chrono$Bacon$AM$bacon_age - chrono$Bacon$AM$bacon_age_uncert_neg,
          #             chrono$Bchron$AM$bchron_age + chrono$Bchron$AM$bchron_age_uncert_pos, chrono$Bchron$AM$bchron_age - chrono$Bchron$AM$bchron_age_uncert_neg,
          #             na.rm=TRUE),
          xlim = x.lim,
          ylim = range(chrono$copRa$gr$copRa_gr, chrono$linInterp$gr$lin_interp_gr, chrono$linReg$gr$lin_reg_gr, chrono$Bacon$gr$bacon_gr, chrono$Bchron$gr$bchron_gr, na.rm = T),
          #ylim = range(chrono$copRa$gr$copRa_gr + chrono$copRa$gr$copRa_gr_uncert_pos, chrono$copRa$gr$copRa_gr - chrono$copRa$gr$copRa_gr_uncert_neg,
          #             #chrono$StalAge$gr$StalAge_gr + chrono$StalAge$gr$StalAge_gr_uncert_pos, chrono$StalAge$gr$StalAge_gr - chrono$StalAge$gr$StalAge_gr_uncert_neg,
          #             chrono$linInterp$gr$lin_interp_gr + chrono$linInterp$gr$lin_interp_gr_uncert_pos, chrono$linInterp$gr$lin_interp_gr - chrono$linInterp$gr$lin_interp_gr_uncert_neg,
          #             chrono$linReg$gr$lin_reg_gr + chrono$linReg$gr$lin_reg_gr_uncert_pos, chrono$linReg$gr$lin_reg_gr - chrono$linReg$gr$lin_reg_gr_uncert_neg,
          #             chrono$Bacon$gr$bacon_gr + chrono$Bacon$gr$bacon_gr_uncert_pos, chrono$Bacon$gr$bacon_gr - chrono$Bacon$gr$bacon_gr_uncert_neg,
          #             chrono$Bchron$gr$bchron_gr + chrono$Bchron$gr$bchron_gr_uncert_pos, chrono$Bchron$gr$bchron_gr - chrono$Bchron$gr$bchron_gr_uncert_neg,
          #             na.rm=TRUE),
          xlab = '', ylab = '')
  mtext(side = 1, line = 2, text = 'Median age [yr BP]')
  mtext(side = 2, line = 2, text = 'Growth rate [mm/yr]')
  
  #legend("topleft",legend=c('copRa median gr','lin. interp. median gr','lin. reg. median gr','Bacon median gr','Bchron median gr'),lty=c(1,1,1,1,1),cex=1,
  #       col = c('royalblue3','skyblue','springgreen2','purple', 'sienna1','',"orange","grey"),ncol=1,title=paste("Information"),bty = "n")
  
  mtext(line = 0, text = 'Median ensemble growth rates')
  
}

plot_gr_iqr <- function(chrono){
  color = c('royalblue3','skyblue','springgreen2','purple', 'sienna1')
  
  n <- rep(0,6)
  if(!plyr::empty(data.frame(chrono$copRa$AM$copRa_age))){n[1]<- 1}
  if(!plyr::empty(data.frame(chrono$linInterp$AM$lin_interp_age))){n[2]<- 2}
  if(!plyr::empty(data.frame(chrono$linReg$AM$lin_reg_age))){n[3]<- 3}
  if(!plyr::empty(data.frame(chrono$Bacon$AM$bacon_age))){n[4]<- 4}
  if(!plyr::empty(data.frame(chrono$Bchron$AM$bchron_age))){n[5]<- 5}
  
  n <- n[which(n !=0)]
  
  matplot(y = cbind(chrono$copRa$gr$copRa_gr_uncert_pos + chrono$copRa$gr$copRa_gr_uncert_neg,
                    chrono$linInterp$gr$lin_interp_gr_uncert_pos + chrono$linInterp$gr$lin_interp_gr_uncert_neg, chrono$linReg$gr$lin_reg_gr_uncert_pos + chrono$linReg$gr$lin_reg_gr_uncert_neg,
                    chrono$Bacon$gr$bacon_gr_uncert_pos*10 + chrono$Bacon$gr$bacon_gr_uncert_neg*10, chrono$Bchron$gr$bchron_gr_uncert_pos + chrono$Bchron$gr$bchron_gr_uncert_neg),
          #y = cbind(chrono$copRa$gr$median.gr, chrono$StalAge$gr$median.gr, chrono$linInterp$gr$median.gr, chrono$linReg$gr$median.gr, chrono$Bacon$gr$median.gr, chrono$Bchron$gr$median.gr),
          x = cbind(chrono$copRa$AM$copRa_age[-1], chrono$linInterp$AM$lin_interp_age[-1], chrono$linReg$AM$lin_reg_age[-1], chrono$Bacon$AM$bacon_age[-1], chrono$Bchron$AM$bchron_age[-1]), 
          col = color[n], lty = 1, type = 'l', lwd = 1.5, 
          ylim = range(chrono$copRa$gr$copRa_gr_uncert_pos + chrono$copRa$gr$copRa_gr_uncert_neg, 
                       chrono$linReg$gr$lin_reg_gr_uncert_pos + chrono$linReg$gr$lin_reg_gr_uncert_neg, chrono$linInterp$gr$lin_interp_gr_uncert_pos + chrono$linInterp$gr$lin_interp_gr_uncert_neg,
                       chrono$Bacon$gr$bacon_gr_uncert_pos*10 + chrono$Bacon$gr$bacon_gr_uncert_neg*10, chrono$Bchron$gr$bchron_gr_uncert_pos + chrono$Bchron$gr$bchron_gr_uncert_neg, na.rm = T),
          xlim = range(chrono$copRa$AM$copRa_age[-1], chrono$linInterp$AM$lin_interp_age[-1], chrono$linReg$AM$lin_reg_age[-1], chrono$Bacon$AM$bacon_age[-1], chrono$Bchron$AM$bchron_age[-1], na.rm = T), 
          #ylim = range(chrono$Bacon$gr$median.gr, na.rm = T),
          xlab = '', ylab = '')
  mtext(side = 1, line = 2, text = 'Median age [yr BP]')
  mtext(side = 2, line = 2, text = 'IQR [mm/yr]')
  
  # legend("topleft",legend=c('copRa median gr','StalAge median gr','lin. interp. median gr','lin. reg. median gr','Bacon median gr','Bchron median gr'),lty=c(1,1,1,1,1,1),cex=1,
  #         col = c('royalblue3','forestgreen','skyblue','springgreen2','purple', 'sienna1','',"orange","grey"),ncol=1,title=paste("Information"),bty = "n")
  
  mtext(line = 0, text = 'Growth rate IQR')
  
}


plot_iqr <- function(chrono,.entity_id, dating., entity_name){
  age <- dating. %>% filter(entity_id == .entity_id) %>% filter( date_used == 'yes') %>%
    mutate(corr_age_uncert = 0.5*(corr_age_uncert_pos + corr_age_uncert_neg)) 
  
  color = c('royalblue3','forestgreen','skyblue','springgreen2','purple', 'sienna1', 'brown2')
  
  n <- rep(0,7)
  if(!plyr::empty(data.frame(chrono$copRa$AM$copRa_age))){n[1]<- 1}
  if(!plyr::empty(data.frame(chrono$StalAge$AM$StalAge_age))){n[2]<- 2}
  if(!plyr::empty(data.frame(chrono$linInterp$AM$lin_interp_age))){n[3]<- 3}
  if(!plyr::empty(data.frame(chrono$linReg$AM$lin_reg_age))){n[4]<- 4}
  if(!plyr::empty(data.frame(chrono$Bacon$AM$bacon_age))){n[5]<- 5}
  if(!plyr::empty(data.frame(chrono$Bchron$AM$bchron_age))){n[6]<- 6}
  if(!plyr::empty(data.frame(chrono$OxCal$AM$OxCal_age))){n[7]<- 7}
  
  n <- n[which(n !=0)]
  
  matplot(x = cbind(chrono$copRa$AM$copRa_age_uncert_pos + chrono$copRa$AM$copRa_age_uncert_neg, chrono$StalAge$AM$StalAge_age_uncert_pos + chrono$StalAge$AM$StalAge_age_uncert_neg, 
                    chrono$linInterp$AM$lin_interp_age_uncert_pos + chrono$linInterp$AM$lin_interp_age_uncert_neg, chrono$linReg$AM$lin_reg_age_uncert_pos + chrono$linReg$AM$lin_reg_age_uncert_neg,
                    chrono$Bacon$AM$bacon_age_uncert_pos + chrono$Bacon$AM$bacon_age_uncert_neg, chrono$Bchron$AM$bchron_age_uncert_pos+chrono$Bchron$AM$bchron_age_uncert_neg,
                    chrono$OxCal$AM$OxCal_age_uncert_pos + chrono$OxCal$AM$OxCal_age_uncert_neg),
          y= chrono$origAM$depth_sample, col = color[n], lty = 1, type = 'l', lwd = 1, 
          xlim = range(chrono$copRa$AM$copRa_age_uncert_pos + chrono$copRa$AM$copRa_age_uncert_neg, chrono$StalAge$AM$StalAge_age_uncert_pos + chrono$StalAge$AM$StalAge_age_uncert_neg, 
                       chrono$linInterp$AM$lin_interp_age_uncert_pos + chrono$linInterp$AM$lin_interp_age_uncert_neg, chrono$linReg$AM$lin_reg_age_uncert_pos + chrono$linReg$AM$lin_reg_age_uncert_neg,
                       chrono$Bacon$AM$bacon_age_uncert_pos + chrono$Bacon$AM$bacon_age_uncert_neg, chrono$Bchron$AM$bchron_age_uncert_pos+chrono$Bchron$AM$bchron_age_uncert_neg,
                       chrono$OxCal$AM$OxCal_age_uncert_pos + chrono$OxCal$AM$OxCal_age_uncert_neg, na.rm=TRUE),
          ylim = c(max(range(chrono$origAM$depth_sample, age$depth_dating, na.rm = T), na.rm = T),min(range(chrono$origAM$depth_sample, age$depth_dating, na.rm = T), na.rm = T)),
          xlab = '', ylab = '', log = 'x')
  mtext(side = 1, line = 2, text = 'IQR [yr]')
  abline(h = age$depth_dating, lty = 3, col = 'black')
  
  if (any(age$date_type=="Event; actively forming")){
    abline(h=age$depth_dating[which(age$date_type=="Event; actively forming")],lty=3,col="orange")
  }
  if (!plyr::empty(data.frame(chrono$hiatus))) {
    abline(h =  chrono$hiatus$depth_sample, col = 'grey', lty = 3)
  }
  
  #legend("right",legend=c('copRa iqr','StalAge iqr','lin. interp. iqr','lin. reg. iqr','Bacon iqr','Bchron iqr','','Dating depth',"Event; actively forming", "Event; hiatus"),lty=c(1,1,1,1,1,1,NA,3,3,3),
  #       cex=0.75,col=c('royalblue3','forestgreen','skyblue','springgreen2','purple', 'sienna1','','black',"orange","grey"),ncol=1,title=paste("Information"),box.col = 'white', bg = 'white')
  
  mtext(line = 0, text = 'Age IQR')
  
}

plot_isotopes <- function(chrono, .entity_id, entity_name, sisal.chrono, x.lim){
  iso<- chrono$proxy
  if(!all(is.na(iso$interp_age))){
   # print('1')
    if(all(is.na(iso$d18O_measurement))){
      matplot(x = iso$interp_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
              #xlim = range(iso$interp_age,na.rm=TRUE),
              xlim = x.lim,
              #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
              ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
      mtext(side = 1, line = 2, text = 'Original AM age [yr BP]')
      axis(side = 4)
      mtext(side = 4, line = 3, text = expression(paste(delta^{13},'C [\u2030]')))
      legend('topleft', legend = c(expression(paste(delta^{13},'C'))), col = c('goldenrod2'),
             lty = 1, title = 'Information',bty = 'n', lwd = 1.5)
      
    } else {
      #print('2')
      matplot(x = iso$interp_age, y= iso$d18O_measurement, col = 'dodgerblue', lty = 1, type = 'l', lwd = 1.5, 
              xlim = x.lim,
              #xlim = range(iso$interp_age,na.rm=TRUE),
              #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
              ylim = range(iso$d18O_measurement, na.rm = T), xlab = '', ylab = '')
      mtext(side = 2, line = 2, text = expression(paste(delta^{18},'O [\u2030]')))
      mtext(side = 1, line = 2, text = 'Original AM age [yr BP]')
      if(!all(is.na(iso$d13C_measurement))){
        par(new = T)
        plot(x = iso$interp_age, y = iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5,
             xlim = x.lim, ylim = range(iso$d13C_measurement, na.rm = T), xlab = NA, ylab = NA, axes = F)
        axis(side = 4)
        mtext(side = 4, line = 3, text = expression(paste(delta^{13},'C [\u2030]')))
        legend('topleft', legend = c(expression(paste(delta^{18},'O')),expression(paste(delta^{13},'C'))), col = c('dodgerblue','goldenrod2'),
               lty = 1, title = 'Information',bty = 'n', lwd = 1.5)
      } else {
        legend('topleft', legend = c(expression(paste(delta^{18},'O'))), col = c('dodgerblue'),
               lty = 1, title = 'Information',bty = 'n', lwd = 1.5)
      }
    }
    
  }else{
    #print('3')
    sc <- sisal.chrono %>% filter(sample_id %in% chrono$proxy$sample_id)
    if(all(is.na(iso$d18O_measurement))){
      
      if(!all(is.na(sc$COPRA_age))){
        if(!plyr::empty(data.frame(chrono$copRa$AM))){
          matplot(x = chrono$copRa$AM$copRa_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
                  xlim = x.lim,
                  #xlim = range(sc$COPRA_age,na.rm=TRUE),
                  #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
                  ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
          mtext(side = 1, line = 2, text = 'copRa median age [yr BP]')
        } else {
          matplot(x = sc$COPRA_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
                  xlim = x.lim,
                  # xlim = range(sc$COPRA_age,na.rm=TRUE),
                  #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
                  ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
          mtext(side = 1, line = 2, text = 'COPRA median age [yr BP]')
        }
        
      } else if(!all(is.na(sc$linear_age))) {
        matplot(x = sc$linear_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
                xlim = x.lim,
                #xlim = range(sc$linear_age,na.rm=TRUE),
                #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
                ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
        mtext(side = 1, line = 2, text = 'Linear median age [yr BP]')
      } else if(!all(is.na(sc$linear_regress_age))) {
        matplot(x = sc$linear_regress_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
                xlim = x.lim,
                #xlim = range(sc$linear_regress_age,na.rm=TRUE),
                #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
                ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
        mtext(side = 1, line = 2, text = 'Linear regression median age [yr BP]')
      } else if(!all(is.na(sc$Bchron_age))) {
        matplot(x = sc$Bchron_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
                xlim = x.lim,
                #xlim = range(sc$Bchron_age,na.rm=TRUE),
                #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
                ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
        mtext(side = 1, line = 2, text = 'Bchron median age [yr BP]')
      } else if(!all(is.na(sc$Bacon_age))) {
        matplot(x = sc$Bacon_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
                xlim = x.lim,
                #xlim = range(sc$Bacon_age,na.rm=TRUE),
                #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
                ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
        mtext(side = 1, line = 2, text = 'Bacon median age [yr BP]')
      }
      axis(side = 4)
      mtext(side = 4, line = 3, text = expression(paste(delta^{13},'C [\u2030]')))
      legend('topleft', legend = c(expression(paste(delta^{13},'C'))), col = c('goldenrod2'),
             lty = 1, title = 'Information',bty = 'n', lwd = 1.5)
    } else {
      if(!all(is.na(sc$COPRA_age))){
        if(!plyr::empty(data.frame(chrono$copRa$AM))){
          matplot(x = chrono$copRa$AM$copRa_age, y= iso$d18O_measurement, col = 'dodgerblue', lty = 1, type = 'l', lwd = 1.5, 
                  xlim = x.lim,
                  #xlim = range(sc$COPRA_age,na.rm=TRUE),
                  #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
                  ylim = range(iso$d18O_measurement, na.rm = T), xlab = '', ylab = '')
          mtext(side = 1, line = 2, text = 'copRa median age [yr BP]')
        } else {
          matplot(x = sc$COPRA_age, y= iso$d18O_measurement, col = 'dodgerblue', lty = 1, type = 'l', lwd = 1.5, 
                  xlim = x.lim,
                  # xlim = range(sc$COPRA_age,na.rm=TRUE),
                  #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
                  ylim = range(iso$d18O_measurement, na.rm = T), xlab = '', ylab = '')
          mtext(side = 1, line = 2, text = 'COPRA median age [yr BP]')
        }
      } else if(!all(is.na(sc$linear_age))) {
        matplot(x = sc$linear_age, y= iso$d18O_measurement, col = 'dodgerblue', lty = 1, type = 'l', lwd = 1.5, 
                xlim = x.lim,
                #xlim = range(sc$linear_age,na.rm=TRUE),
                #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
                ylim = range(iso$d18O_measurement, na.rm = T), xlab = '', ylab = '')
        mtext(side = 1, line = 2, text = 'Linear median age [yr BP]')
      } else if(!all(is.na(sc$linear_regress_age))) {
        matplot(x = sc$linear_regress_age, y= iso$d18O_measurement, col = 'dodgerblue', lty = 1, type = 'l', lwd = 1.5, 
                xlim = x.lim,
                #xlim = range(sc$linear_regress_age,na.rm=TRUE),
                #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
                ylim = range(iso$d18O_measurement, na.rm = T), xlab = '', ylab = '')
        mtext(side = 1, line = 2, text = 'Linear regression median age [yr BP]')
      } else if(!all(is.na(sc$Bchron_age))) {
        matplot(x = sc$Bchron_age, y= iso$d18O_measurement, col = 'dodgerblue', lty = 1, type = 'l', lwd = 1.5, 
                xlim = x.lim,
                #xlim = range(sc$Bchron_age,na.rm=TRUE),
                #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
                ylim = range(iso$d18O_measurement, na.rm = T), xlab = '', ylab = '')
        mtext(side = 1, line = 2, text = 'Bchron median age [yr BP]')
      } else if(!all(is.na(sc$Bacon_age))) {
        matplot(x = sc$Bacon_age, y= iso$d18O_measurement, col = 'dodgerblue', lty = 1, type = 'l', lwd = 1.5, 
                xlim = x.lim,
                #xlim = range(sc$Bacon_age,na.rm=TRUE),
                #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
                ylim = range(iso$d18O_measurement, na.rm = T), xlab = '', ylab = '')
        mtext(side = 1, line = 2, text = 'Bacon median age [yr BP]')
      }
      
      mtext(side = 2, line = 2, text = expression(paste(delta^{18},'O [\u2030]')))
      if(!all(is.na(iso$d13C_measurement))){
        par(new = T)
        if(!all(is.na(sc$COPRA_age))){
          if(!plyr::empty(data.frame(chrono$copRa$AM))){
            plot(x = chrono$copRa$AM$copRa_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
                    xlim = x.lim,
                    #xlim = range(sc$COPRA_age,na.rm=TRUE),
                    #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
                    ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
            mtext(side = 1, line = 2, text = 'copRa median age [yr BP]')
          } else {
            plot(x = sc$COPRA_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
                    xlim = x.lim,
                    # xlim = range(sc$COPRA_age,na.rm=TRUE),
                    #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
                    ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
            mtext(side = 1, line = 2, text = 'COPRA median age [yr BP]')
          }
          
        } else if(!all(is.na(sc$linear_age))) {
          plot(x = sc$linear_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
                  xlim = x.lim,
                  #xlim = range(sc$linear_age,na.rm=TRUE),
                  #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
                  ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
          mtext(side = 1, line = 2, text = 'Linear median age [yr BP]')
        } else if(!all(is.na(sc$linear_regress_age))) {
          plot(x = sc$linear_regress_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
                  xlim = x.lim,
                  #xlim = range(sc$linear_regress_age,na.rm=TRUE),
                  #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
                  ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
          mtext(side = 1, line = 2, text = 'Linear regression median age [yr BP]')
        } else if(!all(is.na(sc$Bchron_age))) {
          plot(x = sc$Bchron_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
                  xlim = x.lim,
                  #xlim = range(sc$Bchron_age,na.rm=TRUE),
                  #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
                  ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
          mtext(side = 1, line = 2, text = 'Bchron median age [yr BP]')
        } else if(!all(is.na(sc$Bacon_age))) {
          plot(x = sc$Bacon_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
                  xlim = x.lim,
                  #xlim = range(sc$Bacon_age,na.rm=TRUE),
                  #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
                  ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
          mtext(side = 1, line = 2, text = 'Bacon median age [yr BP]')
        }
        axis(side = 4)
        mtext(side = 4, line = 3, text = expression(paste(delta^{13},'C [\u2030]')))
        legend('topleft', legend = c(expression(paste(delta^{18},'O')),expression(paste(delta^{13},'C'))), col = c('dodgerblue','goldenrod2'),
               lty = 1, title = 'Information',bty = 'n', lwd = 1.5)
      } else {
        legend('topleft', legend = c(expression(paste(delta^{18},'O'))), col = c('dodgerblue'),
               lty = 1, title = 'Information',bty = 'n', lwd = 1.5)
      }
    }
    
  }
  
  
  mtext(line = 0, text = 'Isotopes')
  
}

plot_dating_detrital <- function(chrono, dating.,.entity_id, entity_name){
  prec <- dating. %>% filter(entity_id == .entity_id) %>% filter(date_used == 'yes' & date_type != 'Event; hiatus') %>% mutate_at(vars(X230Th_232Th_ratio, X230Th_232Th_ratio_uncertainty), as.numeric)
  
  plot(x = prec$dating_id, y= prec$X230Th_232Th_ratio, col = 'black', pch = 1, type = 'p', lwd = 1.5, 
       #xlim = range(iso$interp_age,na.rm=TRUE),
       #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
       ylim = range(prec$X230Th_232Th_ratio + prec$X230Th_232Th_ratio_uncertainty,
                    prec$X230Th_232Th_ratio - prec$X230Th_232Th_ratio_uncertainty, na.rm = T), axes = F, xlab = '', ylab = '')
  arrows(prec$dating_id, prec$X230Th_232Th_ratio-prec$X230Th_232Th_ratio_uncertainty, prec$dating_id, prec$X230Th_232Th_ratio+prec$X230Th_232Th_ratio_uncertainty, length=0.05, angle=90, code=3, col = 'black')
  box(lty = 1)
  axis(side = 1)
  axis(side = 4)
  mtext(side = 4, line = 3, text = expression(paste('(230Th/232Th)')))
  mtext(side = 1, line = 2, text = 'Dating ID')
  
  legend('topleft', legend = 'Dating point', col = c('black'),
         pch = 1, title = 'Information',bty = 'n')
  
  mtext(line = 0, text = 'Detrital 232Th')
}

plot_dating_abs <- function(chrono, dating.,.entity_id, entity_name){
  prec <- dating. %>% filter(entity_id == .entity_id) %>% filter(date_used == 'yes' & date_type != 'Event; hiatus') %>% mutate_at(vars(X230Th_content, X230Th_uncertainty, X238U_content, X232Th_uncertainty), as.numeric)
  
  plot(x = prec$dating_id, y= prec$X230Th_content, col = 'black', pch = 1, type = 'p', lwd = 1.5, 
       ylim = range(prec$X230Th_content + prec$X230Th_uncertainty, prec$X230Th_content - prec$X230Th_uncertainty, na.rm = T), 
       axes = F, xlab = '', ylab = '')
  
  par(new = T)
  plot(x = prec$dating_id, y = prec$X238U_content, col = 'black', pch = 6 , type = 'p', lwd = 1.5,
       ylim = range(prec$X238U_content + prec$X238U_uncertainty, prec$X238U_content - prec$X238U_uncertainty, na.rm = T),
       xlab = NA, ylab = NA, axes = F)
  
  arrows(prec$dating_id, prec$X230Th_content-prec$X230Th_uncertainty, prec$dating_id, prec$X230Th_content+prec$X230Th_uncertainty, length=0.05, angle=90, code=3, col = 'black')
  arrows(prec$dating_id, prec$X238U_content-prec$X238U_uncertainty, prec$dating_id, prec$X238U_content+prec$X238U_uncertainty, length=0.05, angle=90, code=3, col = 'black')
  
  box(lty = 1)
  axis(side = 1)
  axis(side = 4)
  axis(side = 2)
  mtext(side = 4, line = 3, text = 'abs(238U)')
  mtext(side = 1, line = 2, text = 'Dating ID')
  mtext(side = 2, line = 2, text = 'abs(230Th)')
  
  legend('topleft', legend = c('230Th','238U'), col = c('black'),
         pch = c(1,6), title = 'Information',bty = 'n')
  
  mtext(line = 0, text = 'Absolute isotope content')
}


plot_dating <- function(chrono,.entity_id, dating.){
  age <- dating. %>% filter(entity_id == .entity_id) %>% filter( date_used == 'yes' & date_type != 'Event; hiatus' & date_type != 'Event; actively forming') %>%
    mutate_at(vars(corr_age, corr_age_uncert_pos, corr_age_uncert_neg,uncorr_age, uncorr_age_uncert_pos, uncorr_age_uncert_neg), as.numeric) %>% 
    mutate(.pch = if_else(date_type == "MC-ICP-MS U/Th", 2,
                          if_else(date_type == "TIMS", 8,
                                  if_else(date_type == "ICP-MS U/Th Other", 4,
                                          if_else(date_type == "Alpha U/Th", 5,
                                                  if_else(date_type ==  "U/Th unspecified",7, NA_real_))))))
  
  if(!all(is.na(age$uncorr_age))){
    #print('a')
    if(!all(is.na(age$uncorr_age_uncert_pos))){
      plot(x = age$uncorr_age, y= age$corr_age, col = 'black', pch = age$.pch, type = 'p', 
           xlim = range(age$uncorr_age + age$uncorr_age_uncert_pos, age$uncorr_age - age$uncorr_age_uncert_neg,
                        na.rm=TRUE),
           ylim = range(age$corr_age + age$corr_age_uncert_pos,age$corr_age - age$corr_age_uncert_neg,
                        na.rm = T), xlab = '', ylab = '')
      mtext(side = 1, line = 2, text = 'Uncorrected age [yr BP]')
      mtext(side = 2, line = 2, text = 'Corrected age [yr BP]')
      arrows(age$uncorr_age,age$corr_age-age$corr_age_uncert_neg, age$uncorr_age,age$corr_age+age$corr_age_uncert_pos, length=0.05, angle=90, code=3, col = 'black')
      arrows(age$uncorr_age + age$uncorr_age_uncert_pos,age$corr_age, age$uncorr_age - age$uncorr_age_uncert_neg,age$corr_age, length=0.05, angle=90, code=3, col = 'black')
      
      datetype <- (age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(date_type))$date_type
      pch.datetype <- (age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(.pch))$.pch
      
      abline(a = 0, b = 1, lty = 2)
      
      legend("topleft",legend = c(paste(datetype,'age'), 'Angle bissectrice'),
             pch = c(pch.datetype, NA),cex=1, lty = c(rep(NA,length(datetype)), 2),
             col=c(rep('black', length(datetype))),ncol=1,title="Information",bty="n")
    } else {
      plot(x = age$uncorr_age, y= age$corr_age, col = 'black', pch = age$.pch, type = 'p', 
           xlim = range(age$uncorr_age,
                        na.rm=TRUE),
           ylim = range(age$corr_age + age$corr_age_uncert_pos,age$corr_age - age$corr_age_uncert_neg,
                        na.rm = T), xlab = '', ylab = '')
      mtext(side = 1, line = 2, text = 'Uncorrected age [yr BP]')
      mtext(side = 2, line = 2, text = 'Corrected age [yr BP]')
      arrows(age$uncorr_age,age$corr_age-age$corr_age_uncert_neg, age$uncorr_age,age$corr_age+age$corr_age_uncert_pos, length=0.05, angle=90, code=3, col = 'black')
      
      datetype <- (age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(date_type))$date_type
      pch.datetype <- (age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(.pch))$.pch
      
      abline(a = 0, b = 1, lty = 2)
      
      legend("topleft",legend = c(paste(datetype,'age'), 'Angle bissectrice'),
             pch = c(pch.datetype, NA),cex=1, lty = c(rep(NA,length(datetype)), 2),
             col=c(rep('black', length(datetype))),ncol=1,title="Information",bty="n")
    }
    
    mtext(line = 0, text = 'Evaluation of the detrital contribution')
  } else {
    plot.new()
    legend("center", legend = 'No uncorrected age', col = 'black', bty = 'n')
  }
}

final.plot<- function(exp, .entity ,.dating, drop){
  eId <- (exp %>% distinct(entity_id) %>% arrange(entity_id))$entity_id
  sisal.chrono <- read.csv('/home/ariana/SISAL Data/sisalv2/sisal_chronology.csv', header = T, stringsAsFactors = F) %>% mutate_at(vars(everything()), as.numeric)
  #cairo_pdf(paste("~/SISAL/v2/",i,'-',entity_name,'.pdf',sep = ''), 12, 12)
  
  for(i in eId) {
    
    print(i)
    #m <- runFile %>% filter(entity_id == i) 
    DROP <- (drop %>% filter(entity_id == i))$DROP
    
    #if(any(c(m$Bacon,m$Bchron, m$StalAge, m$linInterp, m$copRa, m$linReg, m$OxCal)) & is.na(DROP)){
    entity_name <- (.entity %>% filter(entity_id == i))$entity_name
    graphics.off()
    df_fil <- get(load(paste('/home/ariana/SISAL/v2_upd_2/',i,'-',entity_name, '.RData' ,sep = '')))
      
    cairo_pdf(paste("~/SISAL/v2_upd_2/",i,'-',entity_name,'.pdf',sep = ''), 12, 12)
    layout(matrix(c(1,1,2,3,3,4,5,5,6,7,7,7), nrow=4, ncol=3, byrow = T))
    par(mar = c(4,3,3,4), oma = c(1,3,2,1))
    x.lim <- plot_am(df_fil, i, .dating, entity_name)
    if(any(c(!is_empty(data.frame(df_fil$Bacon$AM)),!is_empty(data.frame(df_fil$Bchron$AM)),!is_empty(data.frame(df_fil$copRa$AM)),!is_empty(data.frame(df_fil$StalAge$AM)),!is_empty(data.frame(df_fil$linInterp$AM)),
             !is_empty(data.frame(df_fil$linReg$AM)),!is_empty(data.frame(df_fil$OxCal$AM)))) & is.na(DROP)){
    #if(any(c(m$Bacon,m$Bchron, m$StalAge, m$linInterp, m$copRa, m$linReg, m$OxCal)) & is.na(DROP)){
      plot_iqr(df_fil, i, .dating, entity_name)
      mtext(paste(i,'-', entity_name, sep = ''), side = 3, line = -2, outer = TRUE)
      if(all(c(is_empty(data.frame(df_fil$Bacon$gr)),is_empty(data.frame(df_fil$Bchron$gr)),is_empty(data.frame(df_fil$copRa$gr)),is_empty(data.frame(df_fil$linInterp$gr)),
               is_empty(data.frame(df_fil$linReg$gr))))){
        plot.new()
        plot.new()
      } else {
        plot_gr(df_fil,i, .dating, entity_name, x.lim)
        plot_gr_iqr(df_fil)
      }
    } else {
      plot.new()
      mtext(paste(i,'-', entity_name, sep = ''), side = 3, line = -2, outer = TRUE)
      plot.new()
      plot.new()
    }
    plot_isotopes(df_fil, i, entity_name, sisal.chrono, x.lim)
    plot_dating(df_fil, i, .dating)
    plot.new()
      
    orig_am <- df_fil$origAM %>% filter(age_model_type != 'NA') %>% distinct(age_model_type)
      
    legend('center',legend=c(paste('original AM:', orig_am),'copRa median age','StalAge mean age','lin. interp. median age','lin. reg. median age','Bacon median age','Bchron median age', 'OxCal median age', "Event; actively forming", "Event; hiatus"),
           lty=c(1,1,1,1,1,1,1,1,3,3),cex=1.5,
           col = c('black','royalblue3','forestgreen','skyblue','springgreen2','purple', 'sienna1', 'brown2',"orange","grey"),title=paste("Age model Information"),bty = "n",  ncol = 3, lwd = 4)
      
    dev.off()
    
    
  }
  
}


proj.sites <- function(site){
  site_pts <- data.frame(cbind(site$longitude, site$latitude))
  names(site_pts) <- c('lon', 'lat')
  
  coordinates(site_pts) <- ~lon+lat
  proj4string(site_pts) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  summary(site_pts)
  
  robin.crs <- CRS("+proj=robin +lon_0=0w")
  site_pts.proj <- spTransform(site_pts, robin.crs)
  
  return(site_pts.proj)
}

plot_world_map <- function(totsite_pts.proj, site){
  # read the projected Robinson shapefiles
  shapepath <- "~/Data/Welt Karte Robinson/derived/glrob_50m/"
  coast_lines_proj <- readShapeLines(paste(shapepath, "glRob_50m_coast_lines.shp", sep=""))
  bb_lines_proj <- readShapeLines(paste(shapepath, "glRob_50m_bb_lines.shp", sep=""))
  bb_poly_proj <- readShapePoly(paste(shapepath, "glRob_50m_bb_poly.shp", sep=""))
  grat30_lines_proj <- readShapeLines(paste(shapepath, "glRob_50m_grat30_lines.shp", sep=""))
  land_proj <- readShapePoly(paste(shapepath,'glRob_50m_land.shp', sep =''))
  
  plot(bb_poly_proj, col="aliceblue", bor="black", lwd=0.1)
  plot(land_proj, col="lightgrey", bor="gray50", lwd=0.4, add = TRUE)
  plot(grat30_lines_proj, col="black", lwd=0.3, add=TRUE)
  plot(coast_lines_proj, col="black", lwd=0.5, add=TRUE)
  plot(bb_lines_proj, col="black", lwd=1.0, add=TRUE)
  
  site_pts.proj <- proj.sites(site)
  plot(site_pts.proj, pch=0, col='indianred', cex=2, lwd=1.6, add=TRUE)
  
}

plot.oxcal<- function(oxcal,.entity ,.dating){
  eId <- (oxcal %>% distinct(entity_id) %>% arrange(entity_id))$entity_id
  for(i in eId) {
    
    print(i)
    m <- oxcal %>% filter(entity_id == i) 
    
    entity_name <- (.entity %>% filter(entity_id == i))$entity_name
    graphics.off()
    pdf(paste("~/SISAL/Oxcal/",i,'-',entity_name,'.pdf',sep = ''), 6,4)
    plot_am_oxcal(m, i, .dating, entity_name)
    dev.off()
    
  }
  
}

plot_am_oxcal <- function(chrono,.entity_id, dating., entity_name){
  not.age <- dating. %>% filter(entity_id == .entity_id) %>% filter(date_used_OxCal == 'no'| date_used_OxCal == 'unknown') %>%
    mutate(corr_age_uncert = 0.5*(corr_age_uncert_pos + corr_age_uncert_neg)) %>%
    mutate(.pch = if_else(date_type == "MC-ICP-MS U/Th", 2,
                          if_else(date_type == "TIMS", 8,
                                  if_else(date_type == "ICP-MS U/Th Other", 4,
                                          if_else(date_type == "Alpha U/Th", 5,
                                                  if_else(date_type ==  "U/Th unspecified",7, 
                                                          if_else(date_type == 'Event; start of laminations' |
                                                                    date_type == 'Event; end of laminations', 3, 
                                                                  if_else(date_type == 'C14', 12, NA_real_)))))))) 
  
  age <- dating. %>% filter(entity_id == .entity_id) %>% filter( date_used_OxCal == 'yes') %>%
    mutate(corr_age_uncert = 0.5*(corr_age_uncert_pos + corr_age_uncert_neg)) %>% 
    mutate(.pch = if_else(date_type == "MC-ICP-MS U/Th", 2,
                          if_else(date_type == "TIMS", 8,
                                  if_else(date_type == "ICP-MS U/Th Other", 4,
                                          if_else(date_type == "Alpha U/Th", 5,
                                                  if_else(date_type ==  "U/Th unspecified",7, 
                                                          if_else(date_type ==  "Event; end of laminations",10,
                                                                  if_else(date_type ==  "Event; start of laminations",12, 11))))))))
  
  #new.age <- read.csv(file.path(file_path, paste(.entity_id,'-',entity_name, sep = ''),'/copRa/new_dating_tb.csv'), header = T) %>% 
  #  filter(corr_age_uncert != age$corr_age_uncert) %>%
  #  mutate(.pch = if_else(date_type == "MC-ICP-MS U/Th", 2,
  #                        if_else(date_type == "TIMS", 8,
  #                                if_else(date_type == "ICP-MS U/Th Other", 4,
  #                                        if_else(date_type == "Alpha U/Th", 5,
  #                                                if_else(date_type ==  "U/Th unspecified",7, NA_real_))))))
  
  matplot(x = chrono$OxCal_age, y= chrono$depth, col = 'deepskyblue4', lty = 1, type = 'l', lwd = 1, 
          xlim = range(not.age$corr_age, not.age$corr_age +not.age$corr_age_uncert, not.age$corr_age -not.age$corr_age_uncert, 
                       age$corr_age, age$corr_age +age$corr_age_uncert, age$corr_age -age$corr_age_uncert, 
                       chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
          #xlim = range(chrono$age,chrono$age+chrono$uncert_pos+5,chrono$age-chrono$uncert_neg-5,na.rm=TRUE),
          ylim = c(max(range(chrono$depth, not.age$depth_dating, age$depth_dating, na.rm = T), na.rm = T),min(range(chrono$depth, not.age$depth_dating, age$depth_dating, na.rm = T), na.rm = T)), 
          xlab = 'OxCal age [yrs BP]', ylab = 'Depth from top [mm]')
  lines(x=chrono$OxCal_age + chrono$OxCal_age_uncert_pos,y=chrono$depth, lty = 2, col = 'red')
  lines(x=chrono$OxCal_age - chrono$OxCal_age_uncert_neg,y = chrono$depth, lty = 2, col = 'red')
  
  if(!plyr::empty(data.frame(not.age))){
    points(x = not.age$corr_age, y=not.age$depth_dating, lty = 2, col = 'red', pch = not.age$.pch)
    arrows(not.age$corr_age-not.age$corr_age_uncert, not.age$depth_dating, not.age$corr_age+not.age$corr_age_uncert, not.age$depth_dating, length=0.05, angle=90, code=3, col = 'red')
  }
  
  points(x = age$corr_age, y=age$depth_dating, lty = 2, col = 'black', pch = age$.pch)
  arrows(age$corr_age-age$corr_age_uncert, age$depth_dating, age$corr_age+age$corr_age_uncert, age$depth_dating, length=0.05, angle=90, code=3, col = 'black')
  
  #if(!plyr::empty(data.frame(new.age))){
  #points(x = new.age$corr_age, y=new.age$depth_dating, lty = 2, col = 'cor', pch = not.age$.pch)
  #  arrows(new.age$corr_age-new.age$corr_age_uncert, new.age$depth_dating, new.age$corr_age+new.age$corr_age_uncert, new.age$depth_dating, length=0.05, angle=90, code=3, col = 'cornflowerblue')
  #}
  
  if (any(age$date_type=="Event; actively forming")){
    abline(h=age$depth_dating[which(age$date_type=="Event; actively forming")],lty=3,col="orange")
  }
  if (!plyr::empty(data.frame(chrono$hiatus))) {
    abline(h =  chrono$hiatus$depth_sample, col = 'grey', lty = 3)
  }
  
  datetype <- (age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(date_type))$date_type
  pch.datetype <- (age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(.pch))$.pch
  
  if(plyr::empty(data.frame(not.age))){
    legend("topright",legend = c(paste(datetype, '- used')),
           pch = c(pch.datetype),cex=1,
           col=c(rep('black', length(datetype))),ncol=1,title="Date Type",bty="n")
  }
  
  
  
  if(!plyr::empty(data.frame(not.age))){
    not.datetype <- (not.age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(date_type))$date_type
    pch.not.datetype <- (not.age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(.pch))$.pch
    legend("topright",legend = c(paste(datetype, '- used'), paste(not.datetype, '- not used/unkown')),
           pch = c(pch.datetype, pch.not.datetype),cex=1,
           col=c(rep('black', length(datetype)), rep('red', length(not.datetype))),ncol=1,title="Date Type",bty="n")
  }
  
  mtext(line = 0, text = 'Age-depth model')
  
}






