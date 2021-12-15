##### quantitative analysis hiatus = 0#####

setwd("/home/ariana/R/SISAL_AM")
source(file = 'run_functions.R')
source(file = 'SISAL_AM_add_functions.R')
source(file = 'SISAL_AM_plot_functions.R')
source(file = "StalAge_1_0_mine.R")
source(file = 'New_methods.R')
source(file = 'lin_interp_reversal.R')
source(file = '~/R/syntheticspeleo-master/R/functions.R')
#source(file = '~/R/syntheticspeleo-master/R/validation_fct.R')

library(zoo)
library(Hmisc)
library(pracma)
library(plyr)
library(tidyverse)


# d.min <- 0
# d.max <- 1000
# t.min <- 0
# t.max <- 10000
# t.n <- 4
# h.length <- c(1000)
# rate <- 3
# t.res <- 2

# n.samp <- 15
# p.min <- 0.1
# p.max <- 1
# p.n <- 5
# n <- 10
# file_path <- '~/synSpeleo/test/'
# nsamp <- 20
#h.length=c(2000,3000,500)

quanAna_AM <- function(precision,t.min=0,t.max=10000,t.n=4,d.min=0,d.max=1000,h.length,file_path1, file_path2 ,rate=3,t.res=2, nsamp=15,n, k){

  set.seed(2020)
  chrono_final <- tibble()


  if(k==0){
    t <- t.max-t.min
    l <- 0
  } else {
    h.n <- h.length[1:k]
    t <- t.max-t.min - sum(h.n)
    l <- length(h.n)
  }


  prec <- precision
  #prec <- round(cumsum(c(0.1,rep((p.max-p.min)/p.n,p.n))),digits=1)
  #if(k ==0){h.n <- tibble(sample_id = numeric(), depth_sample = numeric())}

  # go through precisions
  for (j in prec){
    print(j)

    for (i in 1:n) {
      print(i)

      # create speleo object
      speleo<-NA
      speleo_new <- NA
      v <- NA
      speleo<-generate.speleo(t.top=t.min,t.bottom=t.max,z.top = d.min,z.bottom=d.max,t.res=t.res)


      if(k==0){
        v <- sample(t.n)
        tlist2 <- c(t.min,c(cumsum(c(rep(t/t.n,t.n))[v])))
        rates2 <- rep(rate,t.n)[v]
      } else {
        v <- sample(seq(2,length(h.n)+t.n-1,1))
        tlist2 <- c(t.min,c(cumsum(c(t/t.n,h.n,rep(t/t.n,t.n-1))[c(1,v,length(h.n)+t.n)])))
        rates2 <- c(rate,rep(0,l),rep(rate,t.n-1))[c(1,v,length(h.n)+t.n)]
      }



      print(tlist2)
      print(rates2)

      speleo_new<-growthfcn.changinggamma.new(speleo, tlist = tlist2, rates = rates2)

      speleo_new$datingparams$tdofs<-10
      speleo_new$datingparams$precision<-j
      speleo_new$datingparams$topage<-10
      speleo_new$datingparams$nsample<-nsamp

      speleo_new<-do.dating(speleo_new)

      plotgrowthfcn(speleo_new)
      plotdating(speleo_new$Dtable,add=TRUE)

      pdf(file.path(file_path2,paste('grwthfct-',l,'-',j,'-',i,'.pdf',sep = '')),6,4)
      plotgrowthfcn(speleo_new)
      plotdating(speleo_new$Dtable,add=TRUE)
      graphics.off()

      hiatus <- unique(speleo_new$growthfcn$expected.z[which(diff(speleo_new$growthfcn$expected.z)==0)])


      if(k==0){hiatus <-tibble(sample_id = numeric(), depth_sample = numeric())} else {hiatus <- tibble(sample_id = -i*10-seq(1:length(hiatus)), depth_sample = hiatus)}

      speleo_new<-sample.proxy(speleo = speleo_new,nsamp = 700,z.steps = 1,stype = "sinusoid")
      speleo_new$hiatus <- hiatus
      #dind<-check.Dtable(speleo_class1,add=TRUE)
      save(speleo_new,file = file.path(file_path2, paste('grwthfct-',l,'-',j,'-',i,'.RData',sep='')))
      run_SISAL_chrono(speleo_new, file_path1,paste(l,'-',j,'-',i,sep = '') , hiatus)

      file_name <- paste(l,'-',j,'-',i,sep='')

      if(k==0){
        sample <- speleo_new$proxy$proxymat %>% mutate(entity_id = rep(file_name,length(speleo_new$proxy$proxymat$depth_sample))) %>% select(sample_id, depth_sample, entity_id) %>%
          mutate(sample_id = paste(sample_id,'-',file_name,sep='')) %>%
          arrange(.,depth_sample) %>% mutate_at(vars(entity_id),as.character)

      } else {
        sample <- bind_rows(speleo_new$proxy$proxymat %>% mutate(entity_id = rep(file_name,length(speleo_new$proxy$proxymat$depth_sample))) %>% select(sample_id, depth_sample, entity_id) %>%
                              mutate(sample_id = paste(sample_id,'-',file_name,sep='')),
                            tibble(sample_id = paste(hiatus$sample_id,'-',file_name,sep=''), depth_sample = hiatus$depth_sample, entity_id =rep(file_name,length(hiatus$depth_sample)))) %>%
          arrange(.,depth_sample) %>% mutate_at(vars(entity_id),as.character)

        dt <- speleo_new$Dtable %>% select(dating_id, sampling.depths, sampling.age.est) %>% mutate(h = NA) %>%
          bind_rows(., tibble(dating_id = seq(1:length(hiatus$depth_sample)), sampling.depths = hiatus$depth_sample, h = rep('h',length(hiatus$depth_sample)), sampling.age.est = rep(NA,length(hiatus$depth_sample)))) %>%
          arrange(.,sampling.depths) %>%
          mutate(help = if_else(lead(h)=='h' | lag(h)=='h',1,0)) %>% filter(help == 1)
        dt.age <- (dt %>% mutate(diff = lead(sampling.age.est)-sampling.age.est))[seq(1,dim(dt)[1],2),]
        dt.delete <- dt.age %>% filter(is.na(diff)|diff>2000)
        d.age.na <- NA
        d.age.2000 <- NA
        if(any(is.na(dt.delete$diff))){d.age.na <- (dt.delete %>% filter(is.na(diff)))$sampling.depth} ### hiatus in the end
        if(any(dt.delete$diff>2000)&!is.na(any(dt.delete$diff>2000))){
          print('i')
          d.fil <- (dt.delete %>% filter(diff > 2000))$sampling.depth
          d.age.2000 <- (dt[c(which(dt$sampling.depths %in% d.fil),which(dt$sampling.depths %in% d.fil)+1),])$sampling.depths
        }
      }

      chrono <- merge_synSpeleo_qa_AM(file_name, file_path1, sample) %>% arrange(., depth_sample)

      ind.samples<-sapply(chrono$depth_sample,function(x,growthfcn){which.min(abs(growthfcn$growth.real-x))},growthfcn=speleo_new$growthfcn)
      true.age<-speleo_new$growthfcn$true.age[ind.samples]
      #true.age <- grf$true_age[sapply(chrono$depth_sample,function(x) which.min(abs(x-grf$expected.z)))]

      speleo_new$recorded.age <- true.age



      chrono_new <- chrono %>% mutate(true.age = true.age, precision = rep(j, length(chrono$depth_sample)),
                                      hiatus_nr = rep(k, length(chrono$depth_sample)), iteration = rep(i, length(chrono$depth_sample))) %>%
        mutate(lin_reg_age = if_else(lin_reg_age < -68, NA_real_,lin_reg_age),
               lin_reg_age_uncert_pos = if_else(lin_reg_age < -68, NA_real_,lin_reg_age_uncert_pos),
               lin_reg_age_uncert_neg = if_else(lin_reg_age < -68, NA_real_,lin_reg_age_uncert_neg),
               lin_interp_age = if_else(lin_interp_age < -68, NA_real_,lin_interp_age),
               lin_interp_age_uncert_pos = if_else(lin_interp_age < -68, NA_real_,lin_interp_age_uncert_pos),
               lin_interp_age_uncert_neg = if_else(lin_interp_age < -68, NA_real_,lin_interp_age_uncert_neg),
               copRa_age = if_else(copRa_age < -68, NA_real_,copRa_age),
               copRa_age_uncert_pos = if_else(copRa_age < -68, NA_real_,copRa_age_uncert_pos),
               copRa_age_uncert_neg = if_else(copRa_age < -68, NA_real_,copRa_age_uncert_neg),
               StalAge_age = if_else(StalAge_age < -68, NA_real_,StalAge_age),
               StalAge_age_uncert_pos = if_else(StalAge_age < -68, NA_real_,StalAge_age_uncert_pos),
               StalAge_age_uncert_neg = if_else(StalAge_age < -68, NA_real_,StalAge_age_uncert_neg),
               bacon_age = if_else(bacon_age < -68, NA_real_,bacon_age),
               bacon_age_uncert_pos = if_else(bacon_age < -68, NA_real_,bacon_age_uncert_pos),
               bacon_age_uncert_neg = if_else(bacon_age < -68, NA_real_,bacon_age_uncert_neg))
      if(k!=0){
        if(!is.na(d.age.na)){
          chrono_new <- chrono_new %>% mutate(lin_reg_age = if_else(depth_sample > d.age.na, NA_real_,lin_reg_age),
                                              lin_reg_age_uncert_pos = if_else(depth_sample > d.age.na, NA_real_,lin_reg_age_uncert_pos),
                                              lin_reg_age_uncert_neg = if_else(depth_sample > d.age.na, NA_real_,lin_reg_age_uncert_neg),
                                              lin_interp_age = if_else(depth_sample > d.age.na, NA_real_,lin_interp_age),
                                              lin_interp_age_uncert_pos = if_else(depth_sample > d.age.na, NA_real_,lin_interp_age_uncert_pos),
                                              lin_interp_age_uncert_neg = if_else(depth_sample > d.age.na, NA_real_,lin_interp_age_uncert_neg),
                                              copRa_age = if_else(depth_sample > d.age.na, NA_real_,copRa_age),
                                              copRa_age_uncert_pos = if_else(depth_sample > d.age.na, NA_real_,copRa_age_uncert_pos),
                                              copRa_age_uncert_neg = if_else(depth_sample > d.age.na, NA_real_,copRa_age_uncert_neg),
                                              StalAge_age = if_else(depth_sample > d.age.na, NA_real_,StalAge_age),
                                              StalAge_age_uncert_pos = if_else(depth_sample > d.age.na, NA_real_,StalAge_age_uncert_pos),
                                              StalAge_age_uncert_neg = if_else(depth_sample > d.age.na, NA_real_,StalAge_age_uncert_neg),
                                              bacon_age = if_else(depth_sample > d.age.na, NA_real_,bacon_age),
                                              bacon_age_uncert_pos = if_else(depth_sample > d.age.na, NA_real_,bacon_age_uncert_pos),
                                              bacon_age_uncert_neg = if_else(depth_sample > d.age.na, NA_real_,bacon_age_uncert_neg))
        }

        if(!is.na(d.age.2000)){
          chrono_new <- chrono_new %>% mutate(lin_reg_age = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,lin_reg_age),
                                              lin_reg_age_uncert_pos = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,lin_reg_age_uncert_pos),
                                              lin_reg_age_uncert_neg = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,lin_reg_age_uncert_neg),
                                              lin_interp_age = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,lin_interp_age),
                                              lin_interp_age_uncert_pos = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,lin_interp_age_uncert_pos),
                                              lin_interp_age_uncert_neg = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,lin_interp_age_uncert_neg),
                                              copRa_age = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,copRa_age),
                                              copRa_age_uncert_pos = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,copRa_age_uncert_pos),
                                              copRa_age_uncert_neg = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,copRa_age_uncert_neg),
                                              StalAge_age = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,StalAge_age),
                                              StalAge_age_uncert_pos = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,StalAge_age_uncert_pos),
                                              StalAge_age_uncert_neg = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,StalAge_age_uncert_neg),
                                              bacon_age = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,bacon_age),
                                              bacon_age_uncert_pos = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,bacon_age_uncert_pos),
                                              bacon_age_uncert_neg = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,bacon_age_uncert_neg))

          if(length(d.age.2000) > 2) {
            chrono_new <- chrono_new %>% mutate(lin_reg_age = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,lin_reg_age),
                                                lin_reg_age_uncert_pos = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,lin_reg_age_uncert_pos),
                                                lin_reg_age_uncert_neg = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,lin_reg_age_uncert_neg),
                                                lin_interp_age = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,lin_interp_age),
                                                lin_interp_age_uncert_pos = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,lin_interp_age_uncert_pos),
                                                lin_interp_age_uncert_neg = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,lin_interp_age_uncert_neg),
                                                copRa_age = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,copRa_age),
                                                copRa_age_uncert_pos = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,copRa_age_uncert_pos),
                                                copRa_age_uncert_neg = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,copRa_age_uncert_neg),
                                                StalAge_age = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,StalAge_age),
                                                StalAge_age_uncert_pos = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,StalAge_age_uncert_pos),
                                                StalAge_age_uncert_neg = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,StalAge_age_uncert_neg),
                                                bacon_age = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,bacon_age),
                                                bacon_age_uncert_pos = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,bacon_age_uncert_pos),
                                                bacon_age_uncert_neg = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,bacon_age_uncert_neg))
          }
        }
      }
      chrono_final_0 <- bind_rows(chrono_final,chrono_new)
      save(chrono_new, file = file.path(file_path2, paste('chrono-',l,'-',j,'-',i,'_AM.RData',sep='')))

    }
  }
  return(chrono_final)
}


merge_synSpeleo_qa_AM <- function(file_name, working_directory, se){

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
                                     bacon_age_uncert_neg = rep(NA_real_, length(se$sample_id)))%>%
    mutate_at(vars(entity_id,sample_id),as.character) %>% arrange(., depth_sample)

  setwd(file.path(working_directory, file_name, '/linReg'))
  if(file.exists('linReg_chronology.csv')){
    lR_chrono <- read.csv('linReg_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>%
      distinct(sample_id, .keep_all = T) %>% mutate(sample_id = paste(sample_id,'-',file_name,sep=''))
    names(lR_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
    SISAL_chronology_new <- full_join(lR_chrono, SISAL_chronology_new, by = 'sample_id') %>%
      mutate(lin_reg_age = if_else(sample_id %in% lR_chrono$sample_id, age, lin_reg_age),
             lin_reg_age_uncert_pos = if_else(sample_id %in% lR_chrono$sample_id, uncert_pos, lin_reg_age_uncert_pos),
             lin_reg_age_uncert_neg = if_else(sample_id %in% lR_chrono$sample_id, uncert_neg, lin_reg_age_uncert_neg)) %>%
      select(-age, -uncert_pos, -uncert_neg)
  }

  setwd(file.path(working_directory, file_name, '/linInterp'))
  if(file.exists('linInt_chronology.csv')){
    lI_chrono <- read.csv('linInt_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>%
      distinct(sample_id, .keep_all = T) %>% mutate(sample_id = paste(sample_id,'-',file_name,sep=''))
    names(lI_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
    SISAL_chronology_new <- full_join(lI_chrono, SISAL_chronology_new, by = 'sample_id')  %>%
      mutate(lin_interp_age = if_else(sample_id %in% lI_chrono$sample_id, age, lin_interp_age),
             lin_interp_age_uncert_pos = if_else(sample_id %in% lI_chrono$sample_id, uncert_pos, lin_interp_age_uncert_pos),
             lin_interp_age_uncert_neg = if_else(sample_id %in% lI_chrono$sample_id, uncert_neg, lin_interp_age_uncert_neg)) %>%
      select(-age, -uncert_pos, -uncert_neg)
  }

  setwd(file.path(working_directory, file_name, '/copRa'))
  if(file.exists('copRa_chronology.csv')){
    copRa_chrono <- read.csv('copRa_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>%
      distinct(sample_id, .keep_all = T) %>% mutate(sample_id = paste(sample_id,'-',file_name,sep=''))
    names(copRa_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
    SISAL_chronology_new <- full_join(copRa_chrono, SISAL_chronology_new, by = 'sample_id')  %>%
      mutate(copRa_age = if_else(sample_id %in% copRa_chrono$sample_id, age, copRa_age),
             copRa_age_uncert_pos = if_else(sample_id %in% copRa_chrono$sample_id, uncert_pos, copRa_age_uncert_pos),
             copRa_age_uncert_neg = if_else(sample_id %in% copRa_chrono$sample_id, uncert_neg, copRa_age_uncert_neg)) %>%
      select(-age, -uncert_pos, -uncert_neg)
  }

  setwd(file.path(working_directory, file_name, '/StalAge'))
  if(file.exists('StalAge_chronology.csv')){
    stalage_chrono <- read.csv('StalAge_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>%
      distinct(sample_id, .keep_all = T) %>% mutate(sample_id = paste(sample_id,'-',file_name,sep=''))
    names(stalage_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
    SISAL_chronology_new <- full_join(stalage_chrono, SISAL_chronology_new, by = 'sample_id')  %>%
      mutate(StalAge_age = if_else(sample_id %in% stalage_chrono$sample_id, age, StalAge_age),
             StalAge_age_uncert_pos = if_else(sample_id %in% stalage_chrono$sample_id, uncert_pos, StalAge_age_uncert_pos),
             StalAge_age_uncert_neg = if_else(sample_id %in% stalage_chrono$sample_id, uncert_neg, StalAge_age_uncert_neg)) %>%
      select(-age, -uncert_pos, -uncert_neg)
  }

  setwd(file.path(working_directory, file_name, '/Bacon_runs'))
  if(file.exists('bacon_chronology.csv')){
    bacon_chrono <- read.csv('bacon_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>%
      distinct(sample_id, .keep_all = T) %>% mutate(sample_id = paste(sample_id,'-',file_name,sep=''))
    names(bacon_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
    SISAL_chronology_new <- full_join(bacon_chrono, SISAL_chronology_new, by = 'sample_id') %>%
      mutate(bacon_age = if_else(sample_id %in% bacon_chrono$sample_id, age, bacon_age),
             bacon_age_uncert_pos = if_else(sample_id %in% bacon_chrono$sample_id, uncert_pos, bacon_age_uncert_pos),
             bacon_age_uncert_neg = if_else(sample_id %in% bacon_chrono$sample_id, uncert_neg, bacon_age_uncert_neg)) %>%
      select(-age, -uncert_pos, -uncert_neg)
  }
  return(SISAL_chronology_new)
}

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

#synSpeleo_chrono <- quanAna(p.min=0.1,p.max=300,p.n=7, n=2)
synSpeleo_chrono_AM_3 <- quanAna_AM(file_path1 = '/stacycode/ariana/synSpeleo/AM',file_path2='/stacycode/ariana/synSpeleo/summary/', precision = c(10,50,100,200,300,500), n=100,h.length=c(1000,2000,500), k = 3)

save(synSpeleo_chrono_AM_3, file = '/stacycode/ariana/synSpeleo/summary/synSpeleo_chrono_AM_3.RData')

