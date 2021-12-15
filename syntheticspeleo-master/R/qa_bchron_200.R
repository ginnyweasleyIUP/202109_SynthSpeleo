##### quantitative analysis #####

setwd("/home/ariana/R/SISAL_AM")
source(file = 'run_functions.R')
source(file = 'SISAL_AM_add_functions.R')
source(file = 'SISAL_AM_plot_functions.R')
source(file = "StalAge_1_0_mine.R")
source(file = 'New_methods.R')
source(file = 'lin_interp_reversal.R')
#source(file = '~/R/syntheticspeleo-master/R/functions.R')
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

quanAna_bchron <- function(t.min=0,t.max=10000,t.n=4,d.min=0,d.max=1000,h.length=c(2000,3000,500),file_path1 ='~/synSpeleo/test/', file_path2 = '~/synSpeleo/test/summary' ,rate=3,t.res=2, nsamp=15,p){

  set.seed(2020)
  chrono_final <- tibble()

  for (k in 0:3){
    if(k==0){
      t <- t.max-t.min
      l <- 0
    } else {
      h.n <- h.length[1:k]
      t <- t.max-t.min - sum(h.n)
      l <- length(h.n)
    }



    j <- p
    #if(k ==0){h.n <- tibble(sample_id = numeric(), depth_sample = numeric())}

    # go through precisions

    print(j)

    n <- sample(1:100,10)

    for (i in n) {
      print(i)

      # create speleo object
      load(file.path(file_path2, paste('grwthfct-',l,'-',j,'-',i,'.RData',sep='')))

      hiatus <- speleo_new$hiatus
      #dind<-check.Dtable(speleo_class1,add=TRUE)
      #save(speleo_new,file = file.path(file_path2, paste('grwthfct-',l,'-',j,'-',i,'.RData',sep='')))
      run_SISAL_chrono_bchron(speleo_new, file_path1,paste(l,'-',j,'-',i,sep = '') , hiatus)

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

      chrono <- merge_synSpeleo_qa_bchron(file_name, file_path1, sample) %>% arrange(., depth_sample)

      ind.samples<-sapply(chrono$depth_sample,function(x,growthfcn){which.min(abs(growthfcn$growth.real-x))},growthfcn=speleo_new$growthfcn)
      true.age<-speleo_new$growthfcn$true.age[ind.samples]
      #true.age <- grf$true_age[sapply(chrono$depth_sample,function(x) which.min(abs(x-grf$expected.z)))]

      speleo_new$recorded.age <- true.age



      chrono_new <- chrono %>% mutate(true.age = true.age, precision = rep(j, length(chrono$depth_sample)),
                                      hiatus_nr = rep(k, length(chrono$depth_sample)), iteration = rep(i, length(chrono$depth_sample))) %>%
        mutate(bchron_age = if_else(bchron_age < -68, NA_real_,bchron_age),
               bchron_age_uncert_pos = if_else(bchron_age < -68, NA_real_,bchron_age_uncert_pos),
               bchron_age_uncert_neg = if_else(bchron_age < -68, NA_real_,bchron_age_uncert_neg))
      if(k!=0){
        if(!is.na(d.age.na)){
          chrono_new <- chrono_new %>% mutate(bchron_age = if_else(depth_sample > d.age.na, NA_real_,bchron_age),
                                              bchron_age_uncert_pos = if_else(depth_sample > d.age.na, NA_real_,bchron_age_uncert_pos),
                                              bchron_age_uncert_neg = if_else(depth_sample > d.age.na, NA_real_,bchron_age_uncert_neg))
        }

        if(!is.na(d.age.2000)){
          chrono_new <- chrono_new %>% mutate( bchron_age = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,bchron_age),
                                               bchron_age_uncert_pos = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,bchron_age_uncert_pos),
                                               bchron_age_uncert_neg = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,bchron_age_uncert_neg))

          if(length(d.age.2000) > 2) {
            chrono_new <- chrono_new %>% mutate( bchron_age = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,bchron_age),
                                                 bchron_age_uncert_pos = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,bchron_age_uncert_pos),
                                                 bchron_age_uncert_neg = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,bchron_age_uncert_neg))
          }
        }
      }
      #chrono_final <- bind_rows(chrono_final,chrono_new)
      save(chrono_new, file = file.path(file_path2, paste('chrono-',l,'-',j,'-',i,'_bchron.RData',sep='')))

    }
  }
  #return(chrono_final)
}


merge_synSpeleo_qa_bchron <- function(file_name, working_directory, se){

  SISAL_chronology_new <- data.frame(entity_id = se$entity_id,
                                     sample_id = se$sample_id,
                                     depth_sample = se$depth_sample,
                                     bchron_age = rep(NA_real_, length(se$sample_id)),
                                     bchron_age_uncert_pos = rep(NA_real_, length(se$sample_id)),
                                     bchron_age_uncert_neg = rep(NA_real_, length(se$sample_id)))%>%
    mutate_at(vars(entity_id,sample_id),as.character) %>% arrange(., depth_sample)




  setwd(file.path(working_directory, file_name, '/Bchron'))
  if(file.exists('bchron_chronology.csv')){
    bchron_chrono <- read.csv('bchron_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>%
      distinct(sample_id, .keep_all = T) %>% mutate(sample_id = paste(sample_id,'-',file_name,sep=''))
    names(bchron_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
    SISAL_chronology_new <- full_join(bchron_chrono, SISAL_chronology_new, by = 'sample_id') %>%
      mutate(bchron_age = if_else(sample_id %in% bchron_chrono$sample_id, age, bchron_age),
             bchron_age_uncert_pos = if_else(sample_id %in% bchron_chrono$sample_id, uncert_pos, bchron_age_uncert_pos),
             bchron_age_uncert_neg = if_else(sample_id %in% bchron_chrono$sample_id, uncert_neg, bchron_age_uncert_neg)) %>%
      select(-age, -uncert_pos, -uncert_neg)
  }

  return(SISAL_chronology_new)
}

run_SISAL_chrono_bchron <- function(speleo, working_directory, file_name, hiatus) {
  #saveSynSpeleo(speleo, working_directory, file_name, hiatus)

  err <- NULL

  tryCatch({
    i = 3
    runBchron_new(working_directory, file_name)
    print('Bchron done')},
    error = function(e){err <<- append(err, paste("ERROR in Bchron:",conditionMessage(e), "\n"))
    #error = function(e){print(paste("ERROR in Bchron new:",conditionMessage(e), "\n"))
    #runFile[j,i] <<- F
    bchron <<- F})

  setwd(file.path(working_directory, file_name))
  print('wd changed')

  write.table(err, 'errors_bchron.txt', row.names = F, col.names = F)

  setwd(file.path(working_directory))

  print(paste(file_name, 'done'))
  gc()
}

#synSpeleo_chrono <- quanAna(p.min=0.1,p.max=300,p.n=7, n=2)
quanAna_bchron(file_path1 = '/stacycode/ariana/synSpeleo/AM',file_path2='/stacycode/ariana/synSpeleo/summary/',h.length=c(1000,2000,500), p = 200)

#save(synSpeleo_chrono_bchron_3, file = '/stacycode/ariana/synSpeleo/summary/synSpeleo_chrono_bchron_3.RData')
