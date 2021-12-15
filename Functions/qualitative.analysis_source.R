##### quantitative analysis #####

setwd("/home/ariana/R/SISAL_AM")
source(file = 'run_functions.R')
source(file = 'SISAL_AM_add_functions.R')
source(file = 'SISAL_AM_plot_functions.R')
source(file = "StalAge_1_0_mine.R")
source(file = 'New_methods.R')
source(file = 'lin_interp_reversal.R')
source(file = '~/R/syntheticspeleo-master/R/functions.R')
source(file = '~/R/syntheticspeleo-master/R/validation_fct.R')

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

quanAna <- function(p.min = 0.1, p.max=1,p.n=5,t.min=0,t.max=10000,t.n=4,d.min=0,d.max=1000,h.length=c(2000,3000,500),file_path1 ='~/synSpeleo/test/', file_path2 = '~/synSpeleo/test/summary' ,rate=3,t.res=2, nsamp=15,n=100){
  
  set.seed(2020)
  chrono_final <- tibble()
  for(k in 0:length(h.length)){
    
    
    if(k==0){
      t <- t.max-t.min
      l <- 0
    } else {
      h.n <- h.length[1:k]
      t <- t.max-t.min - sum(h.n)
      l <- length(h.n)
    }
    
    
    
    prec <- round(cumsum(c(0.1,rep((p.max-p.min)/p.n,p.n))),digits=1)
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
          v <- sample(length(h.n)+t.n)
          tlist2 <- c(t.min,c(cumsum(c(h.n,rep(t/t.n,t.n))[v])))
          rates2 <- c(rep(0,l),rep(rate,t.n))[v]
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
        
        chrono <- merge_synSpeleo_qa(file_name, file_path1, sample) %>% arrange(., depth_sample)
        
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
                 bacon_age_uncert_neg = if_else(bacon_age < -68, NA_real_,bacon_age_uncert_neg),
                 bchron_age = if_else(bchron_age < -68, NA_real_,bchron_age),
                 bchron_age_uncert_pos = if_else(bchron_age < -68, NA_real_,bchron_age_uncert_pos),
                 bchron_age_uncert_neg = if_else(bchron_age < -68, NA_real_,bchron_age_uncert_neg))
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
                                                bacon_age_uncert_neg = if_else(depth_sample > d.age.na, NA_real_,bacon_age_uncert_neg),
                                                bchron_age = if_else(depth_sample > d.age.na, NA_real_,bchron_age),
                                                bchron_age_uncert_pos = if_else(depth_sample > d.age.na, NA_real_,bchron_age_uncert_pos),
                                                bchron_age_uncert_neg = if_else(depth_sample > d.age.na, NA_real_,bchron_age_uncert_neg))
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
                                                bacon_age_uncert_neg = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,bacon_age_uncert_neg),
                                                bchron_age = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,bchron_age),
                                                bchron_age_uncert_pos = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,bchron_age_uncert_pos),
                                                bchron_age_uncert_neg = if_else(depth_sample > d.age.2000[1] & depth_sample < d.age.2000[2], NA_real_,bchron_age_uncert_neg))
            
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
                                                  bacon_age_uncert_neg = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,bacon_age_uncert_neg),
                                                  bchron_age = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,bchron_age),
                                                  bchron_age_uncert_pos = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,bchron_age_uncert_pos),
                                                  bchron_age_uncert_neg = if_else(depth_sample > d.age.2000[3] & depth_sample < d.age.2000[4], NA_real_,bchron_age_uncert_neg))
            }
          }
        }
        chrono_final <- bind_rows(chrono_final,chrono_new)
        save(chrono_new, file = file.path(file_path2, paste('chrono-',l,'-',j,'-',i,'.RData',sep='')))
        
      }
    }
  }
  return(chrono_final)
}


merge_synSpeleo_qa <- function(file_name, working_directory, se){
  
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
                                     bchron_age_uncert_neg = rep(NA_real_, length(se$sample_id)))%>%
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

#synSpeleo_chrono <- quanAna(p.min=0.1,p.max=300,p.n=7, n=2)
synSpeleo_chrono <- quanAna(file_path1 = '/stacycode/ariana/synSpeleo/',file_path2='/stacycode/ariana/synSpeleo/summary/', p.min=0.1,p.max=300,p.n=5, n=10,h.length=c(1000,2000,500))

save(synSpeleo_chrono, file = '/stacycode/ariana/synSpeleo/summary/synSpeleo_chrono.RData')


names <- rep(NA_real_,2400)
counter <- 0
for(k in 0:3){
  prec <- c(10,50,100,200,300,500)
  
  for (j in prec){
    
    for (i in 1:100) {
      if(k==0){
        
      }
      counter <- counter+1
      if(k==0){
        name <- paste('/stacycode/ariana/synSpeleo_neu/summary/chrono-',k,'-',j,'-',i,'_AM.RData',sep='')
        names[counter] <- name
      } else {
        name <- paste('/stacycode/ariana/synSpeleo/summary/chrono-',k,'-',j,'-',i,'_AM.RData',sep='')
        names[counter] <- name
      }
    }
  }
}

prec_hiatus <- tibble(enity_id = character(),hiatus_nr = numeric(),prec = numeric(), iteration = numeric() ,lr = numeric(),li=numeric(),copRa=numeric(),bacon=numeric(),stalage=numeric())
for(name in names){
  
  #load(file.path('/stacycode/ariana/synSpeleo/summary',name))
  load(file.path(name))
  
  chrono_sum <- chrono_new %>% summarise(entity_id = unique(entity_id), hiatus_nr = unique(hiatus_nr), iteration = unique(iteration), prec = unique(precision),
                                         lr = sum(true.age > lin_reg_age -lin_reg_age_uncert_neg & true.age < lin_reg_age + lin_reg_age_uncert_pos, na.rm = T)/(sum(!is.na(lin_reg_age))),
                                         li = sum(true.age > lin_interp_age -lin_interp_age_uncert_neg & true.age < lin_interp_age + lin_interp_age_uncert_pos, na.rm = T)/(sum(!is.na(lin_interp_age))),
                                         copRa = sum(true.age > copRa_age - copRa_age_uncert_neg & true.age < copRa_age + copRa_age_uncert_pos, na.rm = T)/(sum(!is.na(copRa_age))),
                                         bacon = sum(true.age > bacon_age -bacon_age_uncert_neg & true.age < bacon_age + bacon_age_uncert_pos, na.rm = T)/(sum(!is.na(bacon_age))),
                                         stalage = sum(true.age > StalAge_age -StalAge_age_uncert_neg & true.age < StalAge_age + StalAge_age_uncert_pos, na.rm = T)/(sum(!is.na(StalAge_age))))
  
  prec_hiatus <- bind_rows(prec_hiatus, chrono_sum)
}

p <- (prec_hiatus %>%  distinct(prec))$prec

mtx_i <- list(matrix(),matrix(),matrix(),matrix(),matrix())
#am <- c('lr','li','copRa','bacon','stalage')
for(f in 1:5){
  mtx <- matrix(NA_real_,60,40)
  j <- 0
  l <- 0
  for (l in 1:6){
    i <- p[l]
    #n <- 1
    for(h in 0:3){
      j <- h+1
      m <- prec_hiatus %>% filter(hiatus_nr == h & prec == i) %>% arrange(.,iteration) %>% pull(f+4)
      mtx_r <- matrix(m,nrow=10,ncol=10,byrow = T)
      mtx[((l-1)*10)+1:10,((j-1)*10)+1:10] <- mtx_r
    }
  }
  
  mtx_i[[f]] <- mtx
}

prec_hiatus_new <- prec_hiatus %>% replace(is.na(.),0) %>% group_by(prec,hiatus_nr) %>% summarise(lr.prec = mean(lr, na.rm = T),
                                                                                                  li.prec = mean(li, na.rm = T),
                                                                                                  copRa.prec = mean(copRa, na.rm = T),
                                                                                                  stalage.prec = mean(stalage, na.rm = T),
                                                                                                  bacon.prec = mean(bacon, na.rm = T))

mean(prec_hiatus_new$lr.prec)
mean(prec_hiatus_new$li.prec)
mean(prec_hiatus_new$copRa.prec)
mean(prec_hiatus_new$stalage.prec)
mean(prec_hiatus_new$bacon.prec)

col_ipcc_sequential_blue <- function(n) {
  colorRampPalette(c(rgb(237,248,251,maxColorValue = 255),
                     rgb(179,205,227,maxColorValue = 255),
                     rgb(140,150,198,maxColorValue = 255),
                     rgb(136,86,167,maxColorValue = 255),
                     rgb(129,15,124,maxColorValue = 255))[c(1,1,1,1,1,1,1,1,2,3,4,5)],space="rgb")(n)
}


# mtx_lr <- mtx
# mtx_li <- mtx
# mtx_bacon <- mtx
# mtx_stalage <- mtx
# mtx_copRa <- mtx
#
#
# mtx_i <- list(mtx_lr,mtx_li,mtx_copRa,mtx_stalage,mtx_bacon)

title <- c('a) Linear Regression','b) Linear Interpolation','c) copRa','d) StalAge','e) Bacon')

if(T){
  pdf('~/synSpeleo/plots/eval_test_new_new.pdf',8,12)
  par(mfrow = c(3,2),mar=c(4,4,2,0.5))
  
  for(i in 1:5){
    image(1:60,1:40,mtx_i[[i]], col = col_ipcc_sequential_blue(230), xlab = '', ylab = '', xaxt = "n", yaxt="n", zlim=c(0,1))
    image(1:60,1:40,array(is.na(mtx_i[[i]]),dim=c(60,40)),add=T,col="black",zlim=c(0.5,1.5))
    #image(1:60,1:40,array(is.na(mtx_i[[i]]),dim=c(60,40)),add=T,col="darkgrey",zlim=c(0.5,1.5))
    mtext(text='Number of Hiatuses', side = 2, line = 2, cex = 0.8)
    mtext(text='Dating Uncertainty', side = 1, line = 3, cex = 0.8)
    axis(side = 1,at = seq(5.5,55.5,10), labels = p,tick = F)
    axis(side = 1,at = seq(10.5,50.5,10), labels = rep('',5),tick = T, lwd.ticks = 2)
    abline(v=10*(1:9)+0.5, lwd = 2)
    axis(side = 2,at = seq(5.5,35.5,10), labels = seq(0,3,1),tick = F)
    axis(side = 2,at = seq(10.5,30.5,10), labels = rep('',3),tick = T, lwd.ticks = 2)
    box(lwd = 2)
    abline(h=10*(1:3)+0.5, lwd=2)
    title(title[i],adj = 0 )
  }
  plot.new()
  fields::image.plot(1:60,1:40,mtx_i[[5]], col = col_ipcc_sequential_blue(230), zlim = c(0,1),
                     legend.only=T,horizontal = T,smallplot = c(0.15,0.85,0.55,0.65))
  text(x = 0.43,y=0.35, labels = expression(bold('Coverage Frequency')),cex=2)
  dev.off()
}



p <- (prec_hiatus %>%  distinct(prec))$prec
mtx <- matrix(NA_real_,70,40)
j <- 0
l <- 0
for (l in 2:8){
  i <- p[l]
  #n <- 1
  for(h in 0:3){
    j <- h+1
    m <- (prec_hiatus %>% filter(hiatus_nr == h & prec == i) %>% arrange(.,iteration) %>% select(copRa))$copRa
    mtx_r <- matrix(m,nrow=10,ncol=10,byrow = T)
    mtx[((l-2)*10)+1:10,((j-1)*10)+1:10] <- mtx_r
  }
}

mtx_copRa <- mtx

fields::image.plot(1:60,1:40,mtx_copRa, col = col_ipcc_sequential_blue(230), xlab = 'Precision', ylab = '#hiatus', xaxt = "n", yaxt="n")
image(1:60,1:40,array(is.na(mtx_copRa),dim=c(60,40)),add=T,col="black",zlim=c(0.5,1.5))
axis(side = 1,at = seq(5.5,75.5,10), labels = p,tick = F)
axis(side = 1,at = seq(10.5,70.5,10), labels = rep('',7),tick = T, lwd.ticks = 2)
abline(v=10*(1:9)+0.5, lwd = 2)
axis(side = 2,at = seq(5.5,35.5,10), labels = seq(0,3,1),tick = F)
axis(side = 2,at = seq(10.5,30.5,10), labels = rep('',3),tick = T, lwd.ticks = 2)
box(lwd = 2)
abline(h=10*(1:3)+0.5, lwd=2)
mtext(side = 3, text = expression(bold('a) Linear Regression')),at = 10 )

l <- 1
j <- 1
i <- 0.1
h <- 0
































