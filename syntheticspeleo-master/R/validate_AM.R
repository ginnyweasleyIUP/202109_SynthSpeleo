source("~/R/syntheticspeleo-master/R/functions.R")
library(zoo)
library(Hmisc)
library(plyr)
library(pracma)

# WORKFLOW
set.seed(0)
speleo<-NA
speleo<-generate.speleo(t.top=0,t.bottom=10000,z.top = 0,z.bottom=1000,t.res=2.5)

### class 1: slight growth rate changes, no hiatus, no 'major' outliers
speleo_class1<-growthfcn.changinggamma.new(speleo, tlist = c(0,2000,5000,10000), rates = c(4,2,1))

speleo_class1$datingparams$tdofs<-10
speleo_class1$datingparams$precision<-100
speleo_class1$datingparams$topage<-10
speleo_class1$datingparams$nsample<-10

speleo_class1<-do.dating(speleo_class1)
speleo_class1$Dtable

plotgrowthfcn(speleo_class1)
plotdating(speleo_class1$Dtable,add=TRUE)
speleo_class1<-sample.proxy(speleo = speleo_class1,nsamp = 70,z.steps = 1,stype = "sinusoid")
dind<-check.Dtable(speleo_class1,add=TRUE)
run_SISAL_chrono(speleo_class1, '~/synSpeleo/', 'class1', tibble(sample_id = numeric(), depth_sample = numeric()))

save(speleo_class1,file = '~/synSpeleo/speleo_class1.RData')

### class 2: slight growth rate changes, no hiatus, no 'major' outliers
speleo_class2<-growthfcn.changinggamma.new(speleo, tlist = c(0,500,4000,8000,10000), rates = c(4,1,0,2))

speleo_class2$datingparams$tdofs<-10
speleo_class2$datingparams$precision<-150
speleo_class2$datingparams$topage<-15
speleo_class2$datingparams$nsample<-10

speleo_class2<-do.dating(speleo_class2)
speleo_class2$Dtable

plotgrowthfcn(speleo_class2)
plotdating(speleo_class2$Dtable,add=TRUE)
speleo_class2<-sample.proxy(speleo = speleo_class2,nsamp = 70,z.steps = 1,stype = "sinusoid")
dind<-check.Dtable(speleo_class2,add=TRUE)

hiatus2 <- unique(speleo_class2$growthfcn$expected.z[which(diff(speleo_class2$growthfcn$expected.z)==0)])

run_SISAL_chrono(speleo_class2, '~/synSpeleo/', 'class2', tibble(sample_id = 150, depth_sample = hiatus))
save(speleo_class2,file = '~/synSpeleo/speleo_class2.RData')
### class 3: slight growth rate changes, no hiatus, no 'major' outliers
speleo_class3<-growthfcn.changinggamma.new(speleo, tlist = c(0,500,1000,2000,2800,4500,7300,10000), rates = c(1,4,0,2,5,0,1))
speleo_class3$datingparams$nfac <-300
speleo_class3$datingparams$tdofs<-5000
speleo_class3$datingparams$precision<-150
speleo_class3$datingparams$topage<-200
speleo_class3$datingparams$nsample<-20

speleo_class3<-do.dating(speleo_class3)
speleo_class3$Dtable

plotgrowthfcn(speleo_class3)
plotdating(speleo_class3$Dtable_new,add=TRUE)
speleo_class3<-sample.proxy(speleo = speleo_class3,nsamp = 70,z.steps = 1,stype = "sinusoid")
dind<-check.Dtable(speleo_class2,add=TRUE)

hiatus3 <- unique(speleo_class3$growthfcn$expected.z[which(diff(speleo_class3$growthfcn$expected.z)==0)])


#speleo_class_prim <- speleo_class3
speleo_class3$Dtable_new <- speleo_class3$Dtable
speleo_class3$Dtable_new[8,] <- NA
speleo_class3$Dtable_new[10,3] <- speleo_class3$Dtable[10,3]+400
speleo_class3$Dtable_new[16,3] <- speleo_class3$Dtable[16,3]-800
speleo_class3$Dtable_new[17,3] <- speleo_class3$Dtable[17,3]-500
speleo_class3$Dtable_new[13,] <- NA
speleo_class3$Dtable_new[6,2] <- speleo_class3$Dtable[6,2]-20
speleo_class3$Dtable_new[3,2] <- speleo_class3$Dtable[3,2]-40
speleo_class3$Dtable_new[20,2] <- speleo_class3$Dtable[20,2]-30
speleo_class3$Dtable_new[21,2] <- speleo_class3$Dtable[21,2]+40



run_SISAL_chrono(speleo_class3, '~/synSpeleo/', 'class3', tibble(sample_id = c(150,151), depth_sample = hiatus))
save(speleo_class3,file = '~/synSpeleo/speleo_class3.RData')


speleo_class4 = speleo_class3
speleo_class4$Dtable <- speleo_class3$Dtable_new %>% filter(!is.na(dating_id))
hiatus4 <- unique(speleo_class4$growthfcn$expected.z[which(diff(speleo_class4$growthfcn$expected.z)==0)])
run_SISAL_chrono(speleo_class4, '~/synSpeleo/', 'class4', tibble(sample_id = c(150,151), depth_sample = hiatus))
save(speleo_class4,file = '~/synSpeleo/speleo_class4.RData')
load('~/synSpeleo/class4.RData')
load('~/synSpeleo/class3.RData')
load('~/synSpeleo/class2.RData')
load('~/synSpeleo/class1.RData')

runFile <- tibble(entity_id = c('class1','class2','class3','class4'))
dating <- bind_rows(speleo_class1$Dtable %>% mutate(entity_id = 'class1'),
                    speleo_class2$Dtable %>% mutate(entity_id = 'class2', dating_id = dating_id+1000),
                    tibble(dating_id = 25000, sampling.depths = hiatus2, sampling.age.est = NA_real_, sampling.age.est.unc = NA_real_, sampling.age.true=NA_real_, entity_id = 'class2'),
                    speleo_class3$Dtable %>% mutate(entity_id = 'class3', dating_id = dating_id+2000),
                    tibble(dating_id = c(35000,35100), sampling.depths = hiatus3, sampling.age.est = c(NA_real_,NA_real_), sampling.age.est.unc = c(NA_real_,NA_real_), sampling.age.true=c(NA_real_,NA_real_), entity_id = c('class3','class3')),
                    speleo_class4$Dtable %>% mutate(entity_id = 'class4', dating_id = dating_id+3000),
                    tibble(dating_id = c(45000,45100), sampling.depths = hiatus3, sampling.age.est = c(NA_real_,NA_real_), sampling.age.est.unc = c(NA_real_,NA_real_), sampling.age.true=c(NA_real_,NA_real_), entity_id = c('class4','class4'))) %>%
  group_by(entity_id) %>% arrange(.,sampling.depths, .by_group=T) %>% mutate(date_type = 'synthetic Age')

sample <- bind_rows(speleo_class1$proxy$proxymat %>% mutate(entity_id = 'class1') %>% select(sample_id, depth_sample, proxy, agetrue, width, entity_id),
                    speleo_class2$proxy$proxymat %>% mutate(entity_id = 'class2', sample_id = sample_id+200) %>% select(sample_id, depth_sample, proxy, agetrue, width,entity_id),
                    tibble(sample_id = 25000, depth_sample = hiatus2, proxy = NA_real_, agetrue = NA_real_, width =NA_real_, entity_id = 'class2'),
                    speleo_class3$proxy$proxymat %>% mutate(entity_id = 'class3', sample_id = sample_id+300) %>% select(sample_id, depth_sample, proxy, agetrue, width, entity_id),
                    tibble(sample_id = c(35000,35100), depth_sample = hiatus3, proxy = c(NA_real_,NA_real_), agetrue = c(NA_real_,NA_real_), width =c(NA_real_,NA_real_), entity_id = c('class3','class3')),
                    speleo_class4$proxy$proxymat %>% mutate(entity_id = 'class4', sample_id = sample_id+400) %>% select(sample_id, depth_sample, proxy, agetrue, width, entity_id),
                    tibble(sample_id = c(45000,45100), depth_sample = hiatus3, proxy = c(NA_real_,NA_real_), agetrue = c(NA_real_,NA_real_), width=c(NA_real_,NA_real_), entity_id = c('class4','class4'))) %>%
  group_by(entity_id) %>% arrange(.,depth_sample, .by_group=T)

hiatus_tb <- bind_rows(tibble(dating_id = 25000, sampling.depths = hiatus2, entity_id = 'class2'),
                       tibble(dating_id = c(35000,35100), sampling.depths = hiatus3, entity_id = c('class3','class3')),
                       tibble(dating_id = c(45000,45100), sampling.depths = hiatus3, entity_id = c('class4','class4')))

load('~/synSpeleo/')

synSpeleo_chrono <-





setwd("/home/ariana/R/SISAL_AM")
source(file = 'run_functions.R')
source(file = 'SISAL_AM_add_functions.R')
source(file = 'SISAL_AM_plot_functions.R')
source(file = "StalAge_1_0_mine.R")
source(file = 'New_methods.R')
source(file = 'lin_interp_reversal.R')




runFile2 <- tibble(entity_id = c('class1','class2','class3','class4'), add = c(0,200,300,400))

synSpeleo_chrono <- SISAL_chronology_new %>% mutate_at(vars(entity_id),as.character) %>% filter(!is.na(entity_id)) %>% group_by(entity_id) %>% arrange(., depth_sample, .by_group = T)
synSpeleo_chrono[1:70,3] <-  speleo_class1$proxy$proxymat$depth_sample

#synSpeleo_chrono_new_new <- SISAL_chronology_new %>% mutate_at(vars(entity_id),as.character()) %>% mutate(depth_sample = if_else(entity_id == 'class1', speleo_class1$proxy$proxymat$depth_sample, depth_sample))

idx2 <- sample %>% filter(entity_id == 'class2' & depth_sample > 545 & depth_sample < 637)
idx3 <- sample %>% filter(entity_id == 'class3' & depth_sample > 809 & depth_sample < 858)
idx4 <- sample %>% filter(entity_id == 'class4' & depth_sample > 809 & depth_sample < 858)

synSpeleo_chrono_new <- synSpeleo_chrono %>% mutate(lin_reg_age = if_else(sample_id %in% c(idx2$sample_id,idx3$sample_id,idx4$sample_id) , NA_real_,lin_reg_age),
                                                    lin_reg_age_uncert_pos = if_else(sample_id %in% c(idx2$sample_id,idx3$sample_id,idx4$sample_id), NA_real_,lin_reg_age_uncert_pos),
                                                    lin_reg_age_uncert_neg = if_else(sample_id %in% c(idx2$sample_id,idx3$sample_id,idx4$sample_id), NA_real_,lin_reg_age_uncert_neg),
                                                    lin_interp_age = if_else(sample_id %in% c(idx2$sample_id,idx3$sample_id,idx4$sample_id), NA_real_,lin_interp_age),
                                                    lin_interp_age_uncert_pos = if_else(sample_id %in% c(idx2$sample_id,idx3$sample_id,idx4$sample_id), NA_real_,lin_interp_age_uncert_pos),
                                                    lin_interp_age_uncert_neg = if_else(sample_id %in% c(idx2$sample_id,idx3$sample_id,idx4$sample_id), NA_real_,lin_interp_age_uncert_neg),
                                                    copRa_age = if_else(sample_id %in% c(idx2$sample_id,idx3$sample_id,idx4$sample_id), NA_real_,copRa_age),
                                                    copRa_age_uncert_pos = if_else(sample_id %in% c(idx2$sample_id,idx3$sample_id,idx4$sample_id), NA_real_,copRa_age_uncert_pos),
                                                    copRa_age_uncert_neg = if_else(sample_id %in% c(idx2$sample_id,idx3$sample_id,idx4$sample_id), NA_real_,copRa_age_uncert_neg),
                                                    StalAge_age = if_else(sample_id %in% c(idx2$sample_id,idx3$sample_id,idx4$sample_id), NA_real_,StalAge_age),
                                                    StalAge_age_uncert_pos = if_else(sample_id %in% c(idx2$sample_id,idx3$sample_id,idx4$sample_id), NA_real_,StalAge_age_uncert_pos),
                                                    StalAge_age_uncert_neg = if_else(sample_id %in% c(idx2$sample_id,idx3$sample_id,idx4$sample_id), NA_real_,StalAge_age_uncert_neg),
                                                    bacon_age = if_else(sample_id %in% c(idx2$sample_id,idx3$sample_id,idx4$sample_id), NA_real_,bacon_age),
                                                    bacon_age_uncert_pos = if_else(sample_id %in% c(idx2$sample_id,idx3$sample_id,idx4$sample_id), NA_real_,bacon_age_uncert_pos),
                                                    bacon_age_uncert_neg = if_else(sample_id %in% c(idx2$sample_id,idx3$sample_id,idx4$sample_id), NA_real_,bacon_age_uncert_neg),
                                                    bchron_age = if_else(sample_id %in% c(idx2$sample_id,idx3$sample_id,idx4$sample_id), NA_real_,bchron_age),
                                                    bchron_age_uncert_pos = if_else(sample_id %in% c(idx2$sample_id,idx3$sample_id,idx4$sample_id), NA_real_,bchron_age_uncert_pos),
                                                    bchron_age_uncert_neg = if_else(sample_id %in% c(idx2$sample_id,idx3$sample_id,idx4$sample_id), NA_real_,bchron_age_uncert_neg)) %>%
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

save(synSpeleo_chrono_new, file ='~/synSpeleo/synSpeleo_chrono.RData')

merge_synSpeleo(runFile, dating, hiatus_tb,'~/synSpeleo/')
correct_am_synSpeleo(synSpeleo_chrono_new, '~/synSpeleo/', '~/synSpeleo/', tibble(entity_id = c('class1','class2','class3','class4')))
correct_ens_synSpeleo(tibble(entity_id = c('class1','class2','class3','class4')), '~/synSpeleo/',hiatus_tb)

speleoList <- list()
speleoList$one <- speleo_class1
speleoList$two <- speleo_class2
speleoList$three <- speleo_class3
speleoList$four <- speleo_class4

plot_synSpeleo(tibble(entity_id = c('class1','class2','class3','class4')), 8,12, hiatus_tb, speleoList)

plotsynSpeleo_method <- function(ens, AM, dating, depth, text,stalage=F){
  if(stalage){
    plot(y = depth, x = AM$age, col = 'darkblue', type = 'l', ylim = c(max(depth, na.rm = T),0), xlab = '', ylab = '')
    lines(y = depth, x = AM$age-AM$age_uncert_neg, col = 'red', type = 'l', lty = 2)
    lines(y = depth, x = AM$age+AM$age_uncert_pos, col = 'red', type = 'l', lty = 2)
  } else {
    matplot(y = ens[,2], x = ens[,3:1000], col = alpha('grey',0.2), type = 'l', ylim = c(max(ens[,2], na.rm = T),0), xlab = '', ylab = '')
    if(is.na(depth)){
      lines(y = ens[,2], x = AM$age, col = 'darkblue', type = 'l')
      lines(y = ens[,2], x = AM$age-AM$age_uncert_neg, col = 'red', type = 'l', lty = 2)
      lines(y = ens[,2], x = AM$age+AM$age_uncert_pos, col = 'red', type = 'l', lty = 2)
    } else {
      lines(y = depth, x = AM$age, col = 'darkblue', type = 'l')
      lines(y = depth, x = AM$age-AM$age_uncert_neg, col = 'red', type = 'l', lty = 2)
      lines(y = depth, x = AM$age+AM$age_uncert_pos, col = 'red', type = 'l', lty = 2)
    }
  }
  #lines(y = proxy$depth_sample, x = proxy$d18O_measurement, col = 'black')
  points(x=dating$sampling.age.est,y=dating$sampling.depths, col = 'black', pch = 19)
  arrows(dating$sampling.age.est-dating$sampling.age.est.unc, dating$sampling.depths,dating$sampling.age.est+dating$sampling.age.est.unc, dating$sampling.depths,length=0.05, angle=90, code=3, col = 'black')
  mtext(text, side = 3, line = 0)
  mtext('Depth [a.u.]', side = 2, line = 2)
  mtext('Years [a.u.]', side = 1, line = 3)
}


plotTrueAge<-function(speleo){
  orig<-zoo::zoo(speleo$growthfcn$expected.z,order.by=speleo$growthfcn$true.age)
  sim<-zoo::zoo(speleo$growthfcn$growth.real,order.by=speleo$growthfcn$true.age)
  orig[diff(speleo$growthfcn$growth.real)==0]<-NA
  sim[diff(speleo$growthfcn$growth.real)==0]<-NA
  lines(sim,col="black")
}

plot_synSpeleo <- function(eID, width, height, hiatus_tb,speleo){
  j <- 0
  for(i in eID$entity_id) {
    j <- j+1
    print(i)
    entity_name <- i

    #print(j)

    speleo_new <- speleo[[j]]
    #print(speleo)

    h <- hiatus_tb %>% filter(entity_id == i)

    file_name <- get(load(file.path('~/synSpeleo/', paste(entity_name, '.RData', sep = ''))))
    dating <- file_name$dating %>% filter(entity_id == i)
    proxy <- file_name$proxy %>% select(sample_id, depth_sample,d18O_measurement)

    proxy_new <- bind_rows(proxy,tibble(entity_id = h$dating_id,depth_sample=h$sampling.depths,d18O_measurement = rep(NA,dim(h)[1]))) %>% arrange(., depth_sample)



    if(i %in% c('class1','class2')) {
      pdf(paste('~/synSpeleo/',i,'.pdf', sep=''), width, height)
      layout(matrix(c(1,2,3,4,5,6), nrow=3, ncol=2, byrow = T))}
    if(i %in% c('class3','class4')) {
      pdf(paste('~/synSpeleo/',i,'.pdf', sep=''), width, height*2/3)
      layout(matrix(c(1,2,3,4), nrow=2, ncol=2, byrow = T))}
    #layout(matrix(c(1,2,3,4,5,6), nrow=3, ncol=2, byrow = T))
    par(mar = c(4,3,3,4), oma = c(1,3,2,1))

    #print('0')


    if(!plyr::empty(data.frame(file_name$linReg$ens))){
      ens <- file_name$linReg$ens
      AM <- file_name$linReg$AM
      colnames(AM) <- c('sample_id', 'age', 'age_uncert_pos', 'age_uncert_neg')
      plotsynSpeleo_method(ens, AM, dating,depth = NA, 'a) LR')
      plotTrueAge(speleo_new)
      #print('geht nicht')
    }

    #print('1')

    if(!plyr::empty(data.frame(file_name$linInterp$ens))){
      ens <- file_name$linInterp$ens
      AM <- file_name$linInterp$AM
      colnames(AM) <- c('sample_id', 'age', 'age_uncert_pos', 'age_uncert_neg')
      plotsynSpeleo_method(ens, AM, dating, depth=NA, 'b) LI')
      plotTrueAge(speleo_new)
    }

    if(!plyr::empty(data.frame(file_name$copRa$ens))){
      ens <- file_name$copRa$ens
      AM <- file_name$copRa$AM
      colnames(AM) <- c('sample_id', 'age', 'age_uncert_pos', 'age_uncert_neg')
      plotsynSpeleo_method(ens, AM, dating, depth=NA, 'c) copRa')
      plotTrueAge(speleo_new)
      }

    if(!plyr::empty(data.frame(file_name$StalAge$AM))){
      AM <- file_name$StalAge$AM
      colnames(AM) <- c('sample_id', 'age', 'age_uncert_pos', 'age_uncert_neg')
     plotsynSpeleo_method(ens=NA, AM, dating, ens[,2], 'd) StalAge',T)
     plotTrueAge(speleo_new)
    }

    if(!plyr::empty(data.frame(file_name$Bacon$ens))){
      ens <- file_name$Bacon$ens
      ens[,2] <- ens[,2]*10
      AM <- file_name$Bacon$AM
      colnames(AM) <- c('sample_id', 'age', 'age_uncert_pos', 'age_uncert_neg')
      plotsynSpeleo_method(ens, AM, dating, depth = NA, 'e) Bacon')
      plotTrueAge(speleo_new)
      }

    if(!plyr::empty(data.frame(file_name$Bchron$ens))){
      ens <- file_name$Bchron$ens
      AM <- file_name$Bchron$AM
      colnames(AM) <- c('sample_id', 'age', 'age_uncert_pos', 'age_uncert_neg')
      plotsynSpeleo_method(ens, AM, dating, depth=NA, 'f) Bchron')
      plotTrueAge(speleo_new)
    }

    graphics.off()
  }
}



