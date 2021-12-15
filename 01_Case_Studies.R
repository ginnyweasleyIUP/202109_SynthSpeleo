# Case Study Analysis

# Case 1: no growthrate change gr = 1, number of hiatus = 1
# Case 2: 1 average growth rate change: gr_1 = 1, gr_2 in (1.5,3), length(period_1) = length(period_2)
# Case 3: 1 average growth rate change: gr_1 = 1, gr_2 in (1.5,3)
#         50% of time: 0.2*tot_length < length(period_1) < 0.4*tot_length
#         50% of time: 0.2*tot_length < length(period_2) < 0.4*tot_length
#         lengths are uniformly distributed
# Case 4: like Case 2 but with hiatus of 5% fixed length at uniformly random position btw 0.1 and 0.9
# Case 5: like Case 3 but with hiatus of 5% fixed length at uniformly random position btw 0.1 and 0.9
# Case 6: gr_i = 1, gr_j in (1.5,3), gr_k = 10
#         33% of times extreme at the beginning
#         33% of times extreme at the end
#         33% of times extreme in the middle
# Case 7: like Case 6, bus gr_k = hiatus


# Count:  e.g. experiment_1_20_57 --> Case 1 with 20 samples with realisation 57
#              + txt file with growth rate and periods and hiatus info
#              each folder contains all Age-Models so this does not need to go into the name! 

source(file = 'Functions/SISAL_AM/run_functions.R')
source(file = 'Functions/SISAL_AM/SISAL_AM_add_functions.R')
source(file = 'Functions/SISAL_AM/SISAL_AM_plot_functions.R')
source(file = "Functions/SISAL_AM/StalAge_1_0_mine.R")
source(file = 'Functions/SISAL_AM/New_methods.R')
source(file = 'Functions/SISAL_AM/lin_interp_reversal.R')
source(file = 'syntheticspeleo-master/R/functions.R')
source(file = 'syntheticspeleo-master/R/validation_fct.R')

library(zoo)
library(Hmisc)
library(pracma)
library(plyr)
library(tidyverse)

# synth speleo fixed values
p = 60
t.min = 0
t.max = 10000
t.n = 4 # what is this?
d.min = 0
d.max = 1000
n = 10 #number of realisations
file_path1 = '/home/ldap-server/ginnyweasley/07_R_Code/202109_SynthSpeleo/Case_Study'
file_path2='/home/ldap-server/ginnyweasley/07_R_Code/202109_SynthSpeleo/Case_Study/summary'
t.res=2 #resolution of the true age
prec = 60



quanAna <- function(p.min = 0.1, p.max=1,p.n=5,
                    t.min=0,t.max=10000,
                    t.res=2,
                    d.min=0,d.max=1000,
                    file_path1 ='test/', 
                    file_path2 = 'test/summary', 
                    n=3, prec = 150){
  
  number.runs = 1
  chrono_final <- tibble()
  
  for(case in 1:7){
    for(nsamp in c(30,20,15,12,6,3)){
      for(n.rep in 1:n){
        
        print(paste0("case: ",case,", nsamp: ",nsamp,", n: ",n.rep,". We're in run ", number.runs, " of ",7*6*n))
        
        # 1) create synthetic speleothem
        speleo<-NA
        speleo_new <- NA
        speleo<-generate.speleo(t.top=t.min,t.bottom=t.max,z.top = d.min,z.bottom=d.max,t.res=t.res)
        
        
        if(case == 1){
          tlist2 = c(0,4500,5500,10000)
          rates2 = c(1,0,1)
        }else if(case == 2){
          tlist2 = c(0,5000,10000)
          gr2 = round(1.5+runif(1)*1.5, digits = 2)
          rates2 = c(1,gr2)
          
        }else if(case == 3){
          #1) decide on short at beginning or at end:
          begin.or.end = round(runif(1))
          #ratio between short and long should be 70:30 and 60:40
          if(begin.or.end == 1){
            length.begin = (0.3+round(runif(1)*0.1, digits = 2))*(t.max-t.min)
          }else{
            length.begin = (0.6+round(runif(1)*0.1, digits = 2))*(t.max-t.min)
          }
          
          tlist2 = c(0,length.begin,t.max)
          gr2 = round(1.5+runif(1)*1.5, digits = 2)
          rates2 = c(1,gr2)
        }else if(case == 4){
          #case 2 + 1 hiatus
          tlist2 = c(0,5000,10000)
          gr2 = round(1.5+runif(1)*1.5, digits = 2)
          h.begin = round(0.1+0.8*runif(1), digits = 2)*(t.max-t.min)
          tlist2 = sort(c(tlist2,h.begin, h.begin+1000))
          if(h.begin>4000 & h.begin<5000){tlist2 = sort(c(t.min,t.max,h.begin, h.begin+1000))}
          rates2 = numeric(length(tlist2)-1)
          for(ii in 1:length(rates2)){
            if(tlist2[ii+1]<=5000){rates2[ii] = 1}
            if(tlist2[ii+1]>5000){rates2[ii] = gr2}
            if(tlist2[ii] == h.begin){rates2[ii] = 0}
          }
        }else if(case == 5){
          #1) decide on short at beginning or at end:
          begin.or.end = round(runif(1))
          #ratio between short and long should be 70:30 and 60:40
          if(begin.or.end == 1){
            length.begin = (0.3+round(runif(1)*0.1, digits = 2))*(t.max-t.min)
          }else{
            length.begin = (0.6+round(runif(1)*0.1, digits = 2))*(t.max-t.min)
          }
          
          tlist2 = c(0,length.begin,t.max)
          gr2 = round(1.5+runif(1)*1.5, digits = 2)
          h.begin = round(0.1+0.8*runif(1), digits = 2)*(t.max-t.min)
          tlist2 = sort(c(tlist2,h.begin, h.begin+1000))
          if(h.begin>length.begin-1000 & h.begin<length.begin){tlist2 = sort(c(t.min,t.max,h.begin, h.begin+1000))}
          rates2 = numeric(length(tlist2)-1)
          for(ii in 1:length(rates2)){
            if(tlist2[ii+1]<=length.begin){rates2[ii] = 1}
            if(tlist2[ii+1]>length.begin){rates2[ii] = gr2}
            if(tlist2[ii] == h.begin){rates2[ii] = 0}
          }
          
        }else if(case == 6){
          gr2 = 2
          begin.or.end = runif(1)
          average.rate = sample(c(1,gr2))
          tlist2 = c(0,3333,6666,10000)
          if(begin.or.end <1/3){
            rates2 = c(1/10,average.rate[1],average.rate[2])*10
          }else if(begin.or.end >= 1/3 & begin.or.end<2/3){
            rates2 = c(average.rate[1],1/10,average.rate[2])*10
          }else{
            rates2 = c(average.rate[1],average.rate[2],1/10)*10
          }
          
        }else if(case == 7){
          average.rate = sample(c(1,2))
          tlist2 = c(0,3333,6666,10000)
          rates2 = c(average.rate[1],0,average.rate[2])*10
        }
        
        print(tlist2)
        print(rates2)
        speleo_new <- NA
        speleo_new<-growthfcn.changinggamma.new(speleo, tlist = tlist2, rates = rates2)
        
        speleo_new$datingparams$tdofs<-10
        speleo_new$datingparams$precision<-prec
        speleo_new$datingparams$topage<-10
        speleo_new$datingparams$nsample<-nsamp
        
        speleo_new<-do.dating(speleo_new)
        
        plotgrowthfcn(speleo_new)
        plotdating(speleo_new$Dtable,add=TRUE)
        
       
        pdf(file.path(file_path2,paste('grwthfct-',case,'-',nsamp,'-',n.rep,'.pdf',sep = '')),6,4)
        plotgrowthfcn(speleo_new)
        plotdating(speleo_new$Dtable,add=TRUE)
        graphics.off()
        
        hiatus <- unique(speleo_new$growthfcn$expected.z[which(diff(speleo_new$growthfcn$expected.z)==0)])
        
        
        if(case %in% c(2,3,6)){hiatus <-tibble(sample_id = numeric(), depth_sample = numeric())} else {hiatus <- tibble(sample_id = -n.rep*10-seq(1:length(hiatus)), depth_sample = hiatus)}
        
        speleo_new<-sample.proxy(speleo = speleo_new,nsamp = 700,z.steps = 1,stype = "sinusoid")
        speleo_new$hiatus <- hiatus
        #dind<-check.Dtable(speleo_class1,add=TRUE)
        save(speleo_new,file = file.path(file_path2, paste('grwthfct-',case,'-',nsamp,'-',n.rep,'.RData',sep='')))
        run_SISAL_chrono(speleo_new, file_path2,paste(case,'-',nsamp,'-',n.rep,sep = '') , hiatus)
        
        file_name <- paste(case,'-',nsamp,'-',n.rep,sep='')
        
        if(case %in% c(2,3,6)){
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
        
        chrono <- merge_synSpeleo_qa(file_name, file_path2, sample) %>% arrange(., depth_sample)
        
        ind.samples<-sapply(chrono$depth_sample,function(x,growthfcn){which.min(abs(growthfcn$growth.real-x))},growthfcn=speleo_new$growthfcn)
        true.age<-speleo_new$growthfcn$true.age[ind.samples]
        #true.age <- grf$true_age[sapply(chrono$depth_sample,function(x) which.min(abs(x-grf$expected.z)))]
        
        speleo_new$recorded.age <- true.age
        
        
        
        chrono_new <- chrono %>% mutate(true.age = true.age, precision = rep(prec, length(chrono$depth_sample)),
                                        nsamp_nr = rep(nsamp, length(chrono$depth_sample)),
                                        case_nr = rep(case, length(chrono$depth_sample)), iteration = rep(n.rep, length(chrono$depth_sample))) %>%
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
        if(case %in% c(1,4,5,7)){
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
        save(chrono_new, file = file.path(file_path2, paste('chrono-',case,'-',nsamp,'-',n.rep,'.RData',sep='')))        
        
        rates_times = list(tlist = tlist2, rates = rates2)
        save(rates_times, file = paste0(file_path2,"/",case,"-",nsamp,"-",n.rep,"/rates_times.RData"))
        
        number.runs = number.runs + 1
      }
    }
  }
  return(chrono_final)
}

synSpeleo_chrono = quanAna(file_path1 = '/home/ldap-server/ginnyweasley/07_R_Code/202109_SynthSpeleo/Case_Study',
                           file_path2='/home/ldap-server/ginnyweasley/07_R_Code/202109_SynthSpeleo/Case_Study/summary3', n = 10)

save(synSpeleo_chrono, file = '/home/ldap-server/ginnyweasley/07_R_Code/202109_SynthSpeleo/Case_Study/synSpeleo_chrono_2.RData')


#Ok so all cases are run with 3 realisations. 
# Let's try Carlas plotting for it here:


############## HERE #########################################################################

#72=n*(p.n+1)*(lenght(h.length)+1)
names <- rep(NA_real_,189)
counter <- 0
for(case in 1:7){
  for(nsamp in c(30,20,15,12,10,8,6,4,3)){
    for(n.rep in 1:3){
      counter <- counter+1
      name <- paste('/home/ldap-server/ginnyweasley/07_R_Code/202109_SynthSpeleo/Case_Study/summary/chrono-',case,'-',nsamp,'-',n.rep,'.RData',sep='')
      names[counter] <- name
    }
  }
}

case_nsamp <- tibble(entity_id = character(),case_nr = numeric(),prec = numeric(), iteration = numeric(), nsamp_nr = numeric(),lr = numeric(),li=numeric(),copRa=numeric(),bacon=numeric(),stalage=numeric())
for(name in names){
  
  #load(file.path('/stacycode/ariana/synSpeleo/summary',name))
  load(file.path(name))
  
  chrono_sum <- chrono_new %>% summarise(entity_id = unique(entity_id), 
                                         case_nr = 0,#case_nr = unique(case_nr),
                                         nsamp_nr = 0,#nsamp_nr = unique(nsamp_nr), 
                                         iteration = unique(iteration), prec = unique(precision),
                                         lr = sum(true.age > lin_reg_age -lin_reg_age_uncert_neg & true.age < lin_reg_age + lin_reg_age_uncert_pos, na.rm = T)/(sum(!is.na(lin_reg_age))),
                                         li = sum(true.age > lin_interp_age -lin_interp_age_uncert_neg & true.age < lin_interp_age + lin_interp_age_uncert_pos, na.rm = T)/(sum(!is.na(lin_interp_age))),
                                         copRa = sum(true.age > copRa_age - copRa_age_uncert_neg & true.age < copRa_age + copRa_age_uncert_pos, na.rm = T)/(sum(!is.na(copRa_age))),
                                         bacon = sum(true.age > bacon_age -bacon_age_uncert_neg & true.age < bacon_age + bacon_age_uncert_pos, na.rm = T)/(sum(!is.na(bacon_age))),
                                         stalage = sum(true.age > StalAge_age -StalAge_age_uncert_neg & true.age < StalAge_age + StalAge_age_uncert_pos, na.rm = T)/(sum(!is.na(StalAge_age))),
                                         bchron = sum(true.age > bchron_age - bchron_age_uncert_neg & true.age < bchron_age + bchron_age_uncert_pos, na.rm = T)/(sum(!is.na(bchron_age))))
  
  case_nsamp <- bind_rows(case_nsamp, chrono_sum)
}

case_nsamp$case_nr = c(rep(1,27),rep(2,27),rep(3,27),rep(4,27),rep(5,27),rep(6,27),rep(7,27))
case_nsamp$nsamp_nr = c(rep(c(30,30,30,20,20,20,15,15,15,12,12,12,10,10,10,8,8,8,6,6,6,4,4,4,3,3,3),7))

samp <- (case_nsamp %>%  distinct(nsamp_nr))$nsamp_nr

mtx_i <- list(matrix(),matrix(),matrix(),matrix(),matrix(), matrix())
#am <- c('lr','li','copRa','bacon','stalage', 'bchron')
for(f in 1:6){
  mtx <- matrix(NA_real_,90,70)
  j <- 0
  l <- 0
  #l goes over samples
  for (l in 1:9){
    i <- samp[l]
    #n <- 1
    #h sind die cases
    for(h in 1:7){
      m <- case_nsamp %>% filter(case_nr == h & nsamp_nr == i) %>% arrange(.,iteration) %>% pull(f+5)
      mtx_r <- matrix(m,nrow=10,ncol=10,byrow = T)
      mtx[((l-1)*10)+1:10,((h-1)*10)+1:10] <- mtx_r
    }
  }
  
  mtx_i[[f]] <- mtx
}

case_nsamp_new <- case_nsamp %>% replace(is.na(.),0) %>% group_by(nsamp_nr,case_nr) %>% summarise(lr.prec = mean(lr, na.rm = T),
                                                                                                  li.prec = mean(li, na.rm = T),
                                                                                                  copRa.prec = mean(copRa, na.rm = T),
                                                                                                  stalage.prec = mean(stalage, na.rm = T),
                                                                                                  bacon.prec = mean(bacon, na.rm = T),
                                                                                                  bchron.prec = mean(bchron, na.rm = T))

mean(case_nsamp_new$lr.prec)
mean(case_nsamp_new$li.prec)
mean(case_nsamp_new$copRa.prec)
mean(case_nsamp_new$stalage.prec)
mean(case_nsamp_new$bacon.prec)
mean(case_nsamp_new$bchron.prec)

col_ipcc_sequential_blue <- function(n) {
  colorRampPalette(c(rgb(237,248,251,maxColorValue = 255),
                     rgb(179,205,227,maxColorValue = 255),
                     rgb(140,150,198,maxColorValue = 255),
                     rgb(136,86,167,maxColorValue = 255),
                     rgb(129,15,124,maxColorValue = 255))[c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,4,5)],space="rgb")(n)
}


# mtx_lr <- mtx
# mtx_li <- mtx
# mtx_bacon <- mtx
# mtx_stalage <- mtx
# mtx_copRa <- mtx
#
#
# mtx_i <- list(mtx_lr,mtx_li,mtx_copRa,mtx_stalage,mtx_bacon)

title <- c('a) Linear Regression','b) Linear Interpolation','c) copRa','d) StalAge','e) Bacon', 'f) Bchron')

layout.m <- matrix(c(1,2,3,4,5,6,7,8),nrow = 4,ncol = 2,byrow = TRUE)

if(T){
  pdf('/home/ldap-server/ginnyweasley/07_R_Code/202109_SynthSpeleo/Plots/FirstCaseStudy.pdf',8,12)
  layout(mat = layout.m,heights = c(0.3,0.3,0.3,0.1))
  
  for(i in 1:6){
    par(mar=c(5,4,2,0.5))
    image(1:90,1:70,mtx_i[[i]], col = col_ipcc_sequential_blue(230), xlab = '', ylab = '', xaxt = "n", yaxt="n", zlim=c(0,1))
    image(1:90,1:70,array(is.na(mtx_i[[i]]),dim=c(90,70)),add=T,col="black",zlim=c(0.5,1.5))
    #image(1:60,1:40,array(is.na(mtx_i[[i]]),dim=c(60,40)),add=T,col="darkgrey",zlim=c(0.5,1.5))
    mtext(text='Case Number', side = 2, line = 2, cex = 0.8)
    mtext(text='Number of age samples', side = 1, line = 3, cex = 0.8)
    axis(side = 1,at = seq(5.5,85.5,length.out = 9), labels = samp,tick = F)
    axis(side = 1,at = seq(10.5,80.5,10), labels = rep('',8),tick = T, lwd.ticks = 2)
    abline(v=10*(1:9)+0.5, lwd = 2)
    axis(side = 2,at = seq(5.5,64.5,length.out = 7), labels = seq(1,7,1),tick = F)
    axis(side = 2,at = seq(10.5,60.5,10), labels = rep('',6),tick = T, lwd.ticks = 2)
    box(lwd = 2)
    abline(h=10*(1:6)+0.5, lwd=2)
    title(title[i],adj = 0 )
  }
  
  par(new = TRUE)
  par(fig = c(0, 1, 0, 0.1))
  fields::image.plot(1:60,1:40,mtx_i[[5]], col = col_ipcc_sequential_blue(230), zlim = c(0,1),
                     legend.only=T,horizontal = T,smallplot = c(0.15,0.85,0.7,1))
  mtext(text = 'Coverage Frequency', line = 2, side = 1)
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
    m <- (prec_hiatus %>% filter(case_nr == h & prec == i) %>% arrange(.,iteration) %>% select(copRa))$copRa
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

