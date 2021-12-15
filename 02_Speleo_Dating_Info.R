# Find out data from the SISAL database

library(zoo)
library(Hmisc)
library(pracma)
library(plyr)
library(tidyverse)

dating_tab = read.csv("/obs/proxydata/speleo/SISAL_v2/dating.csv", stringsAsFactors=FALSE)

entities = unique(dating_tab$entity_id)

dating_info = tibble(entity_id = numeric(length(entities)), tot.depth = numeric(length(entities)), 
                     tot.age = numeric(length(entities)), 
                     tot.gr = numeric(length(entities)),
                     nsamp = numeric(length(entities)), 
                     mean.uncert = numeric(length(entities)),
                     max.uncert = numeric(length(entities)), min.uncert = numeric(length(entities)))



for(ii in 1:length(entities)){
  entity = entities[ii]
  tab = dating_tab %>% filter(entity_id == entity)
  tab$depth_dating = as.numeric(tab$depth_dating)
  tab$uncorr_age = as.numeric(tab$uncorr_age)
  
  tab$depth_dating[is.infinite(tab$depth_dating)] = NA
  tab$uncorr_age[is.infinite(tab$uncorr_age)] = NA
  dating_info$entity_id[ii] = entity
  if(all(is.na(tab$depth_dating))){
    dating_info$tot.depth[ii] <- dating_info$tot.age[ii] <- dating_info$nsamp[ii] <- NA
    dating_info$mean.uncert[ii] <- dating_info$max.uncert[ii] <- dating_info$min.uncert[ii] <- NA
    next
  }else{
    dating_info$tot.depth[ii] = max(tab$depth_dating, na.rm = T)
  }
  if(all(is.na(tab$uncorr_age))){
    dating_info$tot.age[ii] = NA
  }else{
    dating_info$tot.age[ii] =  max(as.numeric(tab$corr_age), na.rm = T)-min(as.numeric(tab$corr_age), na.rm = T)
  }
  
  if(!is.na(dating_info$tot.age[ii])&!is.na(dating_info$tot.depth[ii])){
    dating_info$tot.gr[ii] = dating_info$tot.depth[ii]/dating_info$tot.age[ii]
  }else{dating_info$tot.gr[ii] = NA}
  
  dating_info$nsamp[ii] = length(tab$dating_id)
  dating_info$mean.uncert[ii] = mean(c(as.numeric(tab$uncorr_age_uncert_pos), 
                                       as.numeric(tab$uncorr_age_uncert_neg)), na.rm = T)
  dating_info$max.uncert[ii] = max(c(as.numeric(tab$uncorr_age_uncert_pos), 
                                      as.numeric(tab$uncorr_age_uncert_neg)), na.rm = T)
  dating_info$min.uncert[ii] = min(c(as.numeric(tab$uncorr_age_uncert_pos), 
                                      as.numeric(tab$uncorr_age_uncert_neg)), na.rm = T)
}

#Mean Depth
print(paste("Mean depth:",mean(dating_info$tot.depth, na.rm = T)))
print(paste("Mean age:",mean(dating_info$tot.age, na.rm = T)))

print(paste("Mean depth divided by mean age covered: ",  round(mean(dating_info$tot.depth, na.rm = T)/mean(dating_info$tot.age, na.rm = T)*1000, digits = 2), "mm/ky"))
print("This differs from the mean growth rate over all speleothems, which is:")
print(paste("Mean growth rate:", round(mean(dating_info$tot.gr, na.rm = T)*1000, digits = 2), "mm/ky"))
print(paste("Max growth rate:", round(max(dating_info$tot.gr, na.rm = T)*1000, digits = 2), "mm/ky"))
print(paste("Min growth rate:", round(min(dating_info$tot.gr, na.rm = T)*1000, digits = 2), "mm/ky"))

#1) Whats the medium precision in the database?

# Precision has to be given in relation to total age covered. 
# -> if we want this over the 10000 arbitrary age units, the prec should be over the tot length multiplied by the 10000 a.u.

print(paste("In relation to the 10 000 a.u in age, this is a uncertainty of ", round(mean(dating_info$mean.uncert, na.rm = T)/mean(dating_info$tot.age, na.rm = T)*10000, digits = 2)))


#2) Whats the number of measurements per depth?

print(paste("On average we have ", round(mean(dating_info$nsamp, na.rm = T), digits = 2), "age measurements per speleothem"))
print(paste("This is in relation to the synth speleo:", round(mean(dating_info$nsamp/dating_info$tot.age*10000, na.rm = T), digits = 2), "age measurements per speleothem over 10000y"))
print(paste("The highest resolved speleo has :", round(max(dating_info$nsamp/dating_info$tot.age*10000, na.rm = T), digits = 2), "age measurements per speleothem over 10000y"))
print(paste("The lowest resolved speleo has :", round(min(dating_info$nsamp/dating_info$tot.age*10000, na.rm = T), digits = 2), "age measurements per speleothem over 10000y"))

