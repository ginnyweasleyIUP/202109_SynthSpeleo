
data <- dating_tb_new %>% filter(entity_id == 319 & date_used == 'yes') %>% select(dating_id, corr_age, corr_age_uncert, depth_dating, date_type)
depth_eval <- sample_tb %>% filter(entity_id == 319 & is.na(hiatus)) %>% select(sample_id, depth_sample)
hiatus <- sample_tb %>% filter(entity_id == 319 & hiatus == 'H') %>% select(sample_id, depth_sample)
method <- 'linear'

get_lin_interp <- function(data, depth_eval, method, hiatus) {
  
  library(plyr)
  library(rlist)
  library(Hmisc)
  
  if (!empty(data.frame(hiatus))) {
    data <- add_hiatus(data, hiatus)
  }
  
  age <- unlist(data[,2])
  depth_dating <- unlist(data[,4])
  x_out <- unlist(depth_eval)
  
  e <- approxExtrap(x = depth_dating, y = age, xout = x_out)
  #e2 <- approx(x = depth_dating, y = age, xout = x_out, method = 'linear')
  
  #linInterp2 <- data.frame(depth_eval = e2$x, lin_interp_age = e2$y)
  linInterp <- as_tibble(data.frame(depth_eval = unlist(e$x), lin_interp_age = unlist(e$y)))
  
  if (!empty(data.frame(hiatus))) {
    linInterp <- linInterp %>% rowwise() %>% mutate(lin_interp_age = if_else(depth_eval %in% hiatus$depth_sample, NA_real_, lin_interp_age))
  }
  
  return(linInterp)
}


data <- dating_tb %>% filter(entity_id == 319 & date_used == 'yes' & date_type != 'Event; hiatus') %>% select(dating_id, corr_age, corr_age_uncert, depth_dating, date_type) %>% 
  arrange(., depth_dating)
hiatus <- sample_tb %>% filter(entity_id == 319 & hiatus =='H') %>% select(sample_id, depth_sample) %>% arrange(., depth_sample)

add_hiatus <- function(data, hiatus) {
  
  age <- unlist(data[,2])
  depth_dating <- unlist(data[,4])
  #uncert <- unlist(data[,3])
  x_out <- unlist(hiatus$depth_sample)
  
  e <- approx(x = depth_dating, y = age, xout = x_out, method = 'linear')
  h <- data.frame(dating_id = hiatus$sample_id,corr_age = e$y, corr_age_uncert = rep(NA, length(hiatus$depth_sample)),depth_dating = e$x, date_type = rep('Hiatus',length(hiatus$depth_sample)))
  new <- rbind(data, h)
  new <- new %>% arrange(., corr_age) %>% mutate(corr_age_uncert = if_else(is.na(corr_age_uncert), (lead(corr_age_uncert)+lag(corr_age_uncert))/2, as.double(corr_age_uncert)))
  
  return(new) 
}

