library(Hmisc)
setwd("~/Documents/holo&lgm/48-T8/linInterp")
dating_tb <- read.csv('ages.csv', header = T, stringsAsFactors = F)
depth_sample <- read.csv('depths.csv', header = T, stringsAsFactors = F, colClasses = c('numeric', 'numeric'))

setwd("~/Documents/holo&lgm/48-T8")
original_tb <- read.csv('original_chronology.csv', header = T, stringsAsFactors = F)


dating_tb <- as.tibble(dating_tb)
dating_tb <- add_row(dating_tb, dating_id = 0, corr_age = -48, corr_age_uncert = 10, depth_dating = 0, date_type = 'added by me', .before = 1) # active according to paper
dating_tb_new <- dating_tb %>% slice(-5) %>% add_row(dating_id = 777, corr_age = 604, corr_age_uncert = 9, depth_dating = 101, date_type = 'TIMS changed by me') %>% arrange(., corr_age)# corrected by 50 years

age_1 <- (dating_tb_new %>% filter(corr_age <= 10200))$corr_age
depth_dating_1 <- (dating_tb_new %>% filter(corr_age <= 10200))$depth_dating
x_out_1 <- (original_tb %>% filter(interp_age <= 10200))$interp_age

e_1 <- approxExtrap(x = age_1, y = depth_dating_1, xout = x_out_1)
linInterp <- as_tibble(data.frame(depth_eval = unlist(e_1$y), lin_interp_age = unlist(e_1$x)))

age_2 <- (dating_tb_new %>% filter(corr_age <= 24400 & corr_age >= 12700))$corr_age
depth_dating_2 <- (dating_tb_new %>% filter(corr_age <= 24400 & corr_age >= 12700))$depth_dating
x_out_2 <- (original_tb %>% filter(interp_age <= 24400 & interp_age >= 12700))$interp_age

m <- lm(depth_dating_2 ~ age_2)
e_2<- m$coefficients[[1]] + x_out_2*m$coefficients[[2]]

linReg <- data.frame(depth_eval = e_2, lin_interp_age = x_out_2)

lin_48 <- rbind(linInterp, linReg)

sample_id <- data.frame(sample_id = original_tb$sample_id, lin_interp_age = original_tb$interp_age)

ID48 <- left_join(sample_id, lin_48, by = 'lin_interp_age') %>% select(sample_id, depth_eval, lin_interp_age) 
names(ID48) <- c('sample_id', 'depth_sample', 'interp_age')

ID48_new <- ID48 %>% mutate(diff_dft = lead(depth_sample)-depth_sample,
                            diff_age = lead(interp_age)-interp_age) %>%
  mutate(growth_rate = diff_dft/diff_age)

ave_ID48_young <- ID48_new %>% filter(interp_age < 10200) %>% summarise(avg = median(growth_rate, na.rm = T),
                                                                        q1 = quantile(growth_rate, probs = 0.025, na.rm = T),
                                                                        q2 = quantile(growth_rate, probs = 0.975, na.rm = T),
                                                                        s = sd(growth_rate, na.rm = T))
ave_ID48_old <- ID48_new %>% filter(interp_age > 12700) %>% summarise(avg = mean(growth_rate, na.rm = T),
                                                                      q1 = quantile(growth_rate, probs = 0.025, na.rm = T),
                                                                      q2 = quantile(growth_rate, probs = 0.975, na.rm = T))
setwd("~/R/MWE/Bacon_runs/test")

test_core <- read.csv('test_alt.csv', header = T)
test_core_new <- test_core %>% add_row(ID = 'P4',age = -40, error = 6, depth = 5.0,.before = 1) %>% add_row(ID = 'P6', age = 1120, error = 40, depth = 160.0, .before = 4) #%>% mutate(cc = 0)
write.csv(test_core_new, 'test.csv', row.names = F)

write.table(seq(0,400), 'test_depths.txt', row.names = F, col.names = F)
dft <- read.table('test_depths.txt')


setwd("~/R/MWE")
Bacon('test',depths.file = T, hiatus.depths = c(100,330), postbomb = 4)
