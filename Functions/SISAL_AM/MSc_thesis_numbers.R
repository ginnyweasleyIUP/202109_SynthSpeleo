res <- sisalv2interim$dating %>% filter(date_used == 'yes' & date_type != 'Event; hiatus') %>% group_by(entity_id) %>% 
  summarize(min = min(depth_dating, rm.na = T), max = max(depth_dating, rm.na = T), n = n(), r = n/(max-min), range.t = )

min(res$r, na.rm = T)
max(res$r, na.rm = T)
max(res$n, na.rm=T)
mean(res$n, na.rm=T)



sites <- sisalv2interim$entity %>% group_by(site_id) %>% summarise(n = n())


ave(sites$n, na.rm = T)
max(sites$n, na.rm = T)

reversals <- sisalv2interim$dating %>% filter(date_used == 'yes' & date_type != 'Event; hiatus') %>% group_by(entity_id) %>% arrange(., depth_dating, .by_group=T) %>% 
  mutate_at(vars(corr_age_uncert_pos, corr_age_uncert_neg),as.numeric) %>%
  mutate(diff = if_else(lead(corr_age)-corr_age < 0, T, F), prec = (corr_age_uncert_pos+corr_age_uncert_neg)/(2*(corr_age+1))) %>% 
  summarise(n=sum(diff, na.rm=T), prec_avg = mean(prec, na.rm=T))

max(reversals$n)
mean(reversals$n)
mean(reversals$prec_avg, na.rm=T)*100



hiatus <- left_join(sisalv2interim$sample, sisalv2interim$hiatus, by = 'sample_id') %>% mutate(h = if_else(hiatus == 'H', T, F)) %>% group_by(entity_id) %>% summarise(h = sum(h, na.rm = T))

max(hiatus$h, na.rm=T)
mean(hiatus$h, na.rm = T)
hiatus %>% filter(h > 0) %>% count()
unique(hiatus$h)


s <- SISAL_eval_new_new %>% filter(stalage == 'dropped')

                      

