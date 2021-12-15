prec_eval <- synSpeleo_chrono %>%
  mutate(InOut_lR = if_else(true.age > lin_reg_age -lin_reg_age_uncert_neg & true.age < lin_reg_age +lin_reg_age_uncert_pos,T,F),
         InOut_lI = if_else(true.age > lin_interp_age -lin_interp_age_uncert_neg & true.age < lin_interp_age +lin_interp_age_uncert_pos,T,F),
         InOut_copRa = if_else(true.age > copRa_age -copRa_age_uncert_neg & true.age < copRa_age +copRa_age_uncert_pos,T,F),
         InOut_StalAge = if_else(true.age > StalAge_age -StalAge_age_uncert_neg & true.age < StalAge_age +StalAge_age_uncert_pos,T,F),
         InOut_bchron = if_else(true.age > bchron_age -bchron_age_uncert_neg & true.age < bchron_age +bchron_age_uncert_pos,T,F),
         InOut_bacon = if_else(true.age > bacon_age -bacon_age_uncert_neg & true.age < bacon_age +bacon_age_uncert_pos,T,F)) %>%
  group_by(entity_id) %>%
  mutate(lR = sum(InOut_lR,na.rm=T)/sum(!is.na(InOut_lR)),
         lI = sum(InOut_lI,na.rm=T)/sum(!is.na(InOut_lI)),
         copRa = sum(InOut_copRa,na.rm=T)/sum(!is.na(InOut_copRa)),
         StalAge = sum(InOut_StalAge,na.rm=T)/sum(!is.na(InOut_StalAge)),
         bacon = sum(InOut_bacon,na.rm=T)/sum(!is.na(InOut_bacon)),
         bchron = sum(InOut_bchron,na.rm=T)/sum(!is.na(InOut_bchron))) %>%
  group_byqa_
