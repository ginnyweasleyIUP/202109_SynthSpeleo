### Plot anonymous pics 

runFile <- r %>% filter(!(entity_id %in% c(112,130,131)))

plot_ano <- function(chrono, age, depth, hiatus_tb, entity_id, linReg = F, linInterp = F, copRa = F, StalAge = F, Bacon = F, Bchron = F){
  matplot(x = chrono$age, y= depth$depth_sample, col = 'black', lty = 1, type = 'l', lwd = 1, 
          xlim = c(min(chrono$age, na.rm = T), max(chrono$age, na.rm = T)),
          ylim = c(max(depth$depth_sample, na.rm = T),min(depth$depth_sample, na.rm = T)), xlab = 'Age [yrs BP]', ylab = 'Depth from top [mm]')
  lines(x=chrono$age + chrono$uncert_pos,y=depth$depth_sample, lty = 2, col = 'red')
  lines(x=chrono$age - chrono$uncert_neg,y = depth$depth_sample, lty = 2, col = 'red')
  points(x = age$corr_age, y=age$depth_dating, lty = 2, col = 'orange', pch = 4)
  arrows(age$corr_age-age$corr_age_uncert, age$depth_dating, age$corr_age+age$corr_age_uncert, age$depth_dating, length=0.05, angle=90, code=3, col = 'orange')
  legend("topright", legend = c('date','Median age', 'CI'), pch = c(4, NA, NA), lty = c(NA,1,2), title = 'AM info', col = c('orange', 'black', 'red'), cex = 0.7)
  if (!empty(data.frame(hiatus_tb))) {
    abline(h = hiatus_tb$depth_sample, col = 'grey', lty = 2)
  }
  if(linReg){mtext(line = 2, text = paste(entity_id,' -1'))}
  if(linInterp){mtext(line = 2, text = paste(entity_id,' -2'))}
  if(copRa){mtext(line = 2, text = paste(entity_id,' -3'))}
  if(StalAge){mtext(line = 2, text = paste(entity_id,' -4'))}
  if(Bacon){mtext(line = 2, text = paste(entity_id,' -5'))}
  if(Bchron){mtext(line = 2, text = paste(entity_id,' -6'))}
}

r_test <- runFile %>% filter(entity_id < 10)
r_new <- read_chronology(runFile, file_path = "/home/ariana/Documents/Expert.Ranking", entity. = entity)

read_chronology<- function(runFile,file_path,entity.){
  
  for(i in runFile$entity_id) {
    graphics.off()
    print(i)
    y <- runFile %>% filter(entity_id == i) 
    entity_name <- (entity. %>% filter(entity_id == i))$entity_name
    file_name <- paste(i,"-",entity_name, sep = '')
    dir.create(file.path(file_path,i))
    
    depths <- read.csv(file.path(y$working_directory, file_name, "linInterp/depths.csv"), header = T, stringsAsFactors = F)
    ages <- read.csv(file.path(y$working_directory, file_name, "linReg/ages.csv"), header = T, stringsAsFactors = F)
    hiatus_tb <- read.csv(file.path(y$working_directory, file_name, "hiatus.csv"), header = T, stringsAsFactors = F)
    
    if(y$linReg){
      setwd(file.path(y$working_directory, file_name, '/linReg'))
      if(file.exists('linReg_chronology.csv')){
        chrono <- read.csv('linReg_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% distinct()
        names(chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        if(!all(is.na(chrono$age))){
        pdf(file.path(file_path,i," -1.pdf"), 6,4)
        plot_ano(chrono, ages, depths, hiatus_tb, i, linReg = T)
        dev.off()} else {runFile[i,7] <- F}
        print('ja')
      } else {runFile[i,7] <- F
        print('nein')}
    }
    
    if(y$linInterp){
      setwd(file.path(y$working_directory, file_name, '/linInterp'))
      if(file.exists('linInt_chronology.csv')){
        chrono <- read.csv('linInt_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% distinct()
        names(chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        pdf(file.path(file_path,i," -2.pdf"), 6,4)
        plot_ano(chrono, ages, depths, hiatus_tb, i, linInterp = T)
        dev.off()
        print('ja')
      } else {runFile[i,5] <- F
      print('nein')}
    }
    
    if(y$copRa){
      setwd(file.path(y$working_directory, file_name, '/copRa'))
      if(file.exists('copRa_chronology.csv')){
        chrono <- read.csv('copRa_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% distinct()
        names(chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        pdf(file.path(file_path,i," -3.pdf"), 6,4)
        plot_ano(chrono, ages, depths, hiatus_tb, i, copRa = T)
        dev.off()
        print('ja')
      } else {runFile[i,6] <- F
      print('nein')}
    }
    
    if(y$Bchron){
      setwd(file.path(y$working_directory, file_name, '/Bchron'))
      if(file.exists('bchron_chronology.csv')){
        chrono <- read.csv('bchron_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% distinct()
        names(chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        pdf(file.path(file_path,i," -6.pdf"), 6,4)
        plot_ano(chrono, ages, depths, hiatus_tb, i, Bchron = T)
        dev.off()
        print('ja')
      } else {runFile[i,3] <- F
      print('nein')}
    }
    
    if(y$StalAge){
      setwd(file.path(y$working_directory, file_name, '/StalAge'))
      if(file.exists('StalAge_chronology.csv')){
        chrono <- read.csv('StalAge_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% distinct()
        names(chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        pdf(file.path(file_path,i," -4.pdf"), 6,4)
        plot_ano(chrono, ages, depths, hiatus_tb, i, StalAge = T)
        dev.off()
        print('ja')
      } else {runFile[i,4] <- F
      print('nein')}
    }
    
    if(y$Bacon){
      setwd(file.path(y$working_directory, file_name, '/Bacon_runs'))
      if(file.exists('bacon_chronology.csv')){
        chrono <- read.csv('bacon_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% distinct()
        names(chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        pdf(file.path(file_path,i," -5.pdf"), 6,4)
        plot_ano(chrono, ages, depths, hiatus_tb, i, Bacon = T)
        dev.off()
        print('ja')
      } else {runFile[i,2] <- F
      print('nein')}
    }
    
  }
  return(runFile)
}



if(T){
  
  if(F){
    chrono <- read.csv('linInt_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) %>% distinct()
    names(chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
    pdf(file.path(file_path,i," -2.pdf"), 6,4)
    plot_ano(chrono, ages, depths, hiatus_tb, i, linInterp = T)
    dev.off()
  } else {runFile %>% mutate(linInterp = if_else(entity_id == 1, F, linInterp))}
}

