setwd("~/Documents/CI2019/63-MAW-6/")

o <- read.csv('original_chronology.csv', header = T)
d <- dating_tb %>% filter(entity_id == 63)
u <- read.csv('used_dates.csv', header = T)
nu <- read.csv('not_used_dates.csv', header = T)
h <- read.csv('hiatus.csv', header = T)

used <- u 
used[,c(2,3,4)] <- used[,c(2,3,4)]/1000

setwd("~/Documents/CI2019")
pdf('Abb1.pdf', paper = 'a4')
plot(used[,2], used[,1], ylim = c(210,0), xlim = c(6,18),pch = 4, col = 'black', xlab = '', ylab = '', lwd = 2, axes = FALSE, cex=1.7)
arrows(used[,2]-used[,3], used[,1], used[,2]+used[,3], used[,1], length=0.05, angle=90, code=3, col = 'black', cex= 1.7)
points(nu$corr_age/1000, nu$depth_dating, pch = 25, col = 'red', cex = 1.7)
lines(o$interp_age/1000, o$depth_sample, col = 'black', lwd = 1.5)
lines(o$interp_age/1000 + o$interp_age_uncert_pos/1000, o$depth_sample, col = 'cornflowerblue', lty = 2, lwd = 1.5)
lines(o$interp_age/1000 - o$interp_age_uncert_neg/1000, o$depth_sample, col = 'cornflowerblue', lty = 2, lwd = 1.5)
axis(side = 1,  at = seq(6,18,by=2), labels = seq(6,18,by=2), lwd = 1, lwd.ticks = 1, cex.axis = 1)#, padj = 0.5)
rug(x =  7:18, ticksize = -0.02, side = 1, lwd = 1.2)
mtext(side = 1, text = 'Median age [kyrs BP]', line = 2, cex  =1.3)
segments(x1 = u$corr_age/1000, y0 = u$depth_dating, x0 = 6, col = 'black', lty = 3, lwd = 0.7)
segments(x0 = 6, x1 = 18,y0=h$depth_sample, col = 'orange', lty = 2, lwd = 0.7)
legend("topright",legend = c("Age - used", 'Age - not used'), 
       pch = c(4,25),col = c("black", "red") ,bty='n', ncol = 1, title = "Age-depth model date types", cex = 1, pt.cex = 1.5)
legend("bottomleft",legend = c("Original age-depth model", '90% confidence intervals', 'Hiatus', 'Dating depths'), lwd = 2, lty = c(1,2, 2, 3),col = c("black", "cornflowerblue", 'orange', 'black') ,
       box.lwd = 0,box.col = 'white',box.lty = 0 ,ncol = 1, title = "Additional age-depth model information", cex = 1, pt.cex = 1.5, bg = 'white')
dev.off()

pdf('~/Documents/CI2019/Abb1_3_new_new.pdf', 22, 20)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow = T), widths = c(2,4), heights = 6)
#par(oma = c(1,2,0,2), xpd = T)
par(mar = c(4,3,0,0), oma = c(3,4,0,1), xpd=NA )
plot.new()
plot.window(ylim = c(210,0), xlim = c(0,3.55))
points(rep(2, length(used$depth_dating)), used$depth_dating, pch = 4, col = 'black', cex = 6, lwd = 6)
points(rep(2, length(nu$depth_dating)), nu$depth_dating, pch = 25, col = 'red', cex = 6, lwd = 6)
segments(x0 = 2, x1 = 4, y0 = u$depth_dating, col = 'black', lty=3, lwd=4)
segments(x0=1, x1=3,y0=h$depth_sample, col = 'orange', lty=2, lwd=6)
segments(x0=3, x1=4,y0=h$depth_sample, col = 'orange', lty=2, lwd=4)
axis(side = 2, lwd.ticks = 3, cex.axis = 4, labels = seq(210,0, by = -30), at= seq(210,0, by = -30), lwd = 5)#, padj = -0.5)
rug(x = seq(210,0,-10), ticksize = -0.02, side = 2, lwd = 4)
mtext(side = 2, text = 'Depth from top [mm]', line = 1, cex  =4, outer = T)

plot(used[,2], used[,1], ylim = c(210,0), xlim = c(6,16),pch = 4, col = 'black', xlab = '', ylab = '', lwd = 5, axes = FALSE, cex=6)
arrows(used[,2]-2*used[,3], used[,1], used[,2]+2*used[,3], used[,1], length=0.05, angle=90, code=3, col = 'black', cex= 5)
points(nu$corr_age/1000, nu$depth_dating, pch = 25, col = 'red', cex = 6, lwd = 6)
lines(o$interp_age/1000, o$depth_sample, col = 'black', lwd = 6, cex = 6)
lines(o$interp_age/1000 + o$interp_age_uncert_pos/1000, o$depth_sample, col = 'cornflowerblue', lty = 2, lwd = 6)
lines(o$interp_age/1000 - o$interp_age_uncert_neg/1000, o$depth_sample, col = 'cornflowerblue', lty = 2, lwd = 6)
segments(x1 = u$corr_age/1000, y0 = u$depth_dating, x0 = 6, col = 'black', lty = 3, lwd = 4)
segments(x0 = 6, x1 = 18,y0=h$depth_sample, col = 'orange', lty = 2, lwd = 6)
legend("topright",legend = c("Age - used", 'Age - not used'), 
       pch = c(4,25),col = c("black", "red") ,bty='n', ncol = 1, title = "Age-depth model date types", cex = 3, pt.cex = 5)
legend("bottomleft",legend = c("Original age-depth model", '90% confidence intervals', 'Hiatus', 'Dating depths','Tractable reversal'), lwd = 3, lty = c(1,2, 2, 3, NA),col = c("black", "cornflowerblue", 'orange', 'black', rgb(0.004, 0.67, 0.03)) ,
       box.lwd = 0,box.col = 'white',box.lty = 0 ,ncol = 1, title = "Additional age-depth model information", cex = 3, pt.cex = 4, bg = 'white', pch = c(NA,NA,NA,NA, 1))
axis(side = 1,  at = seq(6,16,by=2), labels = seq(6,16,by=2), lwd = 5, lwd.ticks = 4, cex.axis = 4, padj = 1)
rug(x =  7:18, ticksize = -0.02, side = 1, lwd = 4)
mtext(side = 1, text = 'Median age [kyrs BP]', line = 1.5, cex  =4, outer = T, adj = 0.75)

dev.off()

plotAMC("~/Documents/CI2019", '63-MAW-6', FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)
plotAMC("~/Documents", 'Carole', FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)
working_directory <- setwd("~/Documents/CI2019")
file_name <- '63-MAW-6'
plotAMC <- function(working_directory, file_name, copra, b, bc, stalage , linreg, lininterp, copRa){
  
  setwd(file.path(working_directory,file_name))
  proxy <- read.csv('proxy_data.csv', header = T)
  hiatus <- read.csv('hiatus.csv', header = T)
  original_chrono <- read.csv('original_chronology.csv', header = T)
  #used_age <- read.csv('used_dates.csv', header = T)
  unknown_age <- read.csv('not_used_dates.csv', header = T)
  #used_age <- read.csv('hiatus_dates_interp.csv', header = T)

  depth_sample <- proxy$depth_sample
  #original_age <- original_chrono$interp_age
  
  if(length(hiatus) != 0) {
    used_age <- read.csv('hiatus_dates_interp.csv', header = T)
    #print(hiatus_age)
  } else {
    used_age <- read.csv('used_dates.csv', header = T)
  }
  if(all(is.na(depth_sample))){stop('No depth_sample!!!')}
  
  date_type <- unique(used_age$date_type)
  #age_model <- unique(original_chrono$age_model_type)
  
  if(copra){
    COPRA_age <- original_chrono$interp_age
    COPRA_iqr <- original_chrono$interp_age_uncert_pos + original_chrono$interp_age_uncert_neg
    col_COPRA <- 'royalblue3'
  }
  
  #print('copra')
  
  if(stalage){ # in mm
    setwd(file.path(working_directory, file_name, '/StalAge'))
    StalAge <- read.csv('StalAge_chronology.csv', header = T)
    StalAge_age <- StalAge$StalAge_age
    StalAge_iqr <- StalAge$StalAge_age_uncert_pos + StalAge$StalAge_age_uncert_neg
    col_StalAge <- 'forestgreen'
  }
  #print('stalage')
  
  if(b){ # in cm
    setwd(file.path(working_directory, file_name, '/Bacon_runs'))
    bacon <- read.csv('bacon_chronology.csv', header = T)
    #bacon <- bacon[-500,]
    bacon_age <- bacon$bacon_age
    bacon_iqr <- bacon$bacon_age_uncert_pos + bacon$bacon_age_uncert_neg
    col_bacon <- 'palevioletred4'
    #col_bacon <- 'yellow2'
  }
  #print('bacon')
  
  if(bc){ # in cm 
    setwd(file.path(working_directory, file_name,'/Bchron'))
    bchron <- read.csv('bchron_chronology.csv', header = T)
    bchron_age <- bchron$bchron_age
    bchron_iqr <- bchron$bchron_age_uncert_pos + bchron$bchron_age_uncert_neg
    col_bchron <- 'sienna1'
  }
  #print('bchron')
  
  if(linreg){ # in mm
    setwd(file.path(working_directory, file_name, '/linReg'))
    linReg <- read.csv('linReg_chronology.csv', header = T)
    lin_reg_age <- linReg$lin_reg_age
    lin_reg_iqr <- linReg$lin_reg_age_uncert_pos + linReg$lin_reg_age_uncert_neg
    
    col_linReg <- 'springgreen2'
    
  }
  
  if(lininterp){
    setwd(file.path(working_directory, file_name, '/linInterp'))
    linInt <- read.csv('linInt_chronology.csv', header = T)
    lin_interp_age <- linInt$lin_interp_age
    lin_interp_iqr <- linInt$lin_interp_age_uncert_pos + linInt$lin_interp_age_uncert_neg
    
    #depth_sample_new <- sort(c(depth_sample, hiatus))
    
    #col_linReg <- 'springgreen2'
    col_linInt <- 'skyblue'
  }
  
  if(copRa){
    setwd(file.path(working_directory, file_name, '/copRa'))
    cR <- read.csv('copRa_chronology.csv', header = T)
    cR_age <- cR$copRa_age
    cR_iqr <- cR$copRa_uncert_pos + cR$copRa_uncert_neg
    
    #depth_sample_new <- sort(c(depth_sample, hiatus))
    
    #col_linReg <- 'springgreen2'
    col_cR <- 'royalblue3'
  }
  #print('mc')
  
  age <- cbind(depth_sample)
  iqr <- cbind(depth_sample)
  
  color <- NULL
  color_iqr <- NULL
  
  if(linreg){
    age <- cbind(age, lin_reg_age)
    iqr <- cbind(iqr, lin_reg_iqr)
    
    color <- c(color, col_linReg)
    color_iqr <- c(color_iqr, col_linReg)
  }
  
  if(lininterp){
    age <- cbind(age, lin_interp_age)
    iqr <- cbind(iqr, lin_interp_iqr)
    
    color <- c(color, col_linInt)
    color_iqr <- c(color_iqr, col_linInt)
  }
  
  if(stalage){
    
    age <- cbind(age, StalAge_age)
    iqr <- cbind(iqr, StalAge_iqr)
    
    color <- c(color, col_StalAge)
    color_iqr <- c(color_iqr, col_StalAge)
    
  }
  
  if(b){
    age <- cbind(age, bacon_age)
    iqr <- cbind(iqr, bacon_iqr)
    
    color <- c(color, col_bacon)
    color_iqr <- c(color_iqr, col_bacon)
  }
  
  if(bc){
    age <- cbind(age, bchron_age)
    iqr <- cbind(iqr, bchron_iqr)
    
    color <- c(color, col_bchron)
    color_iqr <- c(color_iqr, col_bchron)
  }
  
  if(copRa){
    age <- cbind(age, cR_age)
    iqr <- cbind(iqr, cR_iqr)
    
    color <- c(color, col_cR)
    color_iqr <- c(color_iqr, col_cR)
  }
  
  if(copra){
    age <- cbind(age, COPRA_age)
    iqr <- cbind(iqr, COPRA_iqr)
    
    color <- c(color, col_COPRA)
    color_iqr <- c(color_iqr, col_COPRA)
  }
  
  leg <- colnames(age)[2:length(colnames(age))]
  leg_iqr <- colnames(iqr)[2:length(colnames(iqr))]
  
  # pdf('am.pdf',6,4)
  # plot_AM_right(age,cbind(unknown_age$corr_age, unknown_age$depth_dating_new),
  #               cbind(used_age$corr_age, used_age$depth_dating_new), legend = leg,
  #               hiatus, date_type, x_lim, y_lim, steps_x, 
  #               steps_y, cex_legend = 0.8, color)
  # dev.off()
  # 
  # pdf('iqr.pdf',6,4)
  # plot_IQR(iqr, age, leg_iqr, x_lim, steps_x, 
  #          steps_y = 2 , cex_legend = 0.4, color_iqr)
  # dev.off()
  #print(leg)
  
  x_lim <- round(max(age[,2:ncol(age)], na.rm = T)/1000,digits = 0)
  if(!empty(data.frame(unknown_age)) && min(unknown_age$corr_age)< min(age[,2:ncol(age)], na.rm = T)){x_lim_bottom <- round(min(unknown_age$corr_age)/1000, digits = 0)} else {x_lim_bottom <- round(min(age[,2:ncol(age)], na.rm = T)/1000, digits = 0)}
  y_lim <- round(max(age[,1])+3, digits = 0)
  if(x_lim/10 > 1) {
    steps_x <- round(x_lim/10)
  } else {steps_x <- round(x_lim/10, digits = 0)}
  #steps_x <- round(x_lim/10)
  steps_y <- round(y_lim/5, digits =0)
  
  setwd(file.path(working_directory, file_name))
  pdf(paste('~/Documents/CI2019/amIQR',file_name,'.pdf', sep = ''),22.3, 13)
  #pdf(paste('iqr_',fname,'.pdf', sep = ''),9.5, 11)
  layout(matrix(c(1,1,2,3), nrow=1, ncol=4, byrow = T))#, widths = c(rep(lcm(9),3)), heights = c(rep(lcm(6.5),3)))
  #layout(matrix(c(1,2), nrow=2, ncol=1, byrow = TRUE), widths = lcm(22), heights = c(rep(lcm(16),1),5))
  par(mar = c(6,4,3,0), oma = c(2,6,0,1))
  plot_AM_right_depth(age,unknown_age, used_age, legend = leg,
                      hiatus, date_type, x_lim_bottom = 0, cex_legend = 0.2, color, file_name)
  
  plot_IQR_depth(iqr, age, used_age, leg_iqr, cex_legend = 0.4, color_iqr, file_name, hiatus)
  
  plot.new()
  #legend("top",legend = c("MC-ICP-MS U/Th age - used", 'Age - not used', 'Hiatus age'), pch = c(4,25,4),col = c("black", "red", 'orange') ,bty='n', ncol = 1, title = "Age-depth model information", cex = 1.3, pt.cex = 1.3, title.adj = 0.5)#, horiz = T)
  #legend('center', legend= c('Original AM: unknown', 'Lin. Reg.', 'Lin. Interp.','StalAge','Bacon','Bchron','COPRA', 'Hiatus depths', 'Dating depths'), 
  #       lwd = 3,col = c(color, 'orange', 'black'), lty = c(rep(1, length(color)),3,3),bty ='n', cex = 1.3, ncol=1) #title = "Additional AM information"
  #legend("top",legend = c("Event: actively dripping (mask)","MC-ICP-MS U/Th age - used", 'Age - not used', 'Hiatus age','' ,'Original AM: unknown', 'Lin. Reg.', 'Lin. Interp.','StalAge','Bacon','Bchron','COPRA', 'Hiatus depth', 'Dating depth'), 
  #       pch = c(4,4,25,4,NA,NA,NA,NA, NA, NA, NA,NA,NA,NA),col = c("blue","black", "red", 'orange',NA,color, 'orange', 'black') ,bty='n', ncol = 1, title = "Age-depth model information", cex = 1.3, pt.cex = 1.3, title.adj = 0.5,
  #       lwd = 2, lty = c(0,0,0,0,0,rep(1, length(color)),2,3))
  legend("top",legend = c("MC-ICP-MS U/Th age - used", 'MC-ICP-MS U/Th age - not used', 'Synthetic hiatus age','', 'Lin.Reg.' ,'Lin. Interp.','StalAge','Bacon','Bchron','COPRA (original chronology)', 'Hiatus depth', 'Dating depth'), 
         pch = c(4,25,4,NA,NA, NA, NA, NA,NA,NA,NA, NA),col = c("black", "red", 'orange',NA,color, 'orange', 'black') ,bty='n', ncol = 1, title = "Age-depth model information", cex = 2.4, pt.cex = 5, title.adj = 0.5,
         lwd = 2, lty = c(0,0,0,0,1,1,1,1,1,1,2,3))
  dev.off()
  
}

plot_AM_right_depth <- function(age, unknown, used, legend, hiatus, date_type, x_lim_bottom, cex_legend, color, file_name) {
  j <- ncol(age)
  
  #age[,1] <- age[,1]*10 # depth to mm
  #unknown[,2] <- unknown[,2]*10 # depth to mm
  #used[,2] <- used[,2]*10 # depth to mm
  
  age[,seq(2,j)] <- age[,seq(2,j)]/1000
  unknown$corr_age <- unknown$corr_age/1000 # to ka
  used[,c(2,3)] <- used[,c(2,3)]/1000 # to ka
  #used[,2] <- used[,2]/1000 # to ka
  
  y_lim <- 130
  #y_lim <- max(age[,1]) + 10
  x_lim <- max(apply(age[,seq(2,j)], 2, function(x){max(x, na.rm = T)})) +1
  #x_lim <- 14
  y_steps <- 25
  #y_steps <- round(y_lim/5)
  x_steps <- round(x_lim/10)
  
  #par(mar = c(4,4,1,0))
  
  matplot(age[, seq(2,j)], age[,1], type = 'l', lty = 1, col = color, xlab = '', ylab = '', lwd = 6, axes = FALSE, ylim = c(210,0), xlim = c(6,16))
  #points(used[,2], used[,1], pch = 4, col = "orange", cex = 2, lwd = 2)
  #points(used[,2], used[,4], pch = 4, col = "black", cex = 2, lwd = 2)
  #arrows(used[,2]-2*used[,3], used[,4], used[,2]+2*used[,3], used[,4], length=0.05, angle=90, code=3, col = 'black')
  #points(used[1,2], used[1,4], pch = 4, col = "blue", cex = 2, lwd = 2)
  #arrows(used[1,2]-used[1,3], used[1,4], used[1,2]+used[1,3], used[1,4], length=0.05, angle=90, code=3, col = 'blue')
  
  h <- used %>% filter(dating_id %in% hiatus$sample_id)
  u <- used %>% filter(!(dating_id %in% hiatus$sample_id))
  points(h[,2], h[,4], pch = 4, col = "orange", cex = 6, lwd = 4)
  arrows(h[,2]-h[,3], h[,4], h[,2]+h[,3], h[,4], length=0.05, angle=90, code=3, col = 'orange')
  
  points(u[,2], u[,4], pch = 4, col = "black", cex = 6, lwd = 4)
  arrows(u[,2]-2*u[,3], u[,4], u[,2]+2*u[,3], u[,4], length=0.05, angle=90, code=3, col = 'black')
  
  segments(x0 = u$corr_age, y0 = u$depth_dating, x1 = x_lim, col = 'black', lty = 3, lwd= 3, outer = T)
  if(length(hiatus$depth_sample) != 0) {
    segments(x0 = h$corr_age, y0 = h$depth_dating, x1 = x_lim, col = 'orange', lty = 2, lwd = 3)
    #abline(h = hiatus$depth_sample, col = 'orange', lty = 3, lwd = 1)
    
    #par(cex=cex_legend)
    #legend('bottomleft', legend= c(legend,'Hiatus'), lwd = 2,col = c(color,'grey'), lty = c(rep(1,j-1),2),bty ='n', title = "Additional AM information")
  } else {
    #par(cex = cex_legend)
    #legend('bottomleft', legend= legend, col = color, lwd = 2,lty = 1, bty ='n', title = "Additional AM information")
  }
  
  if (length(unknown) != 0) {
    points(unknown[,2],unknown[,1], pch = 25, col ="red", cex = 6, lwd = 4)
    #par(cex = cex_legend)
    #legend("topright",legend = c(paste(date_type,"- used"),"U/Th ages - not used"), pch = c(4,25),col = c("orange", "black") ,bty='n', title = "AM date types")
  } else {
    #par(cex = cex_legend)
    #legend("topright",legend = c(paste(date_type,"- used")), pch = c(4),col = c("orange") ,bty = 'n', title = "AM date types")
  }
  
  #axis(side = 1,  lwd = 2, cex.axis = 2.2, at = seq(6,18,by = 2), lwd.ticks = 2,labels = seq(6,18, by = 2), padj = 0.5)
  #axis(side = 2,  lwd = 2, lwd.ticks = 2, cex.axis = 2.2, labels = seq(210,0, by = -10), at= seq(210, 0, by = -10))#, padj = -0.5)
  
  #mtext(side = 3, text = 'Age-depth model', line = 1, cex = 2)
  #mtext(side = 3, text = 'A', line = -1, cex = 2, padj = 1, adj = 0.03)
  #mtext(side = 2, text = 'Depth from top [mm]', line = 3, cex  =1.3)
  #mtext(side = 1, text = 'Median age [kyrs BP]', line = 3, cex  =1.3)
  
  axis(side = 1,  lwd = 5, cex.axis = 5, at = seq(6,16,by = 2), lwd.ticks = 4,labels = seq(6,16, by = 2), padj = 0.5)
  #rug(x =  7:16, ticksize = -0.03, side = 1, lwd = 2)
  axis(side = 2,  lwd = 5, lwd.ticks = 4, cex.axis = 5, labels = seq(210,0, by = -30), at= seq(210, 0, by = -30))#, padj = -0.5)
  #rug(x = seq(210,0,-10), ticksize = -0.03, side = 2, lwd = 2)
  
  #mtext(side = 3, text = 'Age-depth model', line = 1, cex = 2)
  #mtext(side = 3, text = 'A', line = -1, cex = 2, padj = 1, adj = 0.03)
  mtext(side = 2, text = 'Depth from top [mm]', line = 2, cex  =3, outer = T)
  mtext(side = 1, text = 'Median age [kyrs BP]', line = 4, cex  =3, padj = 0.9)#, outer = T)
  
}

plot_IQR_depth <- function(iqr,age,used,legend, cex_legend, color, file_name, hiatus) {
  j <- ncol(iqr)
  
  iqr <- iqr/1000
  #age[,seq(3,j+1)] <- age[,seq(3,j+1)]/1000
  
  y_lim <- 130
  #y_lim <- max(age[,1]) + 10
  #x_lim <- max(apply(iqr[,seq(2,j)], 2, function(x){max(x, na.rm = T)})) +1
  x_lim <- 2.5
  y_steps <- 25
  #y_steps <- round(y_lim/5)
  x_steps <- round(x_lim/10)
  
  u <- used %>% filter(!(dating_id %in% hiatus$sample_id))
  
  matplot(iqr[,seq(2,j)], age[,1], type = 'l', lty = 1, col = color, xlab = '', ylab = '', lwd = 6, axes = FALSE,ylim = c(210,0), xlim = c(0.001, 2), log = 'x')
  #legend('topright', legend= legend, col = color, lwd = 2,lty = 1, bty ='n', title = "Information")
  abline(h = u$depth_dating, col = 'black', lty = 3, lwd = 3)
  abline(h = hiatus$depth_sample, col = 'orange', lty = 2, lwd = 3)
  #points(x = used/1000, y = rep(0, length(used)), pch = 4, col = "orange", cex = 2, lwd = 2)
  
  #if(length(hiatus_ages) != 0) {
  #  abline(v = hiatus_ages, col = 'grey', lty = 2, lwd = 1.5)
  #}
  
  #axis(side = 1,  at = seq(0, x_lim, by = 5), labels = seq(0,x_lim, by = 5), lwd = 2, lwd.ticks = 2, cex.axis = 2.2, padj = 0.5)
  #atz<-c(0.1,0.5,1,2,5,10,15)
  #atz<-axTicks(1)
  #axis(side = 1,  at = atz, labels = atz, lwd = 2, lwd.ticks = 2, cex.axis = 2.2, padj = 0.5)
  axis(side = 1,  at = c(0.001,0.1,0.5,1,2), labels = c(0.001,0.1,0.5,1,2), lwd = 5, lwd.ticks = 4, cex.axis = 5, padj = 0.5)
  
  #axis(side = 2,  lwd = 2, lwd.ticks = 2, cex.axis = 2.2, labels = seq(0, round(y_lim),by = y_steps), at= seq(0,round(y_lim), by = y_steps))#, padj = -0.5)
  #axis(side = 2,  lwd = 2, lwd.ticks = 2, cex.axis = 2, labels = seq(0,1.5, by = 0.25), at= seq(0,1.5, by = 0.25))#, padj = -0.5)
  
  #mtext(side = 3, text = 'Age uncertainties', line = 1, cex = 2)
  #mtext(side = 3, text = 'B', line = -1, cex = 2, padj = 1, adj = 0)
  #mtext(side = 2, text = 'Depth from top [mm]', line = 4, cex = 1.7)
  #mtext(side = 1, text = 'IQR [kyrs]', line = 3, cex  =1.3)
  mtext(side = 1, text = 'IQR [kyrs]', line = 4, cex  = 3, padj = 0.9)#), outer = T)
}
