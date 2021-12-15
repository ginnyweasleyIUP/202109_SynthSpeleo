scan_lin_interp<-function(dating_tb) {
  
  library(Hmisc)
  #library(rpanel)
  library(plotrix)
  
  #x11()
  #errbar(depth, age, age+error, age-error, xlab="Distance from top [mm]", ylab="Age [a]")
  #title(main="Original age data")
  
  dating_id <- dating_tb$dating_id
  age <- dating_tb$corr_age
  depth <- dating_tb$depth_dating
  error <- dating_tb$corr_age_uncert
  date_type <- dating_tb$date_type
  
  hilf<-0
  
  test<-array(NA, c(length(depth), length(depth)))	#Test-Array f?r Screening
  
  
  age<-age[order(depth)]		#Sortiert die Vektoren aufsteigend nach der Tiefe
  error<-error[order(depth)]
  depth<-depth[order(depth)]
  
  
  for (i in 1:length(depth)) {	#Schleife zum Testen auf Inversionen. Wenn Inversion, schreibe 1 in Matrix.
    
    j<-i+1
    
    while (j<=length(depth)) {
      
      if (age[j]+error[j]<age[i]-error[i]) test[i, j]<-1 else test[i, j]<-0	#wenn 2-sigma Fehler nicht ?berlappen, schreibe 1 ins Feld
      
      j<-j+1
      
    }
    
  }
  
  #print(test)
  
  max<-1		#Variable f?r das Maximum der Inversionen
  max_pos<-0	#Variable f?r den Punkt des Maximums
  
  
  repeat {
    
    for (i in 1:length(depth)) {
      
      if (sum(test[i,], na.rm=TRUE)>max) {	#Gehe Zeilen durch ... ## NA z?hlen nicht da na.rm = T 
        
        max<-sum(test[i,], na.rm=TRUE)
        max_pos<-i
        
      }
      
      if (sum(test[,i], na.rm=TRUE)>max) {	#Gehe Spalten durch ...
        
        max<-sum(test[,i], na.rm=TRUE)
        max_pos<-i
        
      }
      
    }
    
    if (max>1) {				#Wenn ein Punkt mit mehr als einem anderen nicht passt, Abfrage, was getan werden soll
      
      #x11()
      #errbar(depth, age, age+error, age-error, xlab="Distance from top [mm]", ylab="Age [a]")
      #title(main="Screening age data for major outliers")
      
      #radius<-(depth[length(depth)]-depth[1])/(age[length(depth)]-age[1])*error[max_pos]
      #draw.circle(x=depth[max_pos], y=age[max_pos], radius=radius, border="red")
      
      #panel<-rp.control(size=c(500, 120), title=paste("Point", max_pos, "is an outlier! What do you want to do?"), max_pos=max_pos)
      
      #rp.button(panel, action=delete, title="Delete!", quitbutton=TRUE, pos=c(0,0,500, 40))
      #rp.button(panel, action=enlarge, title="Enlarge!", quitbutton=TRUE, pos=c(0,41,500, 40))
      #rp.button(panel, action=close, title="Close", quitbutton=TRUE, pos=c(0,81,500, 40))
      
      #rp.block(panel)
      
      #if (choice==1) {			#Auswahl: L?sche Punkt!
      
      #test[max_pos,]<-NA
      #test[,max_pos]<-NA
      #age[max_pos]<-NA
      
      #errbar(depth, age, age+error, age-error, xlab="Distance from top [mm]", ylab="Age [a]")
      #title(main="Screening age data for major outliers")
      
      #}
      
      #if (choice==2) {			#Auswahl: Vergr??ere Fehler!
      
      if (max(test[max_pos,], na.rm=TRUE)==1) { ## was soll es sonst sein????? es muss mindestens 2 einser haben, da max>1
        
        hilf<-age[max_pos]
        
        for (i in (max_pos+1):length(depth)) {
          
          if (is.na(test[max_pos, i])) next ### wann kann das auftreten ???? wir gehen hier reihen durch 
          
          if (test[max_pos, i]==1 && age[i]<hilf) hilf<-age[i]	#; print(hilf)
          
          error[max_pos]<-age[max_pos]-hilf
          
        }
        
      }
      
      if (max(test[, max_pos], na.rm=TRUE)==1) {
        
        hilf<-age[max_pos]
        
        for (i in 1:(max_pos-1)) {
          
          if (is.na(test[i, max_pos])) next
          
          if (test[i, max_pos]==1 && age[i]>hilf) hilf<-age[i]	#; print(hilf)
          
          error[max_pos]<-hilf-age[max_pos]
          
        }
        
      }
      
      test[max_pos,]<-NA 
      test[,max_pos]<-NA
      
      #errbar(depth, age, age+error, age-error, xlab="Distance from top [mm]", ylab="Age [a]")
      #title(main="Screening age data for major outliers")
      
      #}
      
      #if (choice==3) {			#Auswahl: Weder noch -> Hinweis, dass das nicht m?glich ist!
      
      #rp.messagebox("This is not possible. Select another option, please!", title = "STOP!")
      
      #}
      
    }
    
    if (max==1) break
    
    max<-1
    
  }
  
  
  test[is.na(test)]<-0			#Ersetze NA's in Matrix durch Nullen.
  
  
  #x11()
  #errbar(depth, age, age+error, age-error, xlab="Distance from top [mm]", ylab="Age [a]")
  #title(main="Age data screened for major outliers")
  
  
  depth<-depth[!is.na(age)]		#L?sche rausgeworfene Punkte
  error<-error[!is.na(age)]
  age<-age[!is.na(age)]
  
  
  Daten<-data.frame(dating_id = dating_id, corr_age = age, corr_age_uncert = error, depth_dating=depth, date_type = date_type)	#Speichere Punkte in Dataframe
  
  return(Daten)
  
}


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

slope<-function(depth, age, error) {
  
  iter<-200
  count<-0
  
  age_sim<-0
  
  for (i in 1:iter) {
    
    for (k in 1:3) age_sim[k]<-rnorm(1, age[k], error[k]/2)	#Simulation der Alter
    
    fit<-lm(age_sim[1:3] ~ depth [1:3])			#Fit ?ber simulierte Alter
    attr(fit$coefficients, "names")<-NULL
    
    if (fit$coefficients[2]>0) count<-count+1	#wenn Steigung positiv, z?hle eins hoch
    
  }
  
  #print(count/iter)
  
  return(count/iter)
  
}


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

scan_fine_lin_interp<-function(dating_tb){
  
  #attach(Daten)
  dating_id <- dating_tb$dating_id
  age <- dating_tb$corr_age
  depth <- dating_tb$depth_dating
  error <- dating_tb$corr_age_uncert
  date_type <- dating_tb$date_type
  
  if(any(error == 0)){
    knull<- which(error == 0)
    for(m in knull){error[m] <- 0.1}
  } else {knull <- -1}
  
  
  #print(knull)
  
  #print(error)
  #x11()
  #errbar(depth, age, age+error, age-error, xlab="Distance from top [mm]", ylab="Age [a]")			#Plotten der Daten
  #title(main="Screening for minor outliers")
  
  #print(error)
  
  Anteil<-0							#Vektoren f?r das Testen der Fits
  counter<-0
  part<-0
  
  part[1]<-1							#Definition des Beteiligungs-Vektors
  part[length(depth)]<-1
  part[2]<-2
  part[length(depth)-1]<-2
  
  for (i in 3:(length(depth)-2)) part[i]<-3
  
  #print(part)
  
  
  test<-0								#Pruefvariable f?r ?bergeordneten Test
  pruef<-0							#Pruefvariable f?r andere Tests
  
  
  while (test==0) {						#Schleife f?r Iteration bis Fehler alle passen
    
    test<-1								#Setzen der allgemeinen Pr?fvariable auf 1
    
    for (i in 1:length(depth)) Anteil[i]<-0				#Nullsetzen von Anteil
    for (i in 1:length(depth)) counter[i]<-0			#Nullsetzen von counter
    
    
    for (i in 1:(length(depth)-2)) {
      
      #print(i)
      
      
      fit<-lm(age[i:(i+2)] ~ depth [i:(i+2)], weights=1/error[i:(i+2)])	#Fit ?ber 3-Punkt-Intervall
      attr(fit$coefficients, "names")<-NULL
      
      #curve(fit$coefficients[2]*x+fit$coefficients[1], from=depth[i], to=depth[i+2], add=TRUE, col=i)	#Plotten des Fits
      
      age_fit<-0							#Bestimmung des Alters des Fits an den jeweiligen Punkten
      age_fit[i:(i+2)]<-fit$fitted.values
      
      
      pruef<-0							#Null-Setzen der Pr?fvariable
      
      for (k in i:(i+2)) {						#Test ob linearer Fit m?glich
        
        if (age_fit[k]>age[k]+error[k] || age_fit[k]<age[k]-error[k]) pruef<-pruef+1	#wenn Punkt au?erhalb der Fehlergrenzen, z?hle Pruefvariable um 1 hoch
        
      }
      
      #print (pruef)
      
      if (pruef>0) {
        
        counter[i:(i+2)]<-counter[i:(i+2)]+1				#wenn der Fit nicht geht, setze counter um 1 hoch
        
      } else {
        
        if (fit$coefficients[2]<0) {					#wenn die Steigung des Fits negativ ist
          
          if (slope(depth[i:(i+2)], age[i:(i+2)], error[i:(i+2)])<0.2) counter[i:(i+2)]<-counter[i:(i+2)]+1	#wenn weniger als 30% der Fits eine positive Steigung haben, setze counter um 1 hoch
          
        }
        
      }		
      
    }
    
    #print(counter)
    
    Anteil<-counter/part
    
    #print(Anteil)
    
    for (i in 1:length(depth)) {				#gehe Anteil nach 1ern durch
      
      if (Anteil[i]==1) {					#wenn Anteil=1, dann vergr??ere Fehler um 10% und setze allgemeine Pr?fvariable auf 1
        
        error[i]<-error[i]*1.1
        test<-0
        
      }
      
    }
    
    #errbar(depth, age, age+error, age-error, add=TRUE)		#Plotten der Daten
    
  }
  
  #x11()
  #errbar(depth, age, age+error, age-error)			#Plotten der Daten
  #title(main="Age data screened for minor outliers")
  
  if(knull != -1){
    for(m in knull){error[m] <- 0}
  }
  #print(knull)
  #print(error)
  
  Daten<-data.frame(dating_id = dating_id, corr_age = age, corr_age_uncert = error, depth_dating=depth, date_type = date_type)
  
  #detach(Daten)
  #print(Daten)
  
  return(Daten)
  
}

mc_ages_sep_new <- function(linReg = F, linInterp = F, age, age_error, N, c14 = F, working_directory, file_name) { # N number of MC simulations
  #set.seed(2019-02-04)
  library(clam)
  
  number <- 0
  d <- 0
  age_ensemble_final <- NA
  
  if(linReg) {
    if(c14){
      age_ensemble <- apply(cbind(age,age_error), 1, function(x)
        sample(calibrate(x[1], x[2], postbomb = 4, cc = 3, graph = FALSE, yrsteps = 0.05)$calib[,1],N, replace =F, 
               prob = calibrate(x[1], x[2], postbomb = 4, cc = 3, graph = FALSE, yrsteps = 0.05)$calib[,2])) 
      age_ensemble_final <- age_ensemble
      
    } else {
      
      age_ensemble <- apply(cbind(age,age_error), 1, function(x) rnorm(n = N, mean = x[1], sd = x[2])) # calculate N deviates for each age 
      age_ensemble_final <- age_ensemble
      
    }
    return(age_ensemble_final)
  }
  
  
  if(linInterp){
    while(d < N){
      #print('hello')
      k<- N-d
      #print(k)
      
      number <- number+1
      #print(n)
      
      if(c14){
        age_ensemble <- apply(cbind(age,age_error), 1, function(x)
          sample(calibrate(x[1], x[2], postbomb = 4, cc = 3, graph = FALSE, yrsteps = 0.05)$calib[,1],k, replace =F, 
                 prob = calibrate(x[1], x[2], postbomb = 4, cc = 3, graph = FALSE, yrsteps = 0.05)$calib[,2])) 
        # age_ensemble_final <- age_ensemble
        
      } else {
        
        age_ensemble <- apply(cbind(age,age_error), 1, function(x) rnorm(n = k, mean = x[1], sd = x[2])) # calculate N deviates for each age 
        #age_ensemble_linReg <- age_ensemble
        
      }
      
      if(k==1){
        age_ensemble_diff <- diff(age_ensemble)
        run <- T
        
        if (any(age_ensemble_diff < 0 & !is.na(age_ensemble_diff))){
          run <- F
          #print('k=1')
        }
        
      } else {
        age_ensemble_diff <- apply(age_ensemble, 1, diff) # each column contains the derivatives for one MC run -> N columns
        del <- NULL #age_ensemble_copy <- age_ensemble
        run <- T
        #if (k ==2) {print('k=2')}
          for (i in seq(1,k)) {
            #print(age_ensemble_diff[,i])
            if(any(age_ensemble_diff[,i] <0)){
              if (is.null(dim(age_ensemble)[1])) {
                #age_ensemble <- age_ensemble[-i,]
                run <- F
              } else {
                del <- c(del,i)
                #print(age_ensemble)
              }
              
            }
            
          }
        if(!is.null(del)){age_ensemble <- age_ensemble[-del,]}
        #if(k==1 && any(diff(age_ensemble)<0)){age_ensemble<-NULL}
        
      }
      
      if (run) {
        if (k == N){
          age_ensemble_final <- age_ensemble
        } else {
          age_ensemble_final <- rbind(age_ensemble_final, age_ensemble)
        }
      }
      
      run <- T
      d <- dim(age_ensemble_final)[1]
      #diff_final <- apply(age_ensemble_final, 1, diff)
      #print(any(diff_final<0))
      if(is.null(d)){d <- 1}
      if(number>2000 && d > 100){return(age_ensemble_final)
        break
      } else if(number > 2000 && d < 100) {
        setwd(file.path(working_directory, file_name, '/linInterp'))
        write.csv(c(number,d),'mc_fail.csv', row.names = F)
        stop('ERROR: too many iterations, check data!')}
      #print(c('Dim lI:', dim(age_ensemble_linInt_final)[1]))
      #print(c('Dim lR:', dim(age_ensemble_linReg_final)[1]))
    }
    return(age_ensemble_final) # ensemble
  }
  
}

