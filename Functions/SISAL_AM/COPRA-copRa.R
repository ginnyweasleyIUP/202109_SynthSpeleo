setwd("~/Documents/copRa/Test Core/COPRA")

r <- read.table('d_TestCoreCarla_17-Jul-2019.dms', header = F, fill=T)


setwd("~/Documents/copRa/Test Core")
Bacon('test', cc = 0, depths.file = T, hiatus.depths = c(100,300), acc.mean = 10, ask = F)
setwd("~/Documents/copRa/Test Core/Bacon_runs/test")
output <- read.table(paste(info$prefix, ".out", sep = "")) # 6000
log_like <- output[,length(output)]
ac <- Acf(log_like, lag.max = 2000)
n <- which(ac$acf[,1,1]<0.3)[1] -1
j <- 2000*n
print(paste('j=',j))
setwd("~/Documents/copRa/Test Core")
Bacon('test', cc = 0, depths.file = T, hiatus.depths = c(100,300), acc.mean = 10, ask = F, ssize = j)
setwd("~/Documents/copRa/Test Core/Bacon_runs/test")
output <- read.table(paste(info$prefix, ".out", sep = "")) # 6000
thinner(1- 1/n)

j <- 2000
m <- 1
while(n > 1){
  setwd("~/Documents/copRa/Test Core")
  Bacon('test', cc = 0, depths.file = T, hiatus.depths = c(100,300), acc.mean = 10, ask = F, ssize = j)
  if(m>1){thinner(1-1/n)}
  setwd("~/Documents/copRa/Test Core/Bacon_runs/test")
  output <- read.table(paste(info$prefix, ".out", sep = "")) # 6000
  log_like <- output[,length(output)]
  ac <- Acf(log_like, lag.max = 2000)
  n <- which(ac$acf[,1,1]<0.3)[1] -1
  j <- 2000*n
  print(paste('j=',j))
  m <- m+1
}
