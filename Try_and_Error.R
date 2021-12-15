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

p.min = 0.1
p.max=1
p.n=5
t.min=0
t.max=10000
t.n=4
d.min=0
d.max=1000
h.length=c(2000,3000,500)
file_path1 ='~/synSpeleo/test/'
file_path2 = '~/synSpeleo/test/summary'
rate=3
t.res=2
nsamp=15
n=100

# Try out the generate speleothem function

speleo_test<-generate.speleo(t.top=t.min,t.bottom=t.max,z.top = d.min,z.bottom=d.max,t.res=t.res)
tlist <- c(t.min,c(cumsum(c(rep((t.max-t.min)/t.n,t.n))[sample(t.n)])))
rates <- c(1,3,1,3)

speleo_test<-growthfcn.changinggamma.new(speleo_test, tlist = tlist, rates = rates)

speleo_test$datingparams$tdofs<-10
speleo_test$datingparams$precision<-0.1
speleo_test$datingparams$topage<-10
speleo_test$datingparams$nsample<-10

speleo_test<-do.dating(speleo_test)

plotgrowthfcn(speleo_test)
plotdating(speleo_test$Dtable,add=TRUE)

