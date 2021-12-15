source("~/R/syntheticspeleo-master/R/functions.R")
library(zoo)
library(Hmisc)
library(pracma)

# WORKFLOW
set.seed(0)
speleo<-NA
#speleo<-generate.speleo(t.top=-50,t.bottom=30000,z.top = 0,z.bottom=4500,t.res=2.5)
speleo<-generate.speleo(t.top=0,t.bottom=10000,z.top = 0,z.bottom=100,t.res=2.5)
#speleo<-growthfcn.constantgamma(speleo)
#speleo<-growthfcn.changinggamma.new(speleo)
#plotgrowthfcn(speleo)

#speleo<-growthfcn.changinggamma.new(speleo, tlist = c(-50,2000,5000,6000,15000, 17000,21000,25000,25300,30000), rates = c(4,0,3,0,1,0.1,1.3,0,1.5))
speleo<-growthfcn.changinggamma.new(speleo, tlist = c(0,500,2000,5000,10000), rates = c(4,2,0,1))
#speleo <- growthfcn<-growthfcn.constantgamma(speleo)
plotgrowthfcn(speleo)
# there is an inversion if precision is set too wide --> catch inversions by elimination, or automatically increase age uncertainty there
speleo$datingparams$tdofs<-1000
speleo$datingparams$precision<-100
speleo$datingparams$topage<-10
speleo$datingparams$nsample<-30
speleo$datingparams$nfac <- 10000
speleo<-do.dating(speleo)


plotdating(speleo$Dtable,add=TRUE)
speleo<-sample.proxy(speleo = speleo,z.steps = 1,stype = "sinusoid")
dind<-check.Dtable(speleo,add=TRUE)
if (dind==0) {
  speleo$DtableMod

  speleo<-get.simpleagemodel(speleo,add=TRUE,itype = "pchip")
}
if (dind==1) {
  # increase number of realizations & decrease dofs for the t-distr
  ind.rev<-chk.which.rev(speleo$Dtable$sampling.age.est,speleo$Dtable$sampling.depths)

  ignore.dates<-rep(NA,length(speleo$Dtable$sampling.age.est))
  ignore.dates[ind.rev]<-1
  date.tdofs<-rep(10,length(speleo$Dtable$sampling.age.est))
  date.tdofs[(ind.rev-1):(ind.rev+1)]<-2
  speleo<-update.Dtable(speleo,ignore.dates=ignore.dates,date.tdofs = date.tdofs,add=TRUE)
  speleo<-get.simpleagemodel(speleo,add=TRUE,nfac = 200,itype = "spline")


}
if (dind==2) {
  # "radical
  ind.rev<-chk.which.rev(AgeEst=speleo$Dtable$sampling.age.est,DepthS=speleo$Dtable$sampling.depths,speleo$Dtable$sampling.age.est.unc)
  ignore.dates<-rep(NA,length(speleo$Dtable$sampling.age.est))
  ignore.dates[ind.rev]<-1

  speleo<-update.Dtable(speleo,ignore.dates=ignore.dates,add=TRUE)
  get.simpleagemodel(speleo,add=TRUE,nfac=200)
}
export.files(speleo)

# plot a couple of realizations

matplot(speleo$proxy$proxymat$agetrue,speleo$proxy$proxymat$proxy,col="grey",type="l")
lines(apply(array(speleo$proxy$proxymat$agetrue),1,median),speleo$proxy$proxymat$proxy)
tmed<-apply(array(speleo$proxy$proxymat$agetrue),1,median)
prxmed<-approx(speleo$proxy$proxymat$agetrue[1],speleo$proxy$proxymat$proxy,tmed)$y
   lines(tmed,prxmed,col="blue")
lines(speleo$true.age,speleo$proxy$misc$truevals)


matplot(speleo$proxy$ages,speleo$proxy$proxymat$proxy,col="grey",type="l")
lines(apply(array(speleo$proxy$ages),1,median),speleo$proxy$proxymat$proxy)
tmed<-apply(speleo$proxy$ages,1,median)
prxmed<-approx(speleo$proxy$ages[,1],speleo$proxy$proxymat$proxy,tmed)$y
lines(tmed,prxmed,col="blue")
lines(speleo$true.age,speleo$proxy$misc$truevals)

# now construct one proxy table with constant depth (--> sample-specific age uncertainty) and one with regularly sampled time (--> time-specific uncertainty)
# then add the option of having additional proxy uncertainty
