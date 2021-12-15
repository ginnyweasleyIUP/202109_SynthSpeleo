# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
# https://en.wikipedia.org/wiki/Location-scale_family#Converting_a_single_distribution_to_a_location%E2%80%93scale_family

# WORKFLOW
set.seed(0)
speleo<-NA
speleo<-generate.speleo(t.top=0,t.bottom=1000,z.top = 0,z.bottom=1000,t.res=0.33)
speleo<-growthfcn.constantgamma(speleo)

speleo<-growthfcn.changinggamma(speleo)
#growthfcn<-growthfcn.constantgamma()
plotgrowthfcn(speleo)
# there is an inversion if precision is set too wide --> catch inversions by elimination, or automatically increase age uncertainty there
speleo$datingparams$tdofs<-1
speleo$datingparams$precision<-30
speleo$datingparams$topage<-0
speleo$datingparams$nsample<-15
speleo<-do.dating(speleo)


plotdating(speleo$Dtable,add=TRUE)
speleo<-sample.proxy(speleo = speleo,z.steps = 5,stype = "sinusoid")
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

matplot(speleo$proxy$ages,speleo$proxy$proxymat$proxy,col="grey",type="l")
lines(apply(speleo$proxy$ages,1,median),speleo$proxy$proxymat$proxy)
tmed<-apply(speleo$proxy$ages,1,median)
prxmed<-approx(speleo$proxy$ages[,1],speleo$proxy$proxymat$proxy,tmed)$y
lines(tmed,prxmed,col="blue")
lines(speleo$true.age,speleo$proxy$misc$truevals)

# now construct one proxy table with constant depth (--> sample-specific age uncertainty) and one with regularly sampled time (--> time-specific uncertainty)
# then add the option of having additional proxy uncertainty
