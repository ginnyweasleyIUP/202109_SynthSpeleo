#' Location-Scale formulation of the Student t distribution
#' @name TDist_ls
#' @aliases pt_ls
#' @aliases dt_ls
#' @aliases qt_ls
#' @aliases rt_ls
#' @description Location-Scale formulation of the student t-distribution
#' @param x quantiles
#' @param prob vector of quantiles
#' @param df degrees of freedom
#' @param mu mean (location)
#' @param a scale (standard deviation for df=infinity)
#' @import stats
#'
#' @return Analogous to the \code{\link{TDist}} implementation, \code{dt_ls} gives the density, \code{pt_ls} the distribution function, \code{qt_ls} the quantile function and \code{rt_ls} the random deviations
#'
#' @examples
#' curve(exp=dnorm(x,mean=1,sd=3),from=-15,to=15,col="red",xlab="x",ylab="density")
#' curve(expr =dt_ls(x,a=3,mu=1,df=3),from=-15,to=15,add=TRUE)
#' curve(expr =dt_ls(x,a=3,mu=1,df=1),from=-15,to=15,add=TRUE,lty=2)
#' curve(expr =dt_ls(x,a=5,mu=1,df=3),from=-15,to=15,add=TRUE,lty=3)
#' curve(expr =dt_ls(x,a=3,mu=2,df=3),from=-15,to=15,add=TRUE,lty=4)
#' curve(expr=dt_ls(x,a=3,mu=2,df=30),add=TRUE,col="orange",lty=1)
#' legend("topleft",lty=seq(1,4),c("dt_ls:a=3,mu=1,df=3","dt_ls:a=3,mu=1,df=1","dt_ls:a=5,mu=1,df=3","dt_ls:a=3,mu=2,df=3","dnorm:mu=1,sd=3","dt_ls:a=3,mu=2,df=30"),col=c(1,1,1,1,"red","orange"))
#' @rdname TDist_ls
#' @references \url{https://en.wikipedia.org/wiki/Location-scale_family}
#' @export
dt_ls <- function(x, df, mu, a) {
  1/a * dt((x - mu)/a, df)
}
##' @rdname TDist_ls
##' @return \code{pt_ls()} returns probabilities
##' @export
pt_ls <- function(x, df, mu, a) {pt((x - mu)/a, df)}
##' @rdname TDist_ls
##' @return \code{qt_ls()} returns quantiles
##' @export
qt_ls <- function(prob, df, mu, a) {qt(prob, df)*a + mu}
##' @rdname TDist_ls
##' @param n integer
##' @return \code{rt_ls()} returns observations
##' @export
rt_ls <- function(n, df, mu, a) {rt(n,df)*a + mu}



#' Generate synthetic speleothem record
#' @name generate.speleo
#' @param z.top top depth (assumed: in millimeters)
#' @param z.bottom bottom (assumed: in millimeters)
#' @param t.top top age (assumed: in years before present)
#' @param t.bottom bottom age (assumed: in years before present)
#' @param growthfcn growth object (dummy parameter)
#' @param proxy proxy object (dummy parameter)
#' @param t.res temporal resolution of the initial simulation
#' @param Dtable dating table
#' @param name Name of the speleothem
#' @description Generate a synthetic speleothem record
#' @return speleo-object (list)
#' @seealso \code{\link{generate.speleo}} \code{\link{plotgrowthfcn}} \code{\link{do.dating}}
#' @export
#'
#' @examples
#' speleo<-generate.speleo()
#' speleo<-growthfcn.constantgamma(speleo)
#  growthfcn<-growthfcn.constantgamma()
#' plotgrowthfcn(speleo)
#' speleo<-do.dating(speleo)
#' plotdating(speleo$Dtable,add=TRUE)
#' export.files(speleo)
generate.speleo<-function(t.res=1,z.top=0,z.bottom=1000,t.top=0,t.bottom=800,growthfcn=NA,Dtable=NA,proxy=NA,name="Synthetic"){
  speleo<-list()
  speleo$name=name
  speleo$z.top=z.top
  speleo$t.res=t.res
  speleo$z.bottom=z.bottom
  speleo$t.top=t.top
  speleo$t.bottom=t.bottom
  speleo$growthfcn=NULL
  speleo$true.age=seq(t.top,t.bottom,by=t.res)
  speleo$growthparams<-NULL
  speleo$datingparams<-list(nsample=10,tdofs=30,precision=50,nfac=100,use_updated_Dtable_for_agemodel=FALSE,topage=t.top)
  speleo$Dtable=Dtable
  speleo$proxy=proxy
  return(speleo)
}



#' @name growthfcn.constantgamma
#' @title Gamma-distributed growth function with constant mean growth
#' @description Simulate the growth of a constantly accumulating (speleothem) palaeoclimate record.
#' @param skew skewness of the growth distribution
#' @param speleo \code{speleo}-object
#' @return speleo-object \code{speleo}-object with updated \code{growthfcn}
#' @export
#' @seealso \code{\link{generate.speleo}} \code{\link{plotgrowthfcn}} \code{\link{do.dating}}
#' @examples
#' speleo<-generate.speleo()
#' growthfcn<-growthfcn.constantgamma(speleo)
#' plotgrowthfcn(growthfcn)
#'
growthfcn.constantgamma<-function(speleo,skew=1){
  speleo$growthparams<-list()
  speleo$growthparams$growthrate=(speleo$z.bottom-speleo$z.top)/(speleo$t.bottom-speleo$t.top)
  # give the expected mean growth per year, for constant accumulation/growth this is always the same
  speleo$growthrate.exp<-rep(speleo$growthparams$growthrate*speleo$t.res,length(speleo$true.age))
  speleo$growthparams$skew=skew
  speleo$growthparams$shape=4/(speleo$growthparams$skew)^2;
  mn<-speleo$growthparams$growthrate*speleo$t.res
  speleo$growthparams$scale=mn/(speleo$growthparams$shape);
  L=length(speleo$true.age)
  # generate more random variables than necessary
  print(speleo$growthparams$growthrate)
  R=rgamma(L,shape=speleo$growthparams$shape,scale=speleo$growthparams$scale) # vector with gamma-distributed waiting times, with expected mean of the growthrate

  #growth.real<-cumsum(R*speleo$growthrate.exp)*speleo$growthrate
  growth.real<-cumsum(R)
  expected.z<-cumsum(speleo$growthrate.exp)

  speleo$growthfcn<-data.frame(cbind(speleo$true.age,expected.z,growth.real))
  colnames(speleo$growthfcn)[1]<-"true.age"
  return(speleo)
}
# The next step is to sample the age radiometrically at certain depths, and with a certain analytical precision
# To set up "breakpoint" and "hiatus" models I need to figure out how to specify the parameters for those sensibly (maybe either specify time coverage OR growth rates fro segments)



#' @name growthfcn.changinggamma
#' @rdname growthfcn.constantgamma
#' @inheritParams growthfcn.constantgamma
#' @param tlist vector of start and end times for growth rate changes
#' @param rates slopes (must be of length \code{L=length(tlist)-1})
#' @examples
#' speleo<-generate.speleo()
#' growthfcn<-growthfcn.changinggamma(speleo)
#' plotgrowthfcn(growthfcn)
#' @export
growthfcn.changinggamma<-function(speleo,tlist=c(speleo$t.top,(speleo$t.bottom-speleo$t.top)/4,(speleo$t.bottom-speleo$t.top)/2,3*(speleo$t.bottom-speleo$t.top)/4,speleo$t.bottom),rates=c(1,0,2,4),skew=1){
  if (!(length(rates)==(length(tlist)-1))) stop("check defined rats & change times for growth rate")

  speleo$growthparams<-list()

  #hialength<-diff(t.hiat)[seq(1,length(t.hiat)-1,by=2)]

  speleo$growthparams$growthrate=(speleo$z.bottom-speleo$z.top)/(speleo$t.bottom-speleo$t.top)
  # give the expected mean growth per year, for constant accumulation/growth this is always the same

  #timmod<-c(speleo$t.top,speleo$t.bottom-sum(hialength),by=speleo$t.res)

  tlist.ageunits<-speleo$true.age[sapply(tlist,function(x,tpts){which.min(abs(tpts-x))},tpts=speleo$true.age)]


  a.units<-findInterval(speleo$true.age,tlist,all.inside=TRUE) # to which growth rate section does the age belong
  grrts<-rep(NA,length(a.units))
  for (i in 1:length(rates)){
    grrts[a.units==i]<-rates[i]
  }

  fac<-length(grrts)/(length(grrts)-sum(grrts==0))
  grrts<-(grrts/mean(grrts))*speleo$growthparams$growthrate*speleo$t.res*fac #scale to mean growth rate

  speleo$growthrate.exp<-grrts



   speleo$growthparams$skew=skew
   speleo$growthparams$shape=4/(speleo$growthparams$skew)^2;
   speleo$growthparams$scale<-rep(NA,length(grrts))
   R<-rep(NA,length(grrts))
   for (i in 1: length(grrts)){
     mn<-grrts[i]
     speleo$growthparams$scale[i]=mn/(speleo$growthparams$shape);
     L=1
     # # generate more random variables than necessary
     # print(speleo$growthparams$growthrate)
     R[i]<-rgamma(1,shape=speleo$growthparams$shape,scale=speleo$growthparams$scale[i])
     # R=rgamma(L,shape=speleo$growthparams$shape,scale=speleo$growthparams$scale) # vector with gamma-distributed waiting times, with expected mean of the growthrate
   }
  # #growth.real<-cumsum(R*speleo$growthrate.exp)*speleo$growthrate
   growth.real<-cumsum(R)
   expected.z<-cumsum(speleo$growthrate.exp)
  #
   speleo$growthfcn<-data.frame(cbind(speleo$true.age,expected.z,growth.real))
   colnames(speleo$growthfcn)[1]<-"true.age"
  return(speleo)
}

#' @name growthfcn.changinggamma.new
#' @rdname growthfcn.constantgamma.new
#' @inheritParams growthfcn.constantgamma
#' @param tlist vector of start and end times for growth rate changes
#' @param rates slopes relative to each other (must be of length \code{L=length(tlist)-1})
#' @examples
#' speleo<-generate.speleo()
#' growthfcn<-growthfcn.changinggamma(speleo)
#' plotgrowthfcn(growthfcn)
#' @export
growthfcn.changinggamma.new<-function(speleo,tlist=c(speleo$t.top,(speleo$t.bottom-speleo$t.top)/4,(speleo$t.bottom-speleo$t.top)/2,3*(speleo$t.bottom-speleo$t.top)/4,speleo$t.bottom),rates=c(1,0,2,4),skew=1){
  if (!(length(rates)==(length(tlist)-1))) stop("check defined rats & change times for growth rate")

  speleo$growthparams<-list()


  speleo$growthparams$growthrate=(speleo$z.bottom-speleo$z.top)/(speleo$t.bottom-speleo$t.top) # give the expected mean growth per year, for constant accumulation/growth this is always the same


  tlist.ageunits<-speleo$true.age[sapply(tlist,function(x,tpts){which.min(abs(tpts-x))},tpts=speleo$true.age)]


  a.units<-findInterval(speleo$true.age,tlist,all.inside=TRUE) # to which growth rate section does the age belong

  n.rates <- length(rates)
  n.a.units <- length(a.units)

  gr <- (speleo$z.bottom - speleo$z.top)/sum((rates*diff(tlist)))

  grrts<-rep(NA,length(a.units))
  for (i in 1:length(rates)){
    grrts[a.units==i]<-rates[i] * gr
  }

  #fac<-length(grrts)/(length(grrts)-sum(grrts==0))
  grrts<-grrts*speleo$t.res #scale to exp. growth rate

  speleo$growthrate.exp<-grrts



  speleo$growthparams$skew=skew
  speleo$growthparams$shape=4/(speleo$growthparams$skew)^2;
  speleo$growthparams$scale<-rep(NA,length(grrts))
  R<-rep(NA,length(grrts))
  for (i in 1: length(grrts)){
    mn<-grrts[i]
    speleo$growthparams$scale[i]=mn/(speleo$growthparams$shape);
    L=1
    # # generate more random variables than necessary
    # print(speleo$growthparams$growthrate)
    R[i]<-rgamma(1,shape=speleo$growthparams$shape,scale=speleo$growthparams$scale[i])
    # R=rgamma(L,shape=speleo$growthparams$shape,scale=speleo$growthparams$scale) # vector with gamma-distributed waiting times, with expected mean of the growthrate
  }
  # #growth.real<-cumsum(R*speleo$growthrate.exp)*speleo$growthrate
  growth.real<-cumsum(R)
  expected.z<-cumsum(speleo$growthrate.exp)
  #
  speleo$growthfcn<-data.frame(cbind(speleo$true.age,expected.z,growth.real))
  colnames(speleo$growthfcn)[1]<-"true.age"
  return(speleo)
}

#' Plot the growth function of the synthetic speleothem
#'
#' @name plotgrowthfcn
#' @title plot growth function
#' @description Plot the growth function of the synthetic speleothem record
#' @param speleo \code{speleo}-object
#' @return \code{speleo}-object
#' @import zoo
#' @export
plotgrowthfcn<-function(speleo,xlim=F,ylim=F){
  orig<-zoo::zoo(speleo$growthfcn$expected.z,order.by=speleo$growthfcn$true.age)
  sim<-zoo::zoo(speleo$growthfcn$growth.real,order.by=speleo$growthfcn$true.age)
  orig[diff(speleo$growthfcn$growth.real)==0]<-NA
  sim[diff(speleo$growthfcn$growth.real)==0]<-NA
  if(xlim & ylim){
    plot(orig,col="indianred",type="l",xlab="",ylab="", ylim = ylim, xlim = xlim, axes=F)
    lines(sim)
  } else{
    plot(orig,col="indianred",type="l",xlab="years [a.u.]",ylab="depth [a.u.]", ylim = c(round(max(orig[,2], na.rm = T), digits = 0),0))
    lines(sim)
    #legend("topright",col=c("red","black"),c("expected","real"),bty="n",lty=c(1,1))
  }


}




#' Simulate a dating sample on a synthetic speleothem
#'
#' @name do.dating
#' @title do the dating of the synthetic speleothem
#' @param speleo speleo-object, containing \code{datingparams} list with \code{nsample}, the number of point-ages sampled, \code{precision}, the sample age uncertainty for t-distribution (width equivalent to standard deviation if \code{tdofs}=\code{infinity}), and \code{tdofs}, the degrees of freedom for the t-distribution
#'
#' @return speleo speleo-object with updated \code{Dtable} and \code{datingparams}
#'
#' @seealso \code{\link{generate.speleo}} \code{\link{plotgrowthfcn}} \code{\link{do.dating}}
#' @export
#' @import stats
do.dating<-function(speleo){
  #    Dtable<-data.frame(matrix(NA,ncol=3,nrow=10));colnames(Dtable)=c("sampling.depths","sampling.age.est","sampling.age.est.unc")
  #   # To draw from a normal, or long-tailed normal distribution, use the location/scale formulation of student's t-distribution
  #   #https://en.wikipedia.org/wiki/Location-scale_family#Converting_a_single_distribution_to_a_location%E2%80%93scale_family
  #   # location mu (analog to the mean), scale a (analog standard deviation), degrees of freedom df
  #   # for df=infinity the t-distribution is equal to the normal distribution
  #
  #   # dt_ls <- function(x, df, mu, a) 1/a * dt((x - mu)/a, df)
  #   # pt_ls <- function(x, df, mu, a) pt((x - mu)/a, df)
  #   # qt_ls <- function(prob, df, mu, a) qt(prob, df)*a + mu
  if (is.null(speleo$growthfcn)) stop("get growth function first")
  if (!is.list(speleo$datingparams)) stop("add dating parameters")

  sampling.depths<-c(0,seq(speleo$z.top,speleo$z.bottom,length.out = speleo$datingparams$nsample+2)[-c(1,speleo$datingparams$nsample+2)])
  #  which.min(abs(growthfcn$growth.real-sampling.depths))
  ind.samples<-sapply(sampling.depths,function(x,growthfcn){which.min(abs(growthfcn$growth.real-x))},growthfcn=speleo$growthfcn)

  sampling.age.true<-speleo$growthfcn$true.age[ind.samples]


  sampling.age.est<-sapply(sampling.age.true,function(x,tdofs,precision){rt_ls(1,df=tdofs,mu=x,a=precision)},speleo$datingparams$tdofs,speleo$datingparams$precision)
  sampling.age.est.unc<-c(1,rep(2*speleo$datingparams$precision,speleo$datingparams$nsample))
  ignore.dates<-increased.unc<-rep(NA,speleo$datingparams$nsample+1)
  date.tdofs<-c(100,rep(speleo$datingparams$tdofs,speleo$datingparams$nsample))
  dating_id <- seq(1000, 1000+speleo$datingparams$nsample,1)

  speleo$Dtable<-data.frame(cbind(dating_id,sampling.depths,sampling.age.est,sampling.age.est.unc,sampling.age.true,date.tdofs,ignore.dates,increased.unc))
  if (is.na(speleo$datingparams$topage)) speleo$Dtable<-speleo$Dtable[-1,]
  #speleo$datingparams<-list(nsample=speleo$datingparams$nsample,precision=precision,tdofs=tdofs)
  if (any(diff(speleo$Dtable$sampling.age.est)<0)) warning("Inversions!")
  return(speleo)
}

#' Plot the dating on top of the growth function of the synthetic speleothem
#' @name plotdating
#' @description Plot the dating table information on the growth function
#' @param Dtable from \code{speleo}-object and/or growth function
#' @param add logical, if TRUE plot.
#' @seealso \code{\link{generate.speleo}} \code{\link{plotgrowthfcn}} \code{\link{do.dating}}
#' @examples
#' speleo<-generate.speleo()
#' speleo<-growthfcn.constantgamma(speleo)
#  growthfcn<-growthfcn.constantgamma()
#' plotgrowthfcn(speleo)
#' speleo<-do.dating(speleo)
#' plotdating(speleo$Dtable,add=TRUE)
#' export.files(speleo)
#' @export
plotdating<-function(Dtable,add=FALSE){

  if (add==FALSE){
    plot(Dtable$sampling.age.true,Dtable$sampling.depths,bg=adjustcolor("grey",0.75),pch=22)
    points(Dtable$sampling.age.est,Dtable$sampling.depths,bg=adjustcolor("blue",0.75),pch=21)
    arrows(Dtable$sampling.age.true-Dtable$sampling.age.est.unc,Dtable$sampling.depths,Dtable$sampling.age.true+Dtable$sampling.age.est.unc,Dtable$sampling.depths,code=3,angle=90,length=0.05,col="grey")
    arrows(Dtable$sampling.age.est-Dtable$sampling.age.est.unc,Dtable$sampling.depths,Dtable$sampling.age.est+Dtable$sampling.age.est.unc,Dtable$sampling.depths,code=3,angle=90,length=0.05,col="blue")
  }
  if (add==TRUE){
    points(Dtable$sampling.age.true,Dtable$sampling.depths,bg=adjustcolor("grey",0.75),pch=22)
    points(Dtable$sampling.age.est,Dtable$sampling.depths,bg=adjustcolor("blue",0.75),pch=21)
    arrows(Dtable$sampling.age.true-Dtable$sampling.age.est.unc,Dtable$sampling.depths,Dtable$sampling.age.true+Dtable$sampling.age.est.unc,Dtable$sampling.depths,code=3,angle=90,length=0.05,col="grey")
    arrows(Dtable$sampling.age.est-Dtable$sampling.age.est.unc,Dtable$sampling.depths,Dtable$sampling.age.est+Dtable$sampling.age.est.unc,Dtable$sampling.depths,code=3,angle=90,length=0.05,col="blue")

  }
}

#' Write the files for the synthetic speleothem record to disk
#' @param speleo \code{speleo}-object
#' @param writepath valid path
#' @export
#' @importFrom utils write.csv
export.files<-function(speleo,writepath="./"){
  write.csv(speleo$Dtable,file=paste(writepath,"datingtable.csv",sep=""))
  write.csv(speleo$growthfcn,file=paste(writepath,"ageinformation.csv",sep=""))
  write("README",file=paste(writepath,"readme.txt",sep=""))
  textin<-paste("Synthetic speleothem record generated",date(), "\n with the following parameters:\n  Shape:",speleo$datingparams$shape,"\n Skewness",speleo$datingparams$skew,"\n tdofs (DOFs of the t-distribution):",speleo$datingparams$tdofs,"\n Precision",speleo$datingparams$precision)
  cat(textin, file=paste(writepath,"readme.txt",sep=""),append=TRUE)
  cat(paste("z.top:",speleo$z.top,"z.bottom:",speleo$z.bottom,"t.top:",speleo$t.top,"t.bottom:",speleo$t.bottom),file=paste(writepath,"readme.txt",sep=""),append=TRUE)
}



#' Kim & o'Neill temperature dependent calcite fractionation
#'
#' @param temperature temperature in kelvin
#' @param alpha fractionation factor R(A)/R(B)
#' @param delta difference sea water/calcite d18O in permil (Ganssen et al. 2011, referencing to Kim & o'Neil 1997)
#' @param tempdelta inverted delta
#' @references
#' Ganssen et al., (2011), CLim. Past http://dx.doi.org/10.5194/cp-7-1337-2011
#'
#' Kim, S.-T., & O’Neil, J. R. (1997). Equilibrium and nonequilibrium oxygen isotope effects in synthetic carbonates. Geochimica et Cosmochimica Acta, 61(16), 3461–3475. https://doi.org/10.1016/S0016-7037(97)00169-5
#' @return numeric value (if \code{temperature} is given, return \code{alpha}, if \code{alpha} is given, return \code{temperature}; if \code{delta} is given, return \code{temperature} )
#'
#' @examples
#' kim(alpha=kim(temperature=283.15))
#' curve(kim,from=273.15,to=300,xlab="temperature",ylab="alpha")
#' kim(delta=2)
#' # --> The Glacial-Interglacial change of approx. 2 permil would be equivalent to approximately 7K.
#' @export
kim<-function(temperature=NULL,alpha=NULL,delta=NULL, tempdelta=NULL){

  if (!is.null(temperature)){
    # in: temperature (in Kelvin)
    # out: fractionation factor (dimensionless)
    alpha<-exp((18.03*(1000/temperature)-32.42)/1000)
    return(alpha)
  }

  if (!is.null(alpha)){
    # in: fractionation factor [1,inf]
    # out: temperature it corresponds to
    temperature<-18.03*1000/(1000*log(alpha)+32.42)
    return(temperature)
  }
  if(!is.null(delta)){

    # in: difference in permil between host & calcite matrix
    # out: temperature (if this difference is due to fractionation)

    temperature<-16.1-4.64*(delta)+0.09*delta^2
    return(temperature)
  }
  if (!is.null(tempdelta)){
    # in: temperature
    # out: difference in permil between host & calcite matrix
    ret<-c(NA,NA)
    #ret[1]<-4.64/0.18+sqrt((4.64/0.09)^2/4+((tempdelta-16.1)/0.09))
    ret<-4.64/0.18-sqrt((4.64/0.09)^2/4+((tempdelta-16.1)/0.09))
    return(ret)
  }
}

#' Sample a proxy time series
#' @param nsamp \code{integer}
#' @param z.steps type of sampling, following \code{\link{BlockAverageMod}}. If \code{NULL}, consecutive sampling is expected. If equal to single number, this is the depth integrated in one sample. If a vector is given, that is expected to give the depths integrated in each sample individually.
#' @param stype \code{character}, currently only "gaussian" or "sinusoid" (10 cycles), if you supply a vector of length L=\code{length(speleo$true.age)} it is used as signal
#' @param prx.unc \code{numeric} proxy uncertainty
#' @param speleo \code{speleo}-object
#' @return \code{speleo}-object
#' @param ... optional
#' @export
sample.proxy<-function(speleo,nsamp=50,z.steps=NULL,stype="gaussian",prx.unc=1,...){
  proxy<-list(proxymat=NULL,misc=list(nsamp=nsamp,stype=stype,prx.unc=prx.unc,truevals=rep(NA,length(speleo$true.age))),ages=NULL)

  proxymat<-data.frame(matrix(NA,ncol=6,nrow=nsamp))
  colnames(proxymat)<-c('sample_id',"depth_sample","proxy","proxy.sd","agetrue","width")
  if (is.numeric(stype)&length(stype)==length(speleo$true.age)){
    # use supplied time series
    proxy$misc$truevals<-stype
  } else {
  if (stype=="gaussian"){
    # depth-independent non-correlated noise
    proxy$misc$truevals<-rnorm(n=length(speleo$true.age),mean=1,sd=prx.unc)
  }
  if (stype=="sinusoid"){
    per<-(speleo$t.bottom-speleo$t.top)/10
    proxy$misc$truevals<-sin(2*pi*speleo$true.age/per)
  }
  }

  # now sample this vector!
    proxymat$sample_id <-seq(1,nsamp,1)
    proxymat$depth_sample<-seq(min(speleo$growthfcn$growth.real),max(speleo$growthfcn$growth.real),length.out = nsamp+2)[-c(1,nsamp+2)]

    proxymat$proxy<-zoo::coredata(BlockAverageMod(dts = proxymat$depth,tsin = zoo(proxy$misc$truevals,order.by=speleo$growthfcn$growth.real),blockwidth=z.steps))
    #cat(length(proxymat$depth),"\n",length(proxy$misc$truevals),"\n",length(z.steps),"\n")
    proxymat$agetrue<-zoo::coredata(BlockAverageMod(dts = proxymat$depth,tsin = zoo(speleo$true.age,order.by=speleo$growthfcn$growth.real),blockwidth=z.steps))
    #cat(length(speleo$true.age),"\n",length(speleo$growthfcn$growth.real),"\n",length(z.steps))

    proxymat$prx.unc<-rep(NA,length(proxymat$proxy))
    if ((length(z.steps)==1)||(length(z.steps)==length(proxymat$depth))) proxymat$width<-z.steps

  proxy$proxymat<-proxymat
  speleo$proxy<-proxy
  return(speleo)

}


#' check a dating table and classify the difficulty and treatment level
#'
#' @param speleo \code{speleo}-object
#' @param add \code{logical} indicating plotting or not
#' @return \code{integer}, 0 for trivial (monotonous), 1 for slight reversals, 2 for required treatment
#'
#' @export
check.Dtable<-function(speleo,add=TRUE){
  ndat<-length(speleo$Dtable$sampling.age.est)
  ind.rev<-which(diff(speleo$Dtable$sampling.age.est)<0)
  if (chk.overlap.sigma(speleo$Dtable$sampling.age.est,speleo$Dtable$sampling.age.est.unc))
  {
    if (chk.overlap.sigma(speleo$Dtable$sampling.age.est)){
      cat("no reversal at all")
      return(0)
    }
    else {
      cat ("reversal within uncertainties, increase nfac")
      #    speleo$proxy$misc$nfac<-10
      #    speleo$Dtable$ignore.dates[ind.rev]<-"R"
      return(1)
    }
  } else {
    #Dtable can't be modeled like this
    warning("Suggest to do Age model modification!")
    # if (is.null(speleo$DtableMod)) {speleo$DtableMod<-speleo$Dtable
    return(2)
  }
}



#' Update a dating table
#'
#' @param speleo \code{speleo}-object
#' @param ignore.dates \code{vector} of dates. Defaults to \code{NULL}, or a vector of \code{NA}-values. Dates to be ignored are to be set to not-\code{NA}
#' @param increased.unc \code{vector} of multipliers for the dating uncertainty. Set to ones or \code{NULL} for no change.
#' @param date.tdofs \code{vector} of degrees of freedom for the Monte-Carlo point age sampling. Set to large values (e.g. 100) for very gaussian, small values (e.g. 1) for very skewed sampling following a t-Distribution
#' @param add \code{logical} describing whether to plot or not.
#'
#' @return speleo-object with updated DtableMod
#'
#' @export
#' @importFrom grDevices adjustcolor
#' @importFrom graphics arrows legend lines matlines matpoints plot points
update.Dtable<-function(speleo,ignore.dates=NULL,increased.unc=NULL,date.tdofs=NULL,add=FALSE){

  if (is.null(speleo$DtableMod)) speleo$DtableMod<-speleo$Dtable
  if (all(sapply(c(ignore.dates,increased.unc,date.tdofs),is.null))) {cat("nothing to update");return(speleo)}

  if (!is.null(ignore.dates)){
    speleo$DtableMod$ignore.dates<-ignore.dates
    if (add){
      ind.ign<-which(!is.na(ignore.dates))
      points(speleo$DtableMod$sampling.age.est[ind.ign],speleo$DtableMod$sampling.depths[ind.ign],col="red",pch=4,lwd=2)
      }
  }
  if (!is.null(increased.unc)){
    speleo$DtableMod$increased.unc<-increased.unc*speleo$Dtable$sampling.age.est.unc
    if(add){
      #ind.rev<-chk.which.rev(speleo$Dtable$sampling.age.est,speleo$Dtable$sampling.age.est.unc)
      ind.inc<-which(increased.unc>1)
      #points(speleo$DtableMod$sampling.age.est[ind.inc],speleo$DtableMod$sampling.depths[ind.inc],col="orange",pch=4,lwd=2)
      #arrows(speleo$DtableMod$sampling.age.est[(ind.inc+1):(ind.inc-1)]-2*speleo$DtableMod$increased.unc[(ind.rev+1):(ind.rev-1)],speleo$DtableMod$sampling.depths[(ind.rev+1):(ind.rev-1)],speleo$DtableMod$sampling.age.est[(ind.rev+1):(ind.rev-1)]+2*speleo$DtableMod$increased.unc[(ind.rev+1):(ind.rev-1)],speleo$DtableMod$sampling.depths[(ind.rev+1):(ind.rev-1)],code = 3,angle=90,col="orange")
      arrows(speleo$DtableMod$sampling.age.est[ind.inc]-2*speleo$DtableMod$increased.unc[ind.inc],speleo$DtableMod$sampling.depths[ind.inc],speleo$DtableMod$sampling.age.est[ind.inc]+2*speleo$DtableMod$increased.unc[ind.inc],speleo$DtableMod$sampling.depths[ind.inc],code = 3,angle=90,col="orange")
    }

  }
  if (!is.null(date.tdofs)){
    ind.mod.dof<-which(!(date.tdofs==speleo$Dtable$date.tdofs))
    speleo$DtableMod$date.tdofs[ind.mod.dof]<-date.tdofs[ind.mod.dof]
    points(speleo$DtableMod$sampling.age.est[ind.mod.dof],speleo$DtableMod$sampling.depths[ind.mod.dof],col="cyan",pch=5,lwd=2)
  }

  speleo$datingparams$use_updated_Dtable_for_agemodel<-1


  return(speleo)

}



check.for.hiatus<-function(speleo,add=TRUE){
  slps<-(diff(speleo$Dtable$sampling.age.est))/diff(speleo$Dtable$sampling.depths)
  # outlier detection or (alternatively) threshold criterion? If slope is smaller than... and at the same time diff(time) larger than threshold
h.i<-which(sapply(slps,function(x,n,vec){x>n*sd(vec)},4,slps))
if(add) abline(h=speleo$Dtable$sampling.depths[h.i])
return(h.i)
}
# to add the hiatuses in the age model: add hiatus depths; draw hiatus age from uniform between date[before] and date[after]; draw left and right edges from uniform between date[before],hiatus.age and hiatus.age,date[after]
# insert for interpolation

#' Get simple age model
#'
#' @param speleo \code{speleo}-object
#'
#' @param itype interpolation type "linint" or "spline"
#' @param nrealiz \code{integer}-value
#' @param add \code{logical}, default \code{FALSE}, if true plots results
#' @param dtype date type "tdist" or "norm" for "t-distributed" or "gaussian-distributed" date sampling
#' @param nfac \code{integer} value, should be larger than \code{nrealiz}, much larger if there are reversals
#' @name get.simpleagemodel
#' @importFrom signal pchip
#' @export
get.simpleagemodel<-function(speleo,itype="pchip",nrealiz=10,add=FALSE,dtype="tdist",nfac=nrealiz*10){
  if (is.null(speleo$proxy)) stop("add proxy information")
  if (is.null(speleo$Dtable)) stop("add dating information")
  if (speleo$datingparams$use_updated_Dtable_for_agemodel==1){
    # use DtableMod
    DTMod<-speleo$DtableMod
  } else {
    # use Dtable
    DTMod<-speleo$Dtable
  }


  AMrealiz<-NA

  # remove dates for simulation
  if(any(!is.na(DTMod$ignore.dates))) DTMod<-DTMod[-which(!is.na(DTMod$ignore.dates)),]
  DT.t.unc<-DTMod$sampling.age.est.unc
  if (any(!is.na(DTMod$increased.unc))) DT.t.unc[which(!is.na(DTMod$increased.unc))]<-DTMod$sampling.age.est.unc[which(!is.na(DTMod$increased.unc))]*DTMod$increased.unc[which(!is.na(DTMod$increased.unc))]

  if (!chk.overlap.sigma(DTMod$sampling.age.est,DT.t.unc)) stop("update dating table - no monotonous modeling possible")
  DT.z<-DTMod$sampling.depths
  prxd<-speleo$proxy$proxymat$depth
  DT.t<-DTMod$sampling.age.est

  attmpt<-0
  resid<- -1
  catch<-1

  while(resid<0){
    attmpt<-attmpt+1
    cat("Attempt",attmpt,"\n")
    speleo$datingparams$nfac<-round((catch*nfac+nrealiz))*attmpt

    if (dtype=="norm"){
      DTrealiz<-replicate(nfac,mapply(function(x,sdin){rnorm(mean=x,sd=sdin,n=1)},x=DT.t,sdin=DTMod$sampling.age.est.unc))
    }

    else{#t-dist
      DTrealiz<-replicate(nfac,mapply(function(x,sdin,dfsin){rt_ls(mu=x,a=sdin,n=1,df=dfsin)},x=DT.t,sdin=DTMod$sampling.age.est.unc,dfsin=DTMod$date.tdofs))
    }

    colnames(DTrealiz)<-paste("R",seq(1,nfac),sep="")
    rownames(DTrealiz)<-paste("depth",DTMod$sampling.depths)

    rm.ind<-unique(which(diff(DTrealiz)<0,arr.ind=TRUE)[,2])
    catch<-length(rm.ind)/nfac
    resid<-((nfac-round(catch*nfac))-nrealiz)
    if (resid<0)  cat(length(rm.ind)/nfac*100,"percent to be removed due to reversals, updating nfac to larger than",nrealiz+catch*nfac, "by default now", speleo$datingparams$nfac)
    if (attmpt>3){stop("check dating table, too many attempts")}
  }




  if (add){
    matpoints(DTrealiz,DTMod$sampling.depths,col=adjustcolor("grey",0.5),pch=19,cex=0.5)
    if (length(rm.ind)>0)  matpoints(DTrealiz[,rm.ind],DTMod$sampling.depths,col=adjustcolor("indianred",0.5),pch=20,cex=0.5)
  }

  if ((resid>0)&(length(rm.ind)>0)){
    DTrealiz<-DTrealiz[,-rm.ind]
  }

  # http://blog.revolutionanalytics.com/2015/09/interpolation-and-smoothing-functions-in-base-r.html

  #agei<-DTrealiz[,1]
  getlinens<-function(agei,DT.z,prxd){y<-approx(x = DT.z,y=agei,xout=prxd,rule=1);return(y$y)}#,DT.z=DTMod$sampling.depths,prxd=speleo$proxy$proxymat$depth))
  getsplineens<-function(agei,DT.z,prxd){
    # cubic spline interpolation, but do not allow extrapolation
    #y<-spline(x=DT.z,x,xout=prxd)$y +zoo::na.approx(x=DT.z,x,xout=prxd, na.rm = FALSE)
#    yval<-spline(x=DT.z,y=agei,xout=prxd,method="hyman")$y;
    yval<-spline(x=DT.z,y=agei,xout=prxd)$y;
    yval[prxd>1.01*max(DT.z)]<-NA
    yval[prxd<0.25*min(DT.z)]<-NA
    #y<-na.spline(y,x=DT.z,y=x,xout=prxd) + 0*na.approx(vector, na.rm = FALSE)
    return(yval)
  }
  getpchpens<-function(agei,DT.z,prxd){
    # cubic spline interpolation, but do not allow extrapolation
    #y<-spline(x=DT.z,x,xout=prxd)$y +zoo::na.approx(x=DT.z,x,xout=prxd, na.rm = FALSE)
    #    yval<-spline(x=DT.z,y=agei,xout=prxd,method="hyman")$y;
    print(agei)
    print(DT.z)
    yval<-pchip(xi=DT.z,yi=agei,x=prxd)
    yval[prxd>1.01*max(DT.z)]<-NA
    yval[prxd<0.2*min(DT.z)]<-NA
    #y<-na.spline(y,x=DT.z,y=x,xout=prxd) + 0*na.approx(vector, na.rm = FALSE)
    return(yval)
  }

  if (itype=="linint"){
    AMrealiz<-apply(DTrealiz,2,getlinens,DT.z=DT.z,prxd=prxd)
  }
  if (itype=="spline"){
    AMrealiz<-apply(DTrealiz,2,getsplineens,DT.z=DT.z,prxd=prxd)
    # for (i in 1:nrealiz){
    #   getsplineens(x=DTrealiz[,i],DT.z=DTMod$sampling.depths,prxd=speleo$proxy$proxymat$depth)
    # }
  }
  if (itype=="pchip"){
    n <- dim(DTrealiz)[2]
    AMrealiz<-apply(DTrealiz,2,getpchpens,DT.z=DT.z,prxd=prxd)
    # for (i in 1:nrealiz){
    #   getsplineens(x=DTrealiz[,i],DT.z=DTMod$sampling.depths,prxd=speleo$proxy$proxymat$depth)
    # }
  }

  if (add==TRUE){
    matlines(AMrealiz,prxd,col=adjustcolor("grey",0.25),lty=1)
  } else {
    plot(DTMod$sampling.age.est,DTMod$sampling.depths,col="blue",pch=2)
    matpoints(DTrealiz,DT.z,pch=22,cex=0.5,bg=adjustcolor("black",0.1),col=1)
    matlines(AMrealiz,prxd,col=adjustcolor("grey",0.25),lty=1)
  }
  speleo$proxy$ages<-AMrealiz
  return(speleo)
}


#' Check the overlap of consecutive ages
#'
#' @name chk.overlap.sigma
#'
#' @description Check the overlap of consecutive ages within 2 standard deviations, i.a. whether there are reversals or not
#'
#' @param AgeEst vector of ages
#'
#' @param UncEst vector of age uncertainties (1 sd)
#' @return logical
#' @examples
#' AgeEst<-seq(0,1000,by=100);
#' AgeEstUnc<-AgeEst*0.1
#' chk.overlap.sigma(AgeEst)
#' chk.overlap.sigma(AgeEst,AgeEstUnc)
#' AgeEstMod<-AgeEst[c(1,2,3,5,4,6,7,8,9,10)]
#' chk.overlap.sigma(AgeEstMod)
#' @export
chk.overlap.sigma<-function(AgeEst,UncEst=rep(0,length(AgeEst))){
  # check whether a dating table is "feasible" to run an age model for without manipulations
  #
  ndat<-length(AgeEst)
  feasib<-AgeEst[1:(ndat-1)]-2*UncEst[1:(ndat-1)]<(AgeEst[2:ndat]+2*UncEst[2:ndat])
  return(all(feasib))
}


#' Check which ages are reversed, and which one should be removed
#' @name chk.which.rev
#' @inheritParams chk.overlap.sigma
#' @param DepthS vector of depths
#'
#' @return index
#' @export
chk.which.rev<-function(AgeEst,DepthS,UncEst=rep(0,length(AgeEst))){
  # check whether a dating table is "feasible" to run an age model for without manipulations, and which points are reversed
  #print(AgeEst)
  #print(DepthS)
  #print(UncEst)

  ndat<-length(AgeEst)
  feasib<-AgeEst[1:(ndat-1)]-2*UncEst[1:(ndat-1)]<(AgeEst[2:ndat]+2*UncEst[2:ndat])
  ind.rev<-which(!feasib)
  res<-lm(DepthS~AgeEst)


  #abline(res,col="purple")
  cfs <- coef(res)

  for (j in 1:length(ind.rev)){
    if (abs(residuals(res))[ind.rev[j]]<abs(residuals(res))[ind.rev[j]+1]) ind.rev[j]<-ind.rev[j]+1
  }


  return(ind.rev)
}




#' Block average a time series with either consecutive or non-consecutive blocks
#'
#' @param dts Middle points of the time steps for the block averaged process
#' @param tsin Higher resolution process to be block averaged, should be a zoo object
#' @param blockwidth block width for the sampling points dts. Defaults to NULL (consecutive blocks). If a numeric value is given, this is taken as the width (in dts units) around the dts-midpoints, if a vector (of the same length as dts) is given, these are the widthes around dts for the sampling
#'
#' @return \code{zoo}-object with the block-averaged time series
#' @export
#' @import zoo
#' @author Kira Rehfeld (added nonconsecutive option) & Raphaël Hébert (consecutive original)
#' @examples
#' dts<-seq(5,1050,by=55)
#' tsin<-zoo(rnorm(1000),order.by=seq(1,1000,by=1))
#'
#' ts.noncons<-BlockAverageMod(dts=dts,tsin=tsin,blockwidth=20) # an error occurs if dtmin is smaller than blockwidth/2
#' ts.cons<-BlockAverageMod(dts=dts,tsin=tsin)
#' ts.noncons.changing<-BlockAverageMod(dts=dts,tsin=tsin,blockwidth=seq(5,50,length.out=length(dts)))
#' plot(tsin)
#' lines(ts.noncons,col="blue",lwd=2)
#' lines(ts.cons,col="cyan")
#' lines(ts.noncons.changing,col="green")
BlockAverageMod <- function(dts,tsin,blockwidth=NULL){
  if(!is.zoo(tsin)){warning("Input a zoo object, or things may not work correctly.")}
  # Extracts the indices and coredata from the zoo object which the input sim should be.
  indexsim<-index(tsin)
  sim<-coredata(tsin)

  # Creates the bounds of the blocks to be averaged based on the input vector of time steps dts
  # We are assuming that the mid-points between the time steps are the bounds
  # Note that for the first and last points, the we are making symmetrical based on the differences
  # with the second and second-to-last points, respectively.
  ddts<-c(dts[1]+(dts[1]-dts[2])/2,dts+c(diff(dts)/2,(rev(dts)[1]-rev(dts)[2])))

  if (is.null(blockwidth)){

  # This creates the list which assigns each point from the input sim to an interval in ddts
  intervals<-list(findInterval(indexsim,ddts,all.inside = TRUE))

  # We take the average of each block
  aggsim<-aggregate(sim,by=intervals,FUN=mean)$x


  }
  if (!is.null(blockwidth)&length(blockwidth)>=1){
    #blockwidth gives the width of the blocks and the time steps are the centers
    centrs<-dts
    leftedges<-dts-blockwidth/2
    rightedges<-dts+blockwidth/2
#    if (!all(leftedges>ddts[1:length(leftedges)])||(!all(rightedges<ddts[2:(length(ddts))]))) stop()
    if (!(all(leftedges[2:(length(leftedges))]>rightedges[1:(length(leftedges)-1)])&(all(rightedges[2:(length(leftedges))]>leftedges[1:(length(leftedges)-1)])))) stop("overlapping samples!?")
    intervls<-list(findInterval(indexsim,sort(c(leftedges,rightedges)),all.inside = FALSE))
    aggsim.withall<-aggregate(sim,by=intervls,FUN=mean)$x
    #the breaks between left and right edges are included so far --> take them out by taking only every other value
    aggsim<-(aggsim.withall[-seq(1,length(aggsim.withall),by=2)])
  }

  # If the input sim indices indexsim do not cover the blocks ddts, then the length of aggsim will not match
  # the length of dts and some elements will be duplicated in the output, the result would therefore be incorrect
  if(length(aggsim)!=length(dts)){warning("The number of averaged blocks does NOT correspond to the input dts time steps - check left & right edges.")
    cat(length(dts),"\n",length(aggsim),"\n")
  }
#  cat(length(dts),length(aggsim),length(aggsim.withall))
  tsout<-zoo(aggsim,order.by=dts)
  return(tsout)

}


