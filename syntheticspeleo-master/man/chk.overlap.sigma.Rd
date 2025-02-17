% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{chk.overlap.sigma}
\alias{chk.overlap.sigma}
\title{Check the overlap of consecutive ages}
\usage{
chk.overlap.sigma(AgeEst, UncEst = rep(0, length(AgeEst)))
}
\arguments{
\item{AgeEst}{vector of ages}

\item{UncEst}{vector of age uncertainties (1 sd)}
}
\value{
logical
}
\description{
Check the overlap of consecutive ages within 2 standard deviations, i.a. whether there are reversals or not
}
\examples{
AgeEst<-seq(0,1000,by=100);
AgeEstUnc<-AgeEst*0.1
chk.overlap.sigma(AgeEst)
chk.overlap.sigma(AgeEst,AgeEstUnc)
AgeEstMod<-AgeEst[c(1,2,3,5,4,6,7,8,9,10)]
chk.overlap.sigma(AgeEstMod)
}
