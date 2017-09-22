\name{LM.logit-package}
\alias{LM.logit-package}
\alias{LM.logit}
\docType{package}
\title{
\packageTitle{LM.logit}
}
\description{
\packageDescription{LM.logit}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{LM.logit}
\packageIndices{LM.logit}
The likelihood-model of logistic regression is a parametric probability density function of the parameters in the logistic regression. The density parameters are estimated by fitting the density to MCMC draws from each subset data-model likelihood function, and then the fitted densities are recombined.
}
\author{
\packageAuthor{LM.logit}

Maintainer: \packageMaintainer{LM.logit}
}
\references{
\itemize{
 \item \url{http://deltarho.org}
 \item {Qi Liu, Anindya Bhadra,  Bowei Xi, and William S. Cleveland, Likelihood modeling for big data analysis using divide and recombine methods}
}
}

\keyword{likelihood modeling, divide and recombine }
\seealso{
\code{\link{datadr}}
}
\examples{
\dontrun{
set.seed(100)
library(datadr)
library(sn)
library(BayesLogit)
library(MASS)
library(mvtnorm)
library(moments)
x <- matrix(rnorm(1000*5), ncol=5)
ttheta <- rep(1,5)
y <- rbinom(1000, 1, 1 / (1 + exp(- x \%*\% ttheta)))
df <- cbind(y,x)
df <- as.data.frame(df)
names(df) <- c("y", "x1","x2","x3","x4","x5")
df_ddf <- ddf(df)

# in memory backend
 df_div <- divide(df_ddf, by =rrDiv(500))
 rst<- drml(df_div, y~x1+x2, size =1000, burnin =50, approx_method = "SN") 
 pred <- predNew.local(rst, y~x1+x2, df, 1000)
 
 # disc backend
 tmpdir <- "./tmp"
 DiskConn <- localDiskConn(file.path(tmpdir, "KV"), autoYes = TRUE)
 addData(DiskConn, df_div)
 DiskConn <- ddf(DiskConn)
 DiskConn <- updateAttributes(DiskConn)
 rst<- drml(DiskConn, y~x1+x2, size =1000, burnin =50, approx_method = "SN") 
 DiskConn_output <- localDiskConn(file.path(tmpdir, "output1"), autoYes = TRUE)
 pred <- predNew.dr(rst, y~x1+x2, DiskConn, 1000, DiskConn_output)
 head(pred[[1]]$value)
 
 # hdfs backend
 library(Rhipe)
 rhinit()
 seq.file <- list()
 seq.file[[1]] <-  list(2, df[1:500,])
 seq.file[[2]] <-  list(2, df[501:1000,])
 rhwrite(seq.file, file="/tmp/test1", chunk=1, kvpairs=T, verbose=F)
 HDFSconn <- hdfsConn("/tmp/test1", autoYes = TRUE)
 HDFSconn <- ddo(HDFSconn)
 HDFSconn <- updateAttributes(HDFSconn)
 rst<- drml(HDFSconn, y~x1+x2, size =1000, burnin =50, approx_method = "SN") 
 HDFSoutput <- hdfsConn("/tmp/output", autoYes = TRUE)
 pred <- predNew.dr(rst, y~x1+x2, HDFSconn, 1000, HDFSoutput)
 head(pred[[1]]$value)
}
}
