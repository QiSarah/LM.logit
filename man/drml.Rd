% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drml.R
\name{drml}
\alias{drml}
\title{Model Likelihood of Logistic Regression in Divide and Recombine Framework}
\usage{
drml(ddo_object, formula = formula, size, burnin,
  approx_method = approx_method)
}
\arguments{
\item{ddo_object}{a ddo/ddf object (in memory, ) which is obtained by dividing whole data into subsets}

\item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are the same to the formula in lm.}

\item{size}{the number of MCMC iterations saved (target distribution is the posterior distribution of parameters in the logistic regression)}

\item{burnin}{the number of MCMC iterations discarded.}

\item{approx_method}{the method to approximate the posterior distribution such as normal or skew normal, the default one is normal distribution}
}
\value{
\item{norm.mean}{mean parameter of recombined fitted normal distribution}

\item{norm.var}{variance (covariance) of recombined fitted normal distribution}

\item{sn.mod}{mean parameter of normal approximation to recombined fitted skew normal distribution if approx_method is "SN"}

\item{sn.cov}{variance (covariance) of normal approximation to recombined fitted skew normal distribution if approx_method is "SN"}
}
\description{
Model the posterior distribution of parameters in logistic regression by normal distribution or skew normal distribution, where the prior distribution is the uniform distribution.
}
\seealso{
\code{\link{divide}}, \code{\link{recombine}}
}
\examples{
set.seed(100)
library(datadr)
x <- matrix(rnorm(1000*5), ncol=5)
ttheta <- rep(1,5)
y <- rbinom(1000, 1, 1 / (1 + exp(- x \%*\% ttheta)))
df <- cbind(y,x)
df <- as.data.frame(df)
names(df) <- c("y", "x1","x2","x3","x4","x5")
df_ddf <- ddf(df)
df_div <- divide(df_ddf, by =rrDiv(500))
rst<- drml(df_div, y~x1+x2, size =1000, burnin =50, approx_method = "Norm") 
}
\author{
Qi Liu
}
