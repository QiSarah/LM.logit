#' @title Fitted Values at New Data Based on the Fitted Results from Likelihood Modeling
#' 
#' 
#' @description predNew.local is a function to provide the 0.025, 0.5, 0.975 quantiles of the distribution of fitted predict probability based on the fitted density of model parameters.
#' @param fitted_par
#' object returned by drml function
#' @param formula 
#' an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are the same to the formula in lm.

#' @param size 
#' the number of samples drawn from the distribution of model parameters
#' @return value
#' 0.025, 0.5, 0.975 quantiles of the distribution of fitted predict probability
#' @examples
#' set.seed(100)
#' library(datadr)
#' x <- matrix(rnorm(1000*5), ncol=5)
#' ttheta <- rep(1,5)
#' y <- rbinom(1000, 1, 1 / (1 + exp(- x %*% ttheta)))
#' df <- cbind(y,x)
#' df <- as.data.frame(df)
#' names(df) <- c("y", "x1","x2","x3","x4","x5")
#' df_ddf <- ddf(df)
#' df_div <- divide(df_ddf, by =rrDiv(500))
#' rst<- drml(df_div, y~x1+x2, size =1000, burnin =50, approx_method = "SN") 
#' pred <- predNew.local(rst, y~x1+x2, df, 1000)
#'@import MASS
#'@import mvtnorm
#'@import moments
#'@import datadr
#'@export
predNew.local <- function(fitted_par,formula, data, size=1000){
  samp <- mvrnorm(size, mu=fitted_par[[1]], Sigma=fitted_par[[2]])
  x  <- model.frame(formula=formula, data=data)
  sigmod <- function(tt){1/(1+exp(-tt))}
  sigmod_vec <- Vectorize(sigmod)
  y <- apply(x, 1, function(r){quantile(sigmod_vec(r%*%t(samp)), probs=c(0.025,0.5, 0.975))})
    return(t(y))
}