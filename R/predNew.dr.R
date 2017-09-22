#' @title Fitted Values at distributed Datasets Based on the Fitted Results from Likelihood Modeling
#' 
#' 
#' @description predNew.dr is a function to provide the 0.025, 0.5, 0.975 quantiles of the distribution of fitted predict probability based on the fitted density of model parameters.
#' @param fitted_par
#' object returned by drml function
#' @param formula 
#' an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are the same to the formula in lm.

#' @param size 
#' the number of samples drawn from the distribution of model parameters
#' @param output
#' a "kvConnection" object indicating where the output data should reside (see localDiskConn, hdfsConn). 
#' @return value
#' 0.025, 0.5, 0.975 quantiles of the distribution of fitted predict probability
#' @examples
#' \dontrun{
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
#' HDFSconn <- hdfsConn("/tmp/KV", autoYes = TRUE)
#' addData(HDFSconn, df_div)
#' HDFSconn <- ddo(HDFSconn)
#' HDFSconn <- updateAttributes(HDFSconn)
#' rst<- drml(HDFSconn, y~x1+x2, size =1000, burnin =50, approx_method = "SN") 
#' output <- hdfsConn("/tmp/output", autoYes = TRUE)
#' pred <- predNew.dr(rst, y~x1+x2, HDFSconn, 1000, output)
#' get_pred <- updateAttributes(pred)
#' head(pred[[1]]$value)
#' }
#' @seealso
#' \code{\link{predNew.local}}\code{\link{datadr}}
#'@import MASS
#'@import mvtnorm
#'@import moments
#'@import datadr
#'@export
predNew.dr <- function(fitted_par,formula, ddo_object, size=1000, output){
  samp <- mvrnorm(size, mu=fitted_par[[1]], Sigma=fitted_par[[2]])
  sigmod <- function(tt){1/(1+exp(-tt))}
  sigmod_vec <- Vectorize(sigmod)
  map_out <- addTransform(ddo_object, function(v){
    x  <- model.frame(formula=formula, data=v)
    y <- apply(x, 1, function(r){quantile(sigmod_vec(r%*%t(samp)), probs=c(0.025,0.5, 0.975))})
    return(t(y))
  } )
  rst <- recombine(map_out, output=output)
  updateAttributes(rst)
  
}