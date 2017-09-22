#' @title Model Likelihood of Logistic Regression in Memory
#' 
#' 
#' @description Model the posterior distribution of parameters in logistic regression by normal distribution or skew normal distribution, where the prior distribution is the uniform distribution.
#' @param ddo_object
#' a ddo/ddf object which is obtained by dividing whole data into subsets based on different criteria
#' @param formula 
#' an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are the same to the formula in lm.

#' @param size 
#' the number of MCMC iterations saved (target distribution is the posterior distribution of parameters in the logistic regression)
#' @param burnin 
#' the number of MCMC iterations discarded.
#' @param approx_method 
#' the method to approximate the posterior distribution such as normal or skew normal, the default one is normal distribution
#' @author Qi Liu
#' @return norm.mean
#' mean parameter of recombined fitted normal distribution
#' @return norm.var
#' variance (covariance) of recombined fitted normal distribution
#' @return sn.mod
#' mean parameter of normal approximation to recombined fitted skew normal distribution if approx_method is "SN"
#' @return sn.cov
#' variance (covariance) of normal approximation to recombined fitted skew normal distribution if approx_method is "SN"
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
#'@import sn
#'@import BayesLogit
#'@import MASS
#'@import mvtnorm
#'@import moments
#'@import datadr
#'@export


drml<- function(ddo_object, formula=formula, size, burnin, approx_method=approx_method){
  form_term <- model.frame(formula, data=df)
  col_names <-c("intercept",attr(attr(form_term, "terms"),"term.labels"))
  map_out <- addTransform(ddo_object, function(v){
    LMsubset(formula, data =v, size=size, burnin =burnin, conf_level=NULL, approx_method = approx_method)
  } )
  df_recomb <- recombine(map_out)
  df_val <- lapply(df_recomb, '[[',2)
    rst <- list()
    if (approx_method =="Norm"){
      norm.cov <- solve(Reduce("+", lapply(df_val, function(r){solve(r$norm.var)})))
      rst$norm.mean <- as.vector(Reduce("+", lapply(df_val, function(r){solve(r$norm.var, r$norm.mean)}))%*%norm.cov)
      rst$norm.cov <- norm.cov
      names(rst$norm.mean)<- col_names
      colnames(rst$norm.cov) <- col_names
      rownames(rst$norm.cov) <- col_names
    }else{
     sn.omega.inv <- Reduce("+", lapply(df_val, function(r){solve(r$sn.omega)}))
     sn.xi <- Reduce("+", lapply(df_val, function(r){solve(r$sn.omega, r$sn.xi)}))
     sn.xi <- solve(sn.omega.inv, sn.xi)
     p <- length(sn.xi)
     sn.subxi <- lapply(df_val, '[[',"sn.xi")
     if (!any(is.na(sn.subxi))){
     sn.subxi <- matrix(unlist(sn.subxi), nrow=p)
     alpha.wt <- list()
     alpha.wt <- append(alpha.wt, lapply(df_val, function(r){r$sn.alpha/sqrt(diag(r$sn.omega))}))
     alpha.wt <- matrix(unlist(alpha.wt), nrow = p)
     
     SN.mcmc.out <- optim(
       sn.xi,
       function(tt) {
         first <- (-1/2) * t(tt - sn.xi) %*% sn.omega.inv %*% (tt - sn.xi)
         second <- sum(pnorm(apply(alpha.wt * (tt - sn.subxi), 2, sum), log.p = TRUE))
         (first + second)
       },
       control = list(fnscale = -1),
       method = "L-BFGS-B",
       hessian = TRUE
     )
     rst$sn.mod <- SN.mcmc.out$par
     rst$sn.cov <- solve(-SN.mcmc.out$hessian)
     names(rst$sn.mod)<- col_names
     colnames(rst$sn.cov) <- col_names
     rownames(rst$sn.cov) <-col_names
     }else{
       rst$sn.mod <- NULL
       rst$sn.cov <- NULL
       return(print("Skew normal approximation is not applicable for some subsets"))
     }
     
    }
    
   return(rst)
}





