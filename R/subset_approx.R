#' @title Likelihood Modeling for Logistic Regression in subset sense
#' 
#' 
#' @description Model the posterior distribution of parameters in logistic regression by normal distribution or skew normal distribution, where the prior distribution is the uniform distribution.
#' @param x 
#' the model matrix
#' @param y
#' the response variable
#' @param size 
#' the number of MCMC iterations saved (target distribution is the posterior distribution of parameters in the logistic regression)
#' @param burnin 
#' the number of MCMC iterations discarded.
#' @param conf_level 
#' a vector which consists of levels of credible intervals
#' @param approx_method 
#' the method to approximate the posterior distribution such as normal or skew normal, the default one is normal distribution
#' @author Qi Liu
#' @return prob_compare 
#' a dataframe with two columns: approximate probability and true proability. The probability under different credible regions defined by conf_levels is estimated by using Monte Carlo methods
#' @return norm.mean
#' mean parameter of fitted normal distribution
#' @return norm.var
#' variance (covariance) of fitted normal distribution
#' @return sn.xi
#' location parameter of fitted skew normal distribution if approx_method is "SN", Null otherwise
#' @return sn.omega
#' scale parameter of fitted skew normal distributionif approx_method is "SN", Null otherwise
#' @return sn.alpha 
#' shape parameter of fitted skew normal distributionif approx_method is "SN", Null otherwise
#' @examples 
#' 
#'x <- matrix(rnorm(1000*5), ncol=5)
#'y <- rbinom(1000,1,0.5)
#'a <- subset_approx(x,y, size=5000, burnin =500, conf_level = seq(0.05, 0.95, 0.05), approx_method = "Norm")
#'@import sn
#'@import BayesLogit
#'@import MASS
#'@import mvtnorm
#'@import moments
#'@export

subset_approx <- function(x,y, size, burnin, conf_level, approx_method=approx_method){
  p <- dim(x)[2]
  glm_fit <- glm(y~x-1, family=binomial)
  like_mode <- as.numeric(glm_fit$coefficients)
  var_mode <- summary(glm_fit)$cov.unscaled
  
  # log likelihood function for the logistic regression  
  True_loglike <- function(param){
    b <-  apply(x*y, 2, sum) %*% param - sum(log(1+exp(x%*% param)))
    return(b)
  }
  max_true <- as.numeric(True_loglike(like_mode))
  
  #######################################
  # Input: 
  #       mc: samples from approximate distribution
  #       ratio: a vector of ratios to determine the probability within 
  #             contour region
  #Output:
  #       prob: a vector of probability, p_i = p(likelihood(x)/likelihood 
  #             (mode)< ratio_i)
  #######################################
  ratio.fn <- function(mc, ratio ){
    n <- length(ratio)
    prob <- vector(mode = "numeric",length = n)
    for (i in 1:n){
      prob[i] <- mean(apply(mc,1, function(r) ifelse(log(ratio[i]) < True_loglike(r)-max_true ,1,0)))
    }
    prob
  }
  
  ## transform from desired credible regions to ratios
  level2ratio <- function(conflevels){
    n <- length(conflevels)
    rat <- vector(mode = "numeric",length = n)
    for (i in 1:n){
      rat[i] <- exp(-qchisq(conflevels[i], p)/2)
    }
    rat
  }
  
  mc_true <- logit(y, x, samp = size, burn = burnin)$beta
  if (approx_method == "SN"){
    mean.test <- apply(mc_true, 2, mean)
    var.test <- var(mc_true)
    skew.test <- skewness(mc_true)
    d_p <- cp2dp(list(mean = mean.test, var.cov = var.test, gamma1= skew.test), family="SN")
    if(!is.null(d_p)){
      mc_approx <- rmsn(n=size, xi=d_p$beta, Omega=d_p$Omega, alpha=d_p$alpha)
    }else{
      return(print("non-admissible CP"))
    }
  }else{
    mc_approx <- mvrnorm(size, mu=like_mode, Sigma=var_mode)
    d_p <- NULL
  }
  if (!is.null(conf_level)){
  ratio <- level2ratio(conf_level)
  approx_prob <- ratio.fn(mc_approx, ratio)
  true_prob <- ratio.fn(mc_true, ratio)
  rst <- data.frame(approx_prob= approx_prob, true_prob = true_prob)
  }else{
    rst <- NULL
  }
 
  return(list(prob_compare = rst, norm.mean = like_mode, norm.var = var_mode, sn.xi =d_p$beta, sn.omega =d_p$Omega, sn.alpha = d_p$alpha))
}
