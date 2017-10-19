# LM.logit: Model Likelihood of Logistic Regression in Divide and Recombine

[![Build Status](https://github.com/QiSarah/LM.logit.git)](https://github.com/QiSarah/LM.logit.git)

In Divide & Recombine (D&R), big data are divided into subsets, each analytic method is applied to subsets, and the outputs are recombined. The package is the implementation of an innovative D&R procedure to compute likelihood functions of data-model parameters for big data. The likelihood-model is a parametric probability density function of the data-model parameters. The density parameters are estimated by fitting the density to MCMC draws from each subset data-model likelihood function, and then the fitted densities are recombined.

In summary, this is an R package that provides a way to fit logistic regression to big data based 
on likelihood modeling routine under the divide and recombined framework.


## Installation

```r
# from github
devtools::install_github("QiSarah/LM.logit")
```

## Fitting

The main function to call in the package is `drml` function, which can
be used to estimate the parameters of fitted recombined likelihood function of data-model parameters in the logistic regression on the dataset either in the memory or on the HDFS.

```r
 library(datadr)
 # simulate data
 x <- matrix(rnorm(1000*5), ncol=5)
 ttheta <- rep(1,5)
 y <- rbinom(1000, 1, 1 / (1 + exp(- x %*% ttheta)))
 df <- cbind(y,x)
 df <- as.data.frame(df)
 names(df) <- c("y", "x1","x2","x3","x4","x5")
 # divide the data into multiple subsets
 df_ddf <- ddf(df)
 df_div <- divide(df_ddf, by =rrDiv(500))
 # fit data to logistic regression
 rst<- drml(df_div, y~x1+x2, size =1000, burnin =50, approx_method = "SN") 
 ```
## Prediction 

The function for prediction is a`predNew.local` function, which is used to provide the 0.025, 0.5, 0.975 quantiles of the distribution of fitted predict probability using the fitted density of model parameters. The prediction for new dataset
will be conducted in local memory.
```r
df_div <- divide(df_ddf, by =rrDiv(500))
rst<- drml(df_div, y~x1+x2, size =1000, burnin =50, approx_method = "SN") 
pred <- predNew.local(rst, y~x1+x2, df, 1000)
```
The function for prediction is `predNew.dr` function, which is similar to `predNew.local`. The prediction for new dataset
will be conducted on HDFS using MapReduce or on local disk.
```r
HDFSconn <- hdfsConn("/tmp/KV", autoYes = TRUE)
addData(HDFSconn, df_div)
HDFSconn <- ddo(HDFSconn)
HDFSconn <- updateAttributes(HDFSconn)
rst<- drml(HDFSconn, y~x1+x2, size =1000, burnin =50, approx_method = "SN") 
output <- hdfsConn("/tmp/output", autoYes = TRUE)
pred <- predNew.dr(rst, y~x1+x2, HDFSconn, 1000, output)
get_pred <- updateAttributes(pred)
head(pred[[1]]$value)
```

## Subset Approximation

For each subset, the density parameters are estimated by fitting the density to MCMC draws from logistic regression likelihood function. Two likelihood-models have been considered, skew-normal and normal. The function `LMsubset` which is a wrapper of `subset_approx` is used to estimate parameters in one of two likelihood-model (skew-normal or normal) and also provide a list of credible intervals to check the performance of likelihood modeling.
```r
x <- matrix(rnorm(1000*5), ncol=5)
y <- rbinom(1000,1,0.5)
rst <- subset_approx(x,y, size=5000, burnin =500, conf_level = seq(0.05, 0.95, 0.05), approx_method = "Norm")
```
## Acknowledgment

LM.logit development depends on: 
- datadr package by Ryan Hafen
- Rhipe package by Saptarshi Guha
- sn package by Adelchi Azzalini
- BayesLogit package by Nicholas G. Polson, James G. Scott, and Jesse Windle
- MASS package by Brian Ripley, Bill Venables, Douglas M. Bates, Kurt Hornik, Albrecht Gebhardt, David Firth
- moments package by Lukasz Komsta
- mvrnorm function from mvtnorm package in base R by R core Team
