# AdaReg
AdaReg is a R package to robustly estimate linear regression coefficients and detect outliers, analyzed in the paper "AdaReg: Data Adaptive Robust Estimation in Linear Regression with Application in GTEx Gene Expressions". 

Please contact Meng Wang by email <mengw1@stanford.edu> for questions. 

## Installation
`library(devtools)`

`install_github("mwgrassgreen/AdaReg")`

## Usage
`library(AdaReg)`

`adareg.result = AdaReg(design.mx, response.y, gam.seq = seq(0, 3, by=.1), var.gp.id=NULL, tol=10^(-4), step=50)`

`adareg.coef = adareg.result$beta.rob.fit`

`adareg.var = adareg.result$var.sig.gp.fit`

`adareg.res = adareg.result$x.res`

`res.fit = adareg.result$res.info`

### Input
   design.mx --- the design matrix in linear regression
   
   rsponse.y --- the response variable in linear regression
   
   gam.seq --- a sequence of gamma's (default: seq(0, 3, by=.1))
   
   var.gp.id --- a vector of groups indicating the samples in the same variance group (default: all the samples have the same variance)
   
   tol --- tolerance for interations (default: 10^(-4))
   
   step --- step limit (default: 50)
        
### Output
   adareg.coef --- the estimated coefficients
   
   adareg.var --- the estimated sample variance
   
   adareg.res --- a vector of residuals
   
   res.fit --- fitting info for the residuals
   
## Example
To simulate a linear regression dataset from sim.dat.fn() (details to see the paper)

`dat = sim.dat.fn(n=1000, p=20, pi1=0.1, mu1=3, x.mu=0, x.sd=10, e.mu=0, e.sd=1)`

`y.0 = dat$y` (response variable)

`X.0 = dat$X` (design matrix)

`out.ind = dat$ind`(outlier indices)
 
 To apply AdaReg
 
`adareg.result = AdaReg(X.0, y.0, gam.seq = seq(0, 3, by=.1))`

`adareg.coef = adareg.result$beta.rob.fit`

`adareg.res = adareg.result$x.res`

To plot the histogram of the residuals

`hist(adareg.res, freq=FALSE, main="Histogram of residuals")`

`curve(res.fit["pi0.hat"] * dnorm(x, res.fit["mu0.hat"], res.fit["sd0.hat"]), min(adareg.res, na.rm=TRUE), max(adareg.res, na.rm=TRUE), add=TRUE, col="red", lwd=2)` (fitted Gaussian desnity curve)

`rug(adareg.res[out.ind], col="blue", lwd=2)` (true outliers)
