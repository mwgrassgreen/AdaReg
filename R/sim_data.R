#' @title sim.dat.fn
#' @description To simulation linear regression dataset
#' @author Meng Wang
#' \email{mengw1@stanford.edu}
#' @param n Sample size.
#' @param p Number of covariates (including intercept).
#' @param pi1 Proportion of outliers in the response variable.
#' @param mu1 Expection of outliers.
#' @param x.mu Expection of design matrix X (default: 0).
#' @param x.sd Standard deviation of design matrix X (default: 10).
#' @param e.mu Expection of noise (default: 0).
#' @param e.sd Standard deviation of noise (default: 1).
#' @return Simulated data and its parameters.
#' @export
#' 
sim.dat.fn = function (n, p, pi1, mu1, x.mu=0, x.sd=10, e.mu=0, e.sd=1) {
	    beta.00 = c(1, rep(1, p-1))
		X = cbind(rep(1,n), matrix(rnorm(n*(p-1), x.mu, x.sd), n, p-1))
		e.1 = rnorm(n, e.mu, e.sd)
		ind = sample(1:n, n*pi1)
		e.1[ind] = rnorm(n*pi1, mu1,1)
		y = X %*% beta.00 + e.1
		
		return(list(y=y, X=X, beta.00=beta.00, ind=ind))
}
