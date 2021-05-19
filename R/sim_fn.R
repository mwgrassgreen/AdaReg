#' @title sim_fn.R
#' @description Functions for the simulation studies
#' @author Meng Wang
#' \email{mengw1@stanford.edu}



# ------------------------------------------------------------------------------
# data generation function
# mixture distrition: (1 - pi1) * pop.distr + pi1 * out.distr
# to specify pop.distr, out.distr to be 'norm', 't', 'gamma'

dat.gen.fn = function (n, pi1, s, pop.distr='norm', out.distr='norm', mu.0=NULL, sd.0=NULL, t.shift.0=NULL, t.df.0=NULL, t.df.1=NULL, shape.0=NULL, rate.0=NULL, meanlog.0=NULL, sdlog.0=NULL) {

	if (pop.distr == 'norm') {
	   x = rnorm(n, mu.0, sd.0)
	   sd.1 = sd.0		
	   mean.1 = mu.0 + s * sd.1
	}
	if (pop.distr == 't') {
	   x = t.shift.0 + rt(n, df=t.df.0)
	   # to choose t.df.0 > 2
	   sd.1 = sqrt(t.df.0/(t.df.0 - 2))
	   mean.1 = t.shift.0 + s * sd.1
	}
	if (pop.distr == 'gamma') {
	   x = rgamma(n, shape=shape.0, rate=rate.0)
	   sd.1 = sqrt(shape.0/rate.0^2)
	   mean.1 = shape.0/rate.0 + s * sd.1
	}
	if (pop.distr == 'lognorm') {
	   x = rlnorm(n, meanlog=meanlog.0, sdlog=sdlog.0)
	   sd.1 = sqrt((exp(sdlog.0^2) - 1) * exp(sdlog.0^2 + 2*meanlog.0))
	   mean.1 = exp(meanlog.0 + sdlog.0^2/2) + s * sd.1
	}
	if (! pop.distr %in% c('norm', 't', 'gamma', 'lognorm')) {
       stop('to select pop.distr to be norm, t, gamma, lognorm')
    }
	
	if (pi1 > 0) {
		n1 = round(n*pi1)
		o.ind = sample(1:n, n1)
		if (out.distr == 'norm') {
	       x[o.ind] = rnorm(n1, mean.1, sd.1)		
        }
        if (out.distr == 't') {
	       x[o.ind] = mean.1 + rt(n1, df=t.df.1)		
        }    
        if (! out.distr %in% c('norm', 't')) {
        	stop('to select out.distr to be norm or t')
        }
        
	} else {
		o.ind = integer(0)
	}
    return(list(x=x, ind=o.ind))
}



# ------------------------------------------------------------------------------
# to estimate theta	    
theta.trend.est.fn = function (n, step=50, BB=50, gam.seq=seq(0, 3, by=0.5), pop.distr, out.distr, norm.mu0=NULL, norm.sd0=NULL, t.shift0=NULL, t.df0=NULL, t.df1=NULL, gamma.shape0=NULL, gamma.rate0=NULL, lognorm.meanlog0=NULL, lognorm.sdlog0=NULL) {
	
	# parameters for the outlier distribution
	par.01.mx = rbind(c(0, 0), expand.grid(c(0.1, 0.3), c(1,3,5)) )
	colnames(par.01.mx) = c("pi1", "s")
	par.01.mx = par.01.mx[order(par.01.mx[,1]), ]
	dim(par.01.mx)
	#par.01.mx
    par.01.mx.comb = cbind(par.01.mx, matrix(NA, nrow(par.01.mx), 2))
    colnames(par.01.mx.comb)[3:4] = c('mean.1', 'sd.1')
    
	if (pop.distr == 'norm') {
		mean.00 = norm.mu0
		sd.00 = norm.sd0
	}
	if (pop.distr == 't') {
		mean.00 = t.shift0
		sd.00 = sqrt(t.df0/(t.df0-2))
	}
	if (pop.distr == 'gamma') {
		mean.00 = gamma.shape0/gamma.rate0
		sd.00 = sqrt(gamma.shape0/gamma.rate0^2)
	}
	if (pop.distr == 'lognorm') {
	   mean.00 = exp(lognorm.meanlog0 + lognorm.sdlog0^2/2)
	   sd.00 = sqrt((exp(lognorm.sdlog0^2) - 1) * exp(lognorm.sdlog0^2 + 2*lognorm.meanlog0))
	}
    par.01.mx.comb[1, 3:4] = c(mean.00, sd.00) 
			
    # estimators under various mixture model
	fit.crt.array = array(NA, dim=c(length(gam.seq),BB, nrow(par.01.mx)), dimnames=list(gam=gam.seq, B=1:BB, par=paste("(pi1=", par.01.mx[,1], ",mean1=",par.01.mx[,2], ")", sep="")) )
	fit.par.mu.array = fit.par.sd.array = fit.par.pi0.array =fit.par.array = fit.crt.array

	for (i in 1:nrow(par.01.mx)) {
		print(i)
		pi1 = par.01.mx[i, 1]
		s = par.01.mx[i, 2]
		
		if (i > 1) {
			if (pop.distr == 'norm') {
			   sd.1 = norm.sd0	
			   mean.1 = norm.mu0 + s * norm.sd0
			}
			if (pop.distr == 't') {
			   # to choose t.df.0 > 2
			   sd.1 = sqrt(t.df0/(t.df0 - 2))
			   mean.1 = t.shift0 + s * sd.1
			}
			if (pop.distr == 'gamma') {
			   sd.1 = sqrt(gamma.shape0/gamma.rate0^2)
			   mean.1 = gamma.shape0/gamma.rate0 + s * sd.1
			}
			if (pop.distr == 'lognorm') {
			   sd.1 = sqrt((exp(lognorm.sdlog0^2) - 1) * exp(lognorm.sdlog0^2 + 2*lognorm.meanlog0))
			   mean.1 = exp(lognorm.meanlog0 + lognorm.sdlog0^2/2) + s * sd.1
			}
			if (! pop.distr %in% c('norm', 't', 'gamma', 'lognorm')) {
		       stop('to select pop.distr to be norm, t, gamma, lognorm')
		    }
		    par.01.mx.comb[i, 3:4] = c(mean.1, sd.1) 
		}
		
		for (b in 1:BB) {		
			x = dat.gen.fn(n, pi1, s, pop.distr, out.distr, mu.0=norm.mu0, sd.0=norm.sd0, t.shift.0=t.shift0, t.df.0=t.df0, t.df.1=t.df1, shape.0=gamma.shape0, rate.0=gamma.rate0, meanlog.0=lognorm.meanlog0, sdlog.0=lognorm.sdlog0)[['x']]
	        		
			fit = adapt.gam.rob.fit.fn(x, gam.seq, step)
			#names(fit)
			#fit$est.hat
			fit.crt = abs(fit$para.hat.mx[,"efdr0.hat"]-1)
			fit.crt.array[names(fit.crt),b, i] = fit.crt
			
			fit.par = (fit$para.hat.mx[,"mu0.hat"] - mean.00)^2 + (fit$para.hat.mx[,"sd0.hat"] - sd.00)^2
			fit.par.array[names(fit.par),b, i] = fit.par
			
			fit.par.mu = fit$para.hat.mx[,"mu0.hat"] 
			fit.par.mu.array[names(fit.par.mu),b, i] = fit.par.mu
			
			fit.par.sd = fit$para.hat.mx[,"sd0.hat"]
			fit.par.sd.array[names(fit.par.sd),b, i] = fit.par.sd
			
			fit.par.pi0 = pmin(1,fit$para.hat.mx[,"pi0.hat"])
			names(fit.par.pi0) = names(fit$para.hat.mx[,"pi0.hat"])
			fit.par.pi0.array[names(fit.par.pi0),b, i] = fit.par.pi0			
	     }
    }
    
	return(list(par.01.mx.comb=par.01.mx.comb, fit.crt.array=fit.crt.array, fit.par.array=fit.par.array, fit.par.mu.array=fit.par.mu.array, fit.par.sd.array=fit.par.sd.array, fit.par.pi0.array=fit.par.pi0.array))
	
}


# to calculate theoretical E(fdr) under the null
efdr.null.int.fn = function(par.01.mx.comb, pop.distr, out.distr, mu.0=NULL, sd.0=NULL, t.shift.0=NULL, t.df.0=NULL, t.df.1=NULL, shape.0=NULL, rate.0=NULL, meanlog.0=NULL, sdlog.0=NULL){
	efdr.null.int = matrix(NA, nrow(par.01.mx.comb), ncol(par.01.mx.comb)+4)
	for (i in 1:nrow(par.01.mx.comb)) {
	   # theoretical E_0(fdr) function
       theo.efdr.fn = function(x, s=par.01.mx[i,2], pi1=par.01.mx.comb[i,1], mean.1 = par.01.mx.comb[i, 'mean.1'], sd.1 = par.01.mx.comb[i, 'sd.1']) {
		 	   if (pop.distr == 'norm') {
		    		den.0 = dnorm(x, mu.0, sd.0)
		    	}
		    	if (pop.distr == 't') {
		    		den.0 = dt(x-t.shift.0, t.df.0)
		    	}
		    	if (pop.distr == 'gamma') {
		    	   den.0 = dgamma(x, shape.0, rate=rate.0)
		    	}
		        if (pop.distr == 'lognorm') {
		    	   den.0 = dlnorm(x, meanlog.0, sdlog.0)
		    	}
		    		    	
		    	if (out.distr == 'norm') {	                    
		    		den.1 = dnorm(x, mean.1, sd.1)
		    	}
		    	if (out.distr == 't') {
		    		den.1 = dt(x-mean.1, t.df.1)
		    	}	    			    	
		    (1-pi1) * (den.0)^2/((1-pi1) * den.0 + pi1 * den.1)		   
	     }
        
        if (pop.distr == 'norm') {
        	low.lim = qnorm(0.0001, mu.0, sd.0)
    	}
    	if (pop.distr == 't') {
    		low.lim = t.shift.0 + qt(0.0001, t.df.0) 
    	}
    	if (pop.distr == 'gamma') {
    	   low.lim = qgamma(0.0001, shape=shape.0, rate=rate.0)
    	}
    	if (pop.distr == 'lognorm') {
    	   low.lim = qlnorm(0.0001, meanlog.0, sdlog.0)
    	}
    	mean.1 = par.01.mx.comb[i, 'mean.1']
	    sd.1 = par.01.mx.comb[i, 'sd.1']
  		if (out.distr == 'norm') {
  			upp.lim = qnorm(1-0.0001, mean.1, sd.1)
    	}
    	if (out.distr == 't') {
    		upp.lim = mean.1 + qt(1-0.0001, t.df.1)
    	}
		    	
		efdr.int = integrate(theo.efdr.fn, low.lim, upp.lim)
		efdr.null.int[i,] = unlist(c(c(par.01.mx.comb[i, ]), efdr.int$value, low.lim, upp.lim, efdr.int$abs.error))
	}
	colnames(efdr.null.int) = c(colnames(par.01.mx.comb), "efdr.null.int", 'low.limit', 'upp.limit', "int.err")
	return(efdr.null.int)

}

# ------------------------------------------------------------------------------
# to plot criterion trend
plot.crt.trend.fn = function (gam.seq, result, fdr.null.int, mix.model.nm, file.nm=file.trend.nm) {
	# to compare the trend of E_0(|fdr - 1|) to the trend of SE(theta_hat, theta)
	fit.par.array = result[['fit.par.array']]
	par.01.mx = result[['par.01.mx.comb']]
	fit.crt.array = result[['fit.crt.array']]
	
	pdf(file=paste(file.nm,'_part1.pdf',sep=''), width=15, height=7)
	par(mfrow=c(2,4), oma=c(0,0,3,0))
	for (i in 1:4) {
		matplot(gam.seq, fit.par.array[,,i], type="b", col=1, pch=1, main=paste("(pi1=", par.01.mx[,1], ", s=",par.01.mx[,2], ")", sep="")[i], xlab="gamma", ylab="SE(theta_hat)", cex.lab=1.5, cex.main=1.5)
		lines(gam.seq, rowMeans(fit.par.array[,,i], na.rm=TRUE), col="red", lwd=2)
	}
	for (i in 1:4) {
		matplot(gam.seq, fit.crt.array[,,i], type="b", col=1, pch=1, main=paste("(pi1=", par.01.mx[,1], ", s=",par.01.mx[,2], ") \n E_0(fdr) = ", round(fdr.null.int[i, 'efdr.null.int'], 2), sep="")[i], xlab="gamma", ylab="|E_0(fdr).hat - 1|", cex.lab=1.5, cex.main=1.5)
		lines(gam.seq, rowMeans(fit.crt.array[,,i], na.rm=TRUE), col="red", lwd=2)
	}
	title(paste("MSE trend vs. selection trend across gamma's under", mix.model.nm,  "and n = ",n ), outer=TRUE, cex.main=2)
	
	dev.off()
	
	
	pdf(file=paste(file.nm,'_part2.pdf',sep=''), width=15, height=7)
	par(mfrow=c(2,3), oma=c(0,0,3,0))
	for (i in 5:7) {
		matplot(gam.seq, fit.par.array[,,i], type="b", col=1, pch=1, main=paste("(pi1=", par.01.mx[,1], ", s=",par.01.mx[,2], ")", sep="")[i], xlab="gamma", ylab="SE(theta_hat)", cex.lab=1.5, cex.main=1.5)
		lines(gam.seq, rowMeans(fit.par.array[,,i], na.rm=TRUE), col="red", lwd=2)
	}
	for (i in 5:7) {
		matplot(gam.seq, fit.crt.array[,,i], type="b", col=1, pch=1, main=paste("(pi1=", par.01.mx[,1], ", s=",par.01.mx[,2], ") \n E_0(fdr) = ", round(fdr.null.int[i, 'efdr.null.int'], 2), sep="")[i], xlab="gamma", ylab="|E_0(fdr).hat - 1|", cex.lab=1.5, cex.main=1.5)
		lines(gam.seq, rowMeans(fit.crt.array[,,i], na.rm=TRUE), col="red", lwd=2)
	}
	title(paste("MSE trend vs. selection trend across gamma's under", mix.model.nm,  "and n = ",n ), outer=TRUE, cex.main=2)
	
	dev.off()

}

# ------------------------------------------------------------------------------
# to plot theta.hat in various situations
plot.theta.fn = function (gam.seq, result, mix.model.nm, file.nm=file.theta.nm) {
	fit.par.mu.array = result[['fit.par.mu.array']]
	fit.par.sd.array = result[['fit.par.sd.array']]
	fit.par.pi0.array = result[['fit.par.pi0.array']]
	par.01.mx = result[['par.01.mx.comb']]
	mu.00 = par.01.mx[1, 'mean.1']
	sd.00 = par.01.mx[1, 'sd.1']
	
	pdf(file=paste(file.nm,'_part1.pdf',sep=''), width=15, height=7)
	par(mfrow=c(3,4), oma=c(0,0,3,0))
	for (i in 1:4) {
		matplot(gam.seq, fit.par.mu.array[,,i], type="b", col=1, pch=1, main=paste("(pi1=", par.01.mx[,1], ", s=",par.01.mx[,2], ")", sep="")[i], xlab="gamma", ylab="mu0_hat", cex.lab=1.5, cex.main=1.5, ylim=c(min(fit.par.mu.array[,,i], mu.00), max(fit.par.mu.array[,,i], mu.00)))
		abline(h=mu.00, lwd=2, col="blue")
		lines(gam.seq, rowMeans(fit.par.mu.array[,,i], na.rm=TRUE), col="red", lwd=2)
	}
	for (i in 1:4) {
		matplot(gam.seq, fit.par.sd.array[,,i], type="b", col=1, pch=1, main=paste("(pi1=", par.01.mx[,1], ", s=",par.01.mx[,2], ")", sep="")[i], xlab="gamma", ylab="sd0_hat", cex.lab=1.5, cex.main=1.5, ylim=c(min(fit.par.sd.array[,,i], sd.00), max(fit.par.sd.array[,,i], sd.00)))
		abline(h=sd.00, lwd=2, col="blue")
		lines(gam.seq, rowMeans(fit.par.sd.array[,,i], na.rm=TRUE), col="red", lwd=2)
	}
	for (i in 1:4) {
		matplot(gam.seq, 1 - fit.par.pi0.array[,,i], type="b", col=1, pch=1, main=paste("(pi1=", par.01.mx[,1], ", s=",par.01.mx[,2], ")", sep="")[i], xlab="gamma", ylab="pi1_hat", cex.lab=1.5, cex.main=1.5, ylim=c(0, max(1 - fit.par.pi0.array[,,i], par.01.mx[i,1])))
		abline(h=par.01.mx[i,1], lwd=2, col="blue")
		lines(gam.seq, rowMeans(1 - fit.par.pi0.array[,,i], na.rm=TRUE), col="red", lwd=2)
	}
	title(paste("trend of theta_hat across gamma's under", mix.model.nm,  "and n = ",n ), outer=TRUE, cex.main=2)
	dev.off()
	
	pdf(file=paste(file.nm,'_part2.pdf',sep=''), width=15, height=7)
	par(mfrow=c(3,3), oma=c(0,0,3,0))
	for (i in 5:7) {
		matplot(gam.seq, fit.par.mu.array[,,i], type="b", col=1, pch=1, main=paste("(pi1=", par.01.mx[,1], ", s=",par.01.mx[,2], ")", sep="")[i], xlab="gamma", ylab="mu0_hat", cex.lab=1.5, cex.main=1.5, ylim=c(min(fit.par.mu.array[,,i], mu.00), max(fit.par.mu.array[,,i], mu.00)))
		abline(h=mu.00, lwd=2, col="blue")
		lines(gam.seq, rowMeans(fit.par.mu.array[,,i], na.rm=TRUE), col="red", lwd=2)
	}
	for (i in 5:7) {
		matplot(gam.seq, fit.par.sd.array[,,i], type="b", col=1, pch=1, main=paste("(pi1=", par.01.mx[,1], ", s=",par.01.mx[,2], ")", sep="")[i], xlab="gamma", ylab="sd0_hat", cex.lab=1.5, cex.main=1.5, ylim=c(min(fit.par.sd.array[,,i], sd.00), max(fit.par.sd.array[,,i], sd.00)))
		abline(h=sd.00, lwd=2, col="blue")
		lines(gam.seq, rowMeans(fit.par.sd.array[,,i], na.rm=TRUE), col="red", lwd=2)
	}
	for (i in 5:7) {
		matplot(gam.seq, 1 - fit.par.pi0.array[,,i], type="b", col=1, pch=1, main=paste("(pi1=", par.01.mx[,1], ", s=",par.01.mx[,2], ")", sep="")[i], xlab="gamma", ylab="pi1_hat", cex.lab=1.5, cex.main=1.5, ylim=c(0, max(1 - fit.par.pi0.array[,,i], par.01.mx[i,1])))
		abline(h=par.01.mx[i,1], lwd=2, col="blue")
		lines(gam.seq, rowMeans(1 - fit.par.pi0.array[,,i], na.rm=TRUE), col="red", lwd=2)
	}
	title(paste("trend of theta_hat across gamma's under", mix.model.nm,  "and n = ",n ), outer=TRUE, cex.main=2)
	dev.off()
	
}


# ------------------------------------------------------------------------------
# data generation function in the regression model
dat.gen.reg.fn = function (n, p, pi1, s, pop.distr='norm', out.distr='norm', mu.0=0, sd.0=NULL, t.shift.0=0, t.df.0=NULL, t.df.1=NULL, shape.0=NULL, rate.0=1, meanlog.0=0, sdlog.0=NULL) {
	    x.mu = 0
	    x.sd = 10
	    beta.00 = c(1, rep(1, p-1))
		X = cbind(rep(1,n), matrix(rnorm(n*(p-1), x.mu, x.sd), n, p-1))
		dat.e = dat.gen.fn (n, pi1, s, pop.distr, out.distr, mu.0, sd.0, t.shift.0, t.df.0, t.df.1, shape.0, rate.0, meanlog.0, sdlog.0)		 
		e.1 = dat.e[['x']]
		if (pop.distr == 'norm') {
			e.1 = e.1 - mu.0
		}
		if (pop.distr == 't') {
			e.1 = e.1 - t.shift.0
		}
		if (pop.distr == 'gamma') {
			e.1 = e.1 - shape.0/rate.0
		}
		if (pop.distr == 'lognorm') {
			 e.1 = e.1 - exp(meanlog.0 + sdlog.0^2/2)
		}
		if (! pop.distr %in% c('norm', 't', 'gamma', 'lognorm')) {
	       stop('to select pop.distr to be norm, t, gamma, lognorm')
	    }
    	ind = dat.e[['ind']]
		y = X %*% beta.00 + e.1
		
		return(list(y=y, X=X, beta.00=beta.00, ind=ind))
}

# ------------------------------------------------------------------------------
# estimation comparison in the regression model

est.comp.in.reg.model.fn = function (n, p, pop.distr='norm', out.distr='norm', mu.0=0, sd.0=1, t.shift.0=0, t.df.0=NULL, t.df.1=NULL, shape.0=NULL, rate.0=1, meanlog.0=0, sdlog.0=NULL, BB=100, method.ind=0) {
	
	par.01.mx = rbind(c(0, 0), expand.grid(seq(0.1, 0.4, by=0.1), 1:10) )
	colnames(par.01.mx) = c("pi1", "s")
	par.01.mx = par.01.mx[order(par.01.mx[,1]), ]
	dim(par.01.mx)
	
	method.nm = c("ols",  "w.gam=0.5", "w.gam=1",   "w.gam=2",   "w.gam=3",  "w.ada.gam",  "w.huber",   "w.hampel",  "w.bisq", "S" ,"lms"  ,    "lts"  )
	coef.array = array(NA, dim=c(length(method.nm), p, BB, nrow(par.01.mx)), dimnames=list(method=method.nm, fit.coef=paste("coef", 0:(p-1), sep="."), B=1:BB, par=paste("(pi1=", par.01.mx[,1], ", s=",par.01.mx[,2], ")", sep="")) )
	sig.array = array(NA, dim=c(length(method.nm), BB, nrow(par.01.mx)), dimnames=list(method=method.nm,  B=1:BB, par=paste("(pi1=", par.01.mx[,1], ", s=",par.01.mx[,2], ")", sep="")) )
	tpr.array = fpr.array = sig.array
	res.array = array(NA, dim=c(6, 3, BB, nrow(par.01.mx)), dimnames=list(method=method.nm[1:6],   fit.res=c("res.mean.obv", "res.sd.obv", "res.pop.prp.obv"), B=1:BB, par=paste("(pi1=", par.01.mx[,1], ", s=",par.01.mx[,2], ")", sep="")  ) )

	
	rmse.coef.mx = array(NA, dim=c(length(method.nm), nrow(par.01.mx)), dimnames=list(method=method.nm, par=paste("(pi1=", par.01.mx[,1], ", s=",par.01.mx[,2], ")", sep="")) )
	rmse.sig.mx = avg.tpr.mx = avg.fpr.mx = rmse.coef.mx

    if (pop.distr == 'norm') {
    	sd.00 = sd.0
    } 
    if (pop.distr == 't') {
       sd.00 = sqrt(t.df.0/(t.df.0 - 1))
    } 
    if (pop.distr == 'gamma') {
	   sd.00 = sqrt(shape0/rate0^2)
	}
	if (pop.distr == 'lognorm') {
	   sd.00 = sqrt((exp(sdlog.0^2) - 1) * exp(sdlog.0^2 + 2*meanlog.0))
	}
	if (! pop.distr %in% c('norm', 't', 'gamma', 'lognorm')) {
		stop('to select pop.distr to be norm, t, gamma, lognorm')
    }

    
    for (i in 1:nrow(par.01.mx)) {
		print(i)
		pi1 = par.01.mx[i, 1]
		s = par.01.mx[i, 2]
		
		for (b in 1:BB) {
	
			dat = dat.gen.reg.fn (n, p, pi1, s, pop.distr, out.distr, mu.0, sd.0, t.shift.0, t.df.0, t.df.1, shape.0, rate.0, meanlog.0, sdlog.0) 
			y = dat[['y']]
			X = dat[['X']]
			beta.00 = dat[['beta.00']]
			o.ind = dat[['ind']]
			
			if (method.ind == 0) {
				### OLS
				fit.0 = lm(y ~ X - 1)
				#names(fit.0)
				#plot(fit.0)
				fit.0.coef = fit.0$coef
				fit.0.sig = sqrt(sum((fit.0$res)^2)/(n - p))
				fit.0.z = (fit.0$res)/fit.0.sig
				fit.0.out = (1:n)[abs(fit.0.z) > 2.5]
				# length(fit.0.out)
				fit.0.tpr = length(intersect(fit.0.out, o.ind))/length(o.ind)
				fit.0.fpr = length(intersect(fit.0.out, (1:n)[-o.ind]))/(n - length(o.ind))
				fit.0.coef 
				fit.0.sig
				fit.0.tpr
				fit.0.fpr 
				coef.array["ols",,b,i] = fit.0.coef
				sig.array["ols", b, i] = fit.0.sig
				tpr.array["ols", b, i] = fit.0.tpr
				fpr.array["ols", b, i] = fit.0.fpr
				res.array["ols", , b,  i] = c(mean(fit.0.z), sd(fit.0.z), 1)
	        }
					
			##### AdaReg with adapative gamma
			rob.result = AdaReg(X, y, gam.seq = seq(0, 3, by=.1), var.gp.id=rep(1, n), mu.fix.value=NULL, var.fix.value=NULL, tol=10^(-4), step=50) 
			fit.adarob.coef = rob.result[["beta.rob.fit"]]
			rob.result[["res.info"]]
			fit.adarob.sig = sqrt(rob.result[["var.sig.gp.fit"]])
			#hist(rob.result[["x.res"]])
			fit.adarob.out = (1:n)[abs(rob.result[["x.res"]]/fit.adarob.sig) > 2.5]
			fit.adarob.tpr = length(intersect(fit.adarob.out, o.ind))/length(o.ind)
			fit.adarob.fpr = length(intersect(fit.adarob.out, (1:n)[-o.ind]))/(n - length(o.ind))
			fit.adarob.coef 
			fit.adarob.sig
			fit.adarob.tpr
			fit.adarob.fpr 
			coef.array["w.ada.gam",,b,i] = fit.adarob.coef
			sig.array["w.ada.gam", b, i] = fit.adarob.sig
			tpr.array["w.ada.gam", b, i] = fit.adarob.tpr
			fpr.array["w.ada.gam", b, i] = fit.adarob.fpr
			if (!is.null(rob.result[["res.info"]])) {
			    res.array["w.ada.gam", , b,  i] = rob.result[["res.info"]][c("mu0.hat", "sd0.hat", "pi0.hat")]
			}
			
			
		    if (method.ind == 0) {
				#### AdaReg with gamma = 0.5
				rob.result.1 = AdaReg(X, y, gam.seq = 0.5, var.gp.id=rep(1, n), mu.fix.value=NULL, var.fix.value=NULL, tol=10^(-4), step=50) 
				fit.rob1.coef = rob.result.1[["beta.rob.fit"]]
				fit.rob1.sig = sqrt(rob.result.1[["var.sig.gp.fit"]])
				fit.rob1.out = (1:n)[abs(rob.result.1[["x.res"]]/fit.rob1.sig) > 2.5]
				fit.rob1.tpr = length(intersect(fit.rob1.out, o.ind))/length(o.ind)
				fit.rob1.fpr = length(intersect(fit.rob1.out, (1:n)[-o.ind]))/(n - length(o.ind))
				fit.rob1.coef 
				fit.rob1.sig
				fit.rob1.tpr
				fit.rob1.fpr 
				coef.array["w.gam=0.5",,b,i] = fit.rob1.coef
				sig.array["w.gam=0.5", b, i] = fit.rob1.sig
				tpr.array["w.gam=0.5", b, i] = fit.rob1.tpr
				fpr.array["w.gam=0.5", b, i] = fit.rob1.fpr
				if (!is.null(rob.result.1[["res.info"]])) {
				    res.array["w.gam=0.5", , b,  i] = rob.result.1[["res.info"]][c("mu0.hat", "sd0.hat", "pi0.hat")]
				}
				
				#### AdaReg with gamma = 1
				rob.result.2 = AdaReg(X, y, gam.seq = 1, var.gp.id=rep(1, n), mu.fix.value=NULL, var.fix.value=NULL, tol=10^(-4), step=50) 
				fit.rob2.coef = rob.result.2[["beta.rob.fit"]]
				fit.rob2.sig = sqrt(rob.result.2[["var.sig.gp.fit"]])
				fit.rob2.out = (1:n)[abs(rob.result.2[["x.res"]]/fit.rob2.sig) > 2.5]
				fit.rob2.tpr = length(intersect(fit.rob2.out, o.ind))/length(o.ind)
				fit.rob2.fpr = length(intersect(fit.rob2.out, (1:n)[-o.ind]))/(n - length(o.ind))
				fit.rob2.coef 
				fit.rob2.sig
				fit.rob2.tpr
				fit.rob2.fpr 
				coef.array["w.gam=1",,b,i] = fit.rob2.coef
				sig.array["w.gam=1", b, i] = fit.rob2.sig
				tpr.array["w.gam=1", b, i] = fit.rob2.tpr
				fpr.array["w.gam=1", b, i] = fit.rob2.fpr
				if (!is.null(rob.result.2[["res.info"]])) {
				    res.array["w.gam=1", , b,  i] = rob.result.2[["res.info"]][c("mu0.hat", "sd0.hat", "pi0.hat")]
				}	
				
			    #### AdaReg with gamma = 2
				rob.result.3 = AdaReg(X, y, gam.seq = 2, var.gp.id=rep(1, n), mu.fix.value=NULL, var.fix.value=NULL, tol=10^(-4), step=50) 
				fit.rob3.coef = rob.result.3[["beta.rob.fit"]]
				fit.rob3.sig = sqrt(rob.result.3[["var.sig.gp.fit"]])
				fit.rob3.out = (1:n)[abs(rob.result.3[["x.res"]]/fit.rob3.sig) > 2.5]
				fit.rob3.tpr = length(intersect(fit.rob3.out, o.ind))/length(o.ind)
				fit.rob3.fpr = length(intersect(fit.rob3.out, (1:n)[-o.ind]))/(n - length(o.ind))
				fit.rob3.coef 
				fit.rob3.sig
				fit.rob3.tpr
				fit.rob3.fpr 
				coef.array["w.gam=2",,b,i] = fit.rob3.coef
				sig.array["w.gam=2", b, i] = fit.rob3.sig
				tpr.array["w.gam=2", b, i] = fit.rob3.tpr
				fpr.array["w.gam=2", b, i] = fit.rob3.fpr
				if (!is.null(rob.result.3[["res.info"]])) {
				    res.array["w.gam=2", , b,  i] = rob.result.3[["res.info"]][c("mu0.hat", "sd0.hat", "pi0.hat")]
				}	
				
			    #### AdaReg with gamma = 3
				rob.result.4 = AdaReg(X, y, gam.seq = 3, var.gp.id=rep(1, n), mu.fix.value=NULL, var.fix.value=NULL, tol=10^(-4), step=50)
				fit.rob4.coef = rob.result.4[["beta.rob.fit"]]
				fit.rob4.sig = sqrt(rob.result.4[["var.sig.gp.fit"]])
				fit.rob4.out = (1:n)[abs(rob.result.4[["x.res"]]/fit.rob4.sig) > 2.5]
				fit.rob4.tpr = length(intersect(fit.rob4.out, o.ind))/length(o.ind)
				fit.rob4.fpr = length(intersect(fit.rob4.out, (1:n)[-o.ind]))/(n - length(o.ind))
				fit.rob4.coef 
				fit.rob4.sig
				fit.rob4.tpr
				fit.rob4.fpr 
				coef.array["w.gam=3",,b,i] = fit.rob4.coef
				sig.array["w.gam=3", b, i] = fit.rob4.sig
				tpr.array["w.gam=3", b, i] = fit.rob4.tpr
				fpr.array["w.gam=3", b, i] = fit.rob4.fpr
				if (!is.null(rob.result.4[["res.info"]])) {
				    res.array["w.gam=3", , b,  i] = rob.result.4[["res.info"]][c("mu0.hat", "sd0.hat", "pi0.hat")]
				}	
	
				#### LMS
				fit.lms = lqs(X, y, intercept=FALSE, method="lms")
				fit.lms.coef = fit.lms$coef
				fit.lms.sig = fit.lms$scale[1]
				fit.lms.z = (y - X %*% fit.lms.coef)/fit.lms.sig
				fit.lms.out = (1:n)[abs(fit.lms.z) > 2.5]
				fit.lms.tpr = length(intersect(fit.lms.out, o.ind))/length(o.ind)
				fit.lms.fpr = length(intersect(fit.lms.out, (1:n)[-o.ind]))/(n - length(o.ind))
				fit.lms.coef 
				fit.lms.sig
				fit.lms.tpr
				fit.lms.fpr 
				coef.array["lms",,b,i] = fit.lms.coef
				sig.array["lms", b, i] = fit.lms.sig
				tpr.array["lms", b, i] = fit.lms.tpr
				fpr.array["lms", b, i] = fit.lms.fpr
				
				#### LTS		
				fit.lts = lqs(X, y, intercept=FALSE, method="lts")
				fit.lts.coef = fit.lts$coef
				fit.lts.sig = fit.lts$scale[1]
				fit.lts.z = (y - X %*% fit.lts.coef)/fit.lts.sig
				fit.lts.out = (1:n)[abs(fit.lts.z) > 2.5]
				fit.lts.tpr = length(intersect(fit.lts.out, o.ind))/length(o.ind)
				fit.lts.fpr = length(intersect(fit.lts.out, (1:n)[-o.ind]))/(n - length(o.ind))
				fit.lts.coef 
				fit.lts.sig
				fit.lts.tpr
				fit.lts.fpr 
				coef.array["lts",,b,i] = fit.lts.coef
				sig.array["lts", b, i] = fit.lts.sig
				tpr.array["lts", b, i] = fit.lts.tpr
				fpr.array["lts", b, i] = fit.lts.fpr
					
				#### S 
				fit.s =lqs(X, y, intercept=FALSE, method="S")
				fit.s.coef = fit.s$coefficients
				fit.s.sig = fit.s$scale[1]
				fit.s.z = (y - X %*% fit.s.coef)/fit.s.sig
				fit.s.out = (1:n)[abs(fit.s.z) > 2.5]
				fit.s.tpr = length(intersect(fit.s.out, o.ind))/length(o.ind)
				fit.s.fpr = length(intersect(fit.s.out, (1:n)[-o.ind]))/(n - length(o.ind))
				fit.s.coef 
				fit.s.sig
				fit.s.tpr
				fit.s.fpr
				coef.array["S",,b,i] = fit.s.coef
				sig.array["S", b, i] = fit.s.sig
				tpr.array["S", b, i] = fit.s.tpr
				fpr.array["S", b, i] = fit.s.fpr
				
						
				#### Huber
				fit.huber = rlm(y ~ X - 1, psi=psi.huber)
				fit.huber.coef = fit.huber$coef
				fit.huber.sig = fit.huber$s
				fit.huber.z = (y - X %*% fit.huber.coef)/fit.huber.sig
				fit.huber.out = (1:n)[abs(fit.huber.z) > 2.5]
				fit.huber.tpr = length(intersect(fit.huber.out, o.ind))/length(o.ind)
				fit.huber.fpr = length(intersect(fit.huber.out, (1:n)[-o.ind]))/(n - length(o.ind))
				fit.huber.coef 
				fit.huber.sig
				fit.huber.tpr
				fit.huber.fpr 
				coef.array["w.huber",,b,i] = fit.huber.coef
				sig.array["w.huber", b, i] = fit.huber.sig
				tpr.array["w.huber", b, i] = fit.huber.tpr
				fpr.array["w.huber", b, i] = fit.huber.fpr
				
				
				#### Bisquare
				fit.bisq = rlm(y ~ X - 1, psi=psi.bisquare)
				fit.bisq.coef = fit.bisq$coef
				fit.bisq.sig = fit.bisq$s
				fit.bisq.z = (y - X %*% fit.bisq.coef)/fit.bisq.sig
				fit.bisq.out = (1:n)[abs(fit.bisq.z) > 2.5]
				fit.bisq.tpr = length(intersect(fit.bisq.out, o.ind))/length(o.ind)
				fit.bisq.fpr = length(intersect(fit.bisq.out, (1:n)[-o.ind]))/(n - length(o.ind))
				fit.bisq.coef 
				fit.bisq.sig
				fit.bisq.tpr
				fit.bisq.fpr 
				coef.array["w.bisq",,b,i] = fit.bisq.coef
				sig.array["w.bisq", b, i] = fit.bisq.sig
				tpr.array["w.bisq", b, i] = fit.bisq.tpr
				fpr.array["w.bisq", b, i] = fit.bisq.fpr
							
				#### Hampel
				fit.hampel = rlm(y ~ X - 1, psi=psi.hampel)
				fit.hampel.coef = fit.hampel$coef
				fit.hampel.sig = fit.hampel$s
				fit.hampel.z = (y - X %*% fit.hampel.coef)/fit.hampel.sig
				fit.hampel.out = (1:n)[abs(fit.hampel.z) > 2.5]
				fit.hampel.tpr = length(intersect(fit.hampel.out, o.ind))/length(o.ind)
				fit.hampel.fpr = length(intersect(fit.hampel.out, (1:n)[-o.ind]))/(n - length(o.ind))
				fit.hampel.coef 
				fit.hampel.sig
				fit.hampel.tpr
				fit.hampel.fpr 
				coef.array["w.hampel",,b,i] = fit.hampel.coef
				sig.array["w.hampel", b, i] = fit.hampel.sig
				tpr.array["w.hampel", b, i] = fit.hampel.tpr
				fpr.array["w.hampel", b, i] = fit.hampel.fpr
			}
		}
			
			
		for (j in 1:length(method.nm)) {
			rmse.coef.mx[j,i] = sqrt(mean(colSums( (coef.array[j,,,i] - beta.00)^2 ),na.rm=TRUE ))
			rmse.sig.mx[j,i] = sqrt(mean( (sig.array[j,,i] - sd.00)^2 ,na.rm=TRUE ))
			avg.tpr.mx[j,i] = mean( tpr.array[j,,i],na.rm=TRUE )
			avg.fpr.mx[j,i] = mean( fpr.array[j,,i],na.rm=TRUE )		
		}
				
	}
	
	#### diagnosis: to filter the cases when the residual mean > 1 for 'w.ada.gam'
	rmse.coef.mx.f = rmse.coef.mx
	rmse.sig.mx.f = rmse.sig.mx
	avg.tpr.mx.f = avg.tpr.mx
	avg.fpr.mx.f = avg.fpr.mx
	thres = 1
	res.test.mx = matrix(0, nrow(par.01.mx), 6)
	for (i in 1:nrow(par.01.mx)) {
	      res.test.mx[i, ] = rowSums(abs(res.array[, "res.mean.obv", ,  i]) >= thres, na.rm=TRUE)	
	}
	rownames(res.test.mx ) = paste("(pi1=", par.01.mx[,1], ", s=",par.01.mx[,2], ")", sep="")
	colnames(res.test.mx) = method.nm[1:6]
	#res.test.mx
	
	pbm.par = names(res.test.mx[res.test.mx[,"w.ada.gam"] > 0 ,"w.ada.gam"])
	if (length(pbm.par) > 0) {
		for (i in 1:length(pbm.par)) {
			pbm.par.1 = pbm.par[i] 
			ind = which(abs(res.array["w.ada.gam", "res.mean.obv", , pbm.par.1 ]) >= thres)
			res.array["w.ada.gam", "res.mean.obv", ind, pbm.par.1 ]
		    coef.array["w.ada.gam",, ind, pbm.par.1] 
		    sig.array["w.ada.gam",ind, pbm.par.1] 
		    tpr.array["w.ada.gam", ind,pbm.par.1]
		    fpr.array["w.ada.gam", ind,pbm.par.1]    
		    
			coef.array.f = coef.array 
			coef.array.f["w.ada.gam",, ind, pbm.par.1]  = NA 
		    sig.array.f = sig.array 
			sig.array.f["w.ada.gam",ind, pbm.par.1]  = NA 
			tpr.array.f = tpr.array 
			tpr.array.f["w.ada.gam",ind, pbm.par.1]  = NA 
			fpr.array.f = fpr.array 
			fpr.array.f["w.ada.gam",ind, pbm.par.1]  = NA 
			
			for (j in 1:length(method.nm)) {
				rmse.coef.mx.f[j,pbm.par.1] = sqrt(mean(colSums( (coef.array.f[j,,,pbm.par.1] - beta.00)^2 ),na.rm=TRUE ))
				rmse.sig.mx.f[j,pbm.par.1] = sqrt(mean( (sig.array.f[j,,pbm.par.1] - sd.00)^2 ,na.rm=TRUE ))
				avg.tpr.mx.f[j,pbm.par.1] = mean( tpr.array.f[j,,pbm.par.1],na.rm=TRUE )
				avg.fpr.mx.f[j,pbm.par.1] = mean( fpr.array.f[j,,pbm.par.1],na.rm=TRUE )		
		    }
			
		}
	}

	
    return(list(par.01.mx=par.01.mx, method.nm=method.nm, coef.array=coef.array, sig.array=sig.array, tpr.array=tpr.array, fpr.array=fpr.array, rmse.coef.mx=rmse.coef.mx, rmse.sig.mx=rmse.sig.mx, avg.tpr.mx=avg.tpr.mx, avg.fpr.mx=avg.fpr.mx, res.array= res.array, res.test.mx=res.test.mx, pbm.par=pbm.par, rmse.coef.mx.f = rmse.coef.mx.f, rmse.sig.mx.f = rmse.sig.mx.f, avg.tpr.mx.f = avg.tpr.mx.f, avg.fpr.mx.f = avg.fpr.mx.f, BB=BB, n=n, p=p))
}


# ------------------------------------------------------------------------------
# to plot RMSE on coefficients in regression model from various methods
plot.rmse.coef.fn = function(result, file.rmse.coef.nm, noise.model.nm) {
	my.col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50) # "Spectral"
    
    method.nm = result[['method.nm']]
    BB = result[['BB']]
    n = result[['n']]
    p = result[['p']]
	rmse.coef.mx.f = result[['rmse.coef.mx.f']]	
	rmse.coef.mx.1 = rmse.coef.mx.f
	rmse.coef.mx.1[rmse.coef.mx.1 > 2 & !is.na(rmse.coef.mx.1)] = 2
	rmse.coef.mx.1[is.na(rmse.coef.mx.1)] = 2
	
	pdf(file=file.rmse.coef.nm, height= 5.5, width=12)
	par(mfrow=c(2,6), mai=c(0.6, 0.6, 0.5, 0.1), oma=c(0,0,3,0))
	for(j in 1:length(method.nm)) {
		nm = method.nm[j]
		mx = rbind(rep(rmse.coef.mx.1[nm, 1], 5), cbind(rep(rmse.coef.mx.1[nm, 1], 10), matrix(rmse.coef.mx.1[nm,-1], 10, 4)))
		colnames(mx) = paste("pi1=", seq(0, 0.4, by=0.1), sep="")
		rownames(mx) = paste("s=", 0:10, sep="")

		
		image(mx, col =my.col, xaxt="n", yaxt="n", main=paste(nm), zlim=range(rmse.coef.mx.1), xlab="s", ylab="pi.1", cex.main=1.5, cex.lab=1.5) # terrain.colors(100) # xlab='s'
		axis(1, at=seq(0, 1, len=11), 0:10)
		axis(2, at=seq(0, 1, len=5), seq(0, 0.4, by=0.1))
		contour(round(mx,2), levels=seq(0, 2, by=0.1),  add = TRUE, col = "black", lwd=1, labcex=0.8) # levels = seq(min(round(sqrt(mx),2)), max(round(sqrt(mx),2)), by=.1)
	}
	mtext(paste("comparisons in RMSE for beta.0 under n = ",n, ", p = ", p, ", B = ", BB, '\n noise ~ ', noise.model.nm ,sep=""), outer=TRUE, cex=1.2)
	dev.off()
	
}

# ------------------------------------------------------------------------------
# to plot RMSE on sigma.0 in regression model from various methods
plot.rmse.sigma.fn = function(result, file.rmse.sigma.nm, noise.model.nm) {
	my.col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50) # "Spectral"
    
    method.nm = result[['method.nm']]
    BB = result[['BB']]
    n = result[['n']]
    p = result[['p']]
	rmse.sig.mx.f = result[['rmse.sig.mx.f']]	
	rmse.sig.mx.1 = rmse.sig.mx.f
	rmse.sig.mx.1[rmse.sig.mx.1 > 2 & !is.na(rmse.sig.mx.1)] = 2
	rmse.sig.mx.1[is.na(rmse.sig.mx.1)] = 2

	pdf(file=file.rmse.sigma.nm, height= 5.5, width=12)
	par(mfrow=c(2,6), mai=c(0.6, 0.6, 0.5, 0.1), oma=c(0,0,3,0))
	for(j in 1:length(method.nm)) {
		nm = method.nm[j]
		mx = rbind(rep(rmse.sig.mx.1[nm, 1], 5), cbind(rep(rmse.sig.mx.1[nm, 1], 10), matrix(rmse.sig.mx.1[nm,-1], 10, 4)))
		colnames(mx) = paste("pi1=", seq(0, 0.4, by=0.1), sep="")
		rownames(mx) = paste("s=", 0:10, sep="")
		
		image(mx, col =my.col, xaxt="n", yaxt="n", main=paste(nm), zlim=range(rmse.sig.mx.1), xlab="s", ylab="pi.1", cex.main=1.5, cex.lab=1.5) # terrain.colors(100)
		axis(1, at=seq(0, 1, len=11), 0:10)
		axis(2, at=seq(0, 1, len=5), seq(0, 0.4, by=0.1))
		contour(round(mx,2), levels=seq(0, 2, by=0.1),  add = TRUE, col = "black", lwd=1, labcex=0.8) # levels = seq(min(round(sqrt(mx),2)), max(round(sqrt(mx),2)), by=.1)
	}
	mtext(paste("comparisons in RMSE for sigma.0 under n = ",n, ", p = ", p, ", B = ", BB, '\n noise ~ ', noise.model.nm ,sep=""), outer=TRUE, cex=1.2)
	dev.off()
	
}


# ------------------------------------------------------------------------------
# to plot FPR.TPR in regression model from various methods
plot.fpr.tpr.fn = function(result, file.fpr.tpr.nm, noise.model.nm) {
	my.col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50) # "Spectral"
    
    method.nm = result[['method.nm']]
    BB = result[['BB']]
    n = result[['n']]
    p = result[['p']]
    avg.fpr.mx = result[['avg.fpr.mx']]
    avg.tpr.mx = result[['avg.tpr.mx']]
            
    pdf(file=file.fpr.tpr.nm, height= 5.5, width=12)
    par(mfrow=c(2,6), mai=c(0.6, 0.6, 0.5, 0.1), oma=c(0,0,3,0))
	col.mu1 = rainbow(10)
	col.ls = rep(col.mu1, 4)
	pch.ls = rep(1:4, each=10)
	for(j in 1:length(method.nm)) {
		nm=method.nm[j]
		plot(avg.fpr.mx[nm, -1], avg.tpr.mx[nm, -1], xlim=c(0,max(avg.fpr.mx[-1], na.rm=TRUE)+0.15), ylim=c(0,1), main=nm, col=col.ls, pch=pch.ls, xlab="FPR", ylab="TPR", cex.lab=1.5, cex.main=1.5, cex=1.5)
	
		if (j %in% c(1, 2,3,4)) {
			legend("topright", colnames(avg.fpr.mx)[-1][((j-1)*10+1):(j*10)], col=col.ls[((j-1)*10+1):(j*10)], pch=pch.ls[((j-1)*10+1):(j*10)], cex=0.8)
		}
	}
	
	mtext(paste("comparisons in (FPR, TPR) under n = ",n, ", p = ", p, ", B = ", BB, '\n noise ~ ', noise.model.nm, sep=""), outer=TRUE, cex=1.2)
	dev.off()
	
}

