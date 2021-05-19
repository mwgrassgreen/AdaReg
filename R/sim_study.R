#' @title sim_fn.R
#' @description Simulation studies
#' @author Meng Wang
#' \email{mengw1@stanford.edu}


source("AdaReg.R")
source('sim_fn.R')

# ---------------------#
#    norm-norm model   #
# ---------------------#
n = 200
step = 50
BB = 50
pop.distr = 'norm'
out.distr = 'norm'
mu0 = 0
sd0 = 1
gam.seq = seq(0, 3, by=0.5)
result = theta.trend.est.fn (n, step, BB, gam.seq, pop.distr, out.distr, norm.mu0=mu0, norm.sd0=sd0)
fdr.null.int = efdr.null.int.fn(par.01.mx.comb=result[['par.01.mx.comb']], pop.distr, out.distr, mu.0=mu0, sd.0=sd0)

# please specify 'my_dir' to run the code
mix.model.nm = paste('(1 - pi_1) N(', mu0,', ', sd0, ') + pi_1 N(mu_1, ', sd0, ')', sep='')
file.trend.nm = paste(my_dir, '/trend_comparison_', pop.distr, "_mixed_", out.distr, "_n", n, sep='')
file.theta.nm = paste(my_dir, '/theta_comparison_', pop.distr, "_mixed_", out.distr, "_n", n, sep='')
plot.crt.trend.fn(gam.seq, result, fdr.null.int, mix.model.nm, file.nm=file.trend.nm) 
plot.theta.fn(gam.seq, result, mix.model.nm, file.nm=file.theta.nm) 

n = 2000
step = 50
BB = 50
pop.distr = 'norm'
out.distr = 'norm'
mu0 = 0
sd0 = 1
gam.seq = seq(0, 3, by=0.5)
result = theta.trend.est.fn (n, step, BB, gam.seq, pop.distr, out.distr, norm.mu0=mu0, norm.sd0=sd0)
fdr.null.int = efdr.null.int.fn(par.01.mx.comb=result[['par.01.mx.comb']], pop.distr, out.distr, mu.0=mu0, sd.0=sd0)

mix.model.nm = paste('(1 - pi_1) N(', mu0,', ', sd0, ') + pi_1 N(mu_1, ', sd0, ')', sep='')
file.trend.nm = paste(my_dir, '/trend_comparison_', pop.distr, "_mixed_", out.distr, "_n", n, sep='')
file.theta.nm = paste(my_dir, '/theta_comparison_', pop.distr, "_mixed_", out.distr, "_n", n, sep='')
plot.crt.trend.fn(gam.seq, result, fdr.null.int, mix.model.nm, file.nm=file.trend.nm) 
plot.theta.fn(gam.seq, result, mix.model.nm, file.nm=file.theta.nm) 


# ---------------------#
#     norm-t model     #
# ---------------------#
n = 2000
step = 50
BB = 50
pop.distr = 'norm'
out.distr = 't'
mu0 = 0
sd0 = 1
df1 = 3 #5
gam.seq = seq(0, 3, by=0.5)
result = theta.trend.est.fn (n, step, BB, gam.seq, pop.distr, out.distr, norm.mu0=mu0, norm.sd0=sd0, t.df1=df1)
fdr.null.int = efdr.null.int.fn(par.01.mx.comb=result[['par.01.mx.comb']], pop.distr, out.distr, mu.0=mu0, sd.0=sd0, t.df.1=df1)

mix.model.nm = paste('(1 - pi_1) N(', mu0,', ', sd0, ') + pi_1 (mu_1 + t(', df1, '))', sep='')
file.trend.nm = paste(my_dir, '/trend_comparison_', pop.distr, "_mixed_", out.distr, "_n", n, sep='')
file.theta.nm = paste(my_dir, '/theta_comparison_', pop.distr, "_mixed_", out.distr, "_n", n, sep='')
plot.crt.trend.fn(gam.seq, result, fdr.null.int, mix.model.nm, file.nm=file.trend.nm) 
plot.theta.fn(gam.seq, result, mix.model.nm, file.nm=file.theta.nm) 


# ---------------------#
#     t-norm model     #
# ---------------------#
n = 2000
step = 50
BB = 50
pop.distr = 't'
out.distr = 'norm'
shift0 = 0
df0 = 3
gam.seq = seq(0, 3, by=0.5)
result = theta.trend.est.fn (n, step, BB, gam.seq, pop.distr, out.distr, t.shift0=shift0, t.df0=df0)
fdr.null.int = efdr.null.int.fn(par.01.mx.comb=result[['par.01.mx.comb']], pop.distr, out.distr, t.shift.0=shift0, t.df.0=df0)

mix.model.nm = paste('(1 - pi_1) t(', df0, ') + pi_1 N(mu_1, sigma(t))', sep='')
file.trend.nm = paste(my_dir, '/trend_comparison_', pop.distr, "_mixed_", out.distr, "_n", n, sep='')
file.theta.nm = paste(my_dir, '/theta_comparison_', pop.distr, "_mixed_", out.distr, "_n", n, sep='')
plot.crt.trend.fn(gam.seq, result, fdr.null.int, mix.model.nm, file.nm=file.trend.nm) 
plot.theta.fn(gam.seq, result, mix.model.nm, file.nm=file.theta.nm) 


# ---------------------#
#   gamma-norm model   #
# ---------------------#
n = 2000
step = 50
BB = 50
pop.distr = 'gamma'
out.distr = 'norm'
shape0 = 5 # 7 
rate0 = 1
gam.seq = seq(0, 3, by=0.5)
result = theta.trend.est.fn (n, step, BB, gam.seq, pop.distr, out.distr, gamma.shape0=shape0, gamma.rate0=rate0)
fdr.null.int = efdr.null.int.fn(par.01.mx.comb=result[['par.01.mx.comb']], pop.distr, out.distr, shape.0=shape0, rate.0=rate0)

mix.model.nm = paste('(1 - pi_1) gamma(', shape0, ',', rate0,') + pi_1 N(mu_1, sigma(gamma))', sep='')
file.trend.nm = paste(my_dir, '/trend_comparison_', pop.distr, "_mixed_", out.distr, "_n", n, sep='')
file.theta.nm = paste(my_dir, '/theta_comparison_', pop.distr, "_mixed_", out.distr, "_n", n, sep='')
plot.crt.trend.fn(gam.seq, result, fdr.null.int, mix.model.nm, file.nm=file.trend.nm) 
plot.theta.fn(gam.seq, result, mix.model.nm, file.nm=file.theta.nm) 


# ---------------------#
#  lognorm-norm model  #
# ---------------------#
n = 2000
step = 50
BB = 50
pop.distr = 'lognorm'
out.distr = 'norm'
logmean0 = 0
logsd0 = 0.25 # skewness increasing with logsd^2
gam.seq = seq(0, 3, by=0.5)
result = theta.trend.est.fn (n, step, BB, gam.seq, pop.distr, out.distr,lognorm.meanlog0=logmean0, lognorm.sdlog0=logsd0)
fdr.null.int = efdr.null.int.fn(par.01.mx.comb=result[['par.01.mx.comb']], pop.distr, out.distr, meanlog.0=logmean0, sdlog.0=logsd0)

mix.model.nm = paste('(1 - pi_1) lognorm(', logmean0, ',', logsd0,') + pi_1 N(mu_1, sigma(lognorm))', sep='')
file.trend.nm = paste(my_dir, '/trend_comparison_', pop.distr, "_mixed_", out.distr, "_n", n, sep='')
file.theta.nm = paste(my_dir, '/theta_comparison_', pop.distr, "_mixed_", out.distr, "_n", n, sep='')
plot.crt.trend.fn(gam.seq, result, fdr.null.int, mix.model.nm, file.nm=file.trend.nm) 
plot.theta.fn(gam.seq, result, mix.model.nm, file.nm=file.theta.nm) 



#------------------- REGRESSION MODEL -------------------# 
library(WRS2)
library(MASS)
library(RColorBrewer)


# --------------------------#
#      regression model     #
#    with norm-norm noise   #
# --------------------------#

n = 200 
p = 2
pop.distr='norm'
out.distr='norm'
sd.0 = 1
BB = 100
noise.model.nm = paste('(1 - pi_1) N(0,', sd.0,') + pi_1 N(mu_1,',  sd.0, ')', sep='')
result = est.comp.in.reg.model.fn(n, p, pop.distr, out.distr, sd.0=sd.0, BB=BB)

####
n = 2000
p = 2
pop.distr='norm'
out.distr='norm'
sd.0 = 1
BB = 100
noise.model.nm = paste('(1 - pi_1) N(0,', sd.0,') + pi_1 N(mu_1,',  sd.0, ')', sep='')
result = est.comp.in.reg.model.fn(n, p, pop.distr, out.distr, sd.0=sd.0, BB=BB)


####
n = 2000
p = 20
pop.distr='norm'
out.distr='norm'
sd.0 = 1
BB = 100
noise.model.nm = paste('(1 - pi_1) N(0,', sd.0,') + pi_1 N(mu_1,',  sd.0, ')', sep='')
result = est.comp.in.reg.model.fn(n, p, pop.distr, out.distr, sd.0=sd.0, BB=BB)


#### to plot results above

file.rmse.coef.nm = paste(my_dir, '/rmse_coef_comparison_n', n, "_p", p, '_under_noise_', pop.distr, "_mixed_", out.distr, ".pdf", sep='')
file.rmse.sigma.nm = paste(my_dir, '/rmse_sigma_comparison_n', n, "_p", p, '_under_noise_', pop.distr, "_mixed_", out.distr, ".pdf", sep='')
file.fpr.tpr.nm = paste(my_dir, '/fpr_tpr_comparison_n', n, "_p", p, '_under_noise_', pop.distr, "_mixed_", out.distr, ".pdf", sep='')

plot.rmse.coef.fn(result, file.rmse.coef.nm, noise.model.nm)
plot.rmse.sigma.fn(result, file.rmse.sigma.nm, noise.model.nm)
plot.fpr.tpr.fn(result, file.fpr.tpr.nm, noise.model.nm)


# --------------------------#
#      regression model     #
#    with norm-t noise      #
# --------------------------#

n = 2000 
p = 2
pop.distr='norm'
out.distr='t'
sd.0 = 1
t.df.1 = 3
BB = 100
noise.model.nm = paste('(1 - pi_1) N(0,', sd.0,') + pi_1 (mu_1 + t(',  t.df.1, '))', sep='')
result = est.comp.in.reg.model.fn(n, p, pop.distr, out.distr, sd.0=sd.0, t.df.1=t.df.1, BB=BB)


# --------------------------#
#      regression model     #
#    with t-norm noise      #
# --------------------------#

n = 2000
p = 2
pop.distr='t'
out.distr='norm'
t.df.0 = 5 # 3
BB = 100
noise.model.nm = paste('(1 - pi_1) t(', t.df.0,') + pi_1 N(mu_1, sigma(t))', sep='')
result = est.comp.in.reg.model.fn(n, p, pop.distr, out.distr, t.df.0=t.df.0, BB=BB)


# --------------------------#
#      regression model     #
#   with gamma-norm noise   #
# --------------------------#

n = 2000
p = 2
pop.distr = 'gamma'
out.distr = 'norm'
shape0 = 5 # 7
rate0 = 1
BB = 100
noise.model.nm = paste('(1 - pi_1) gamma(', shape0, ',', rate0,') + pi_1 N(mu_1, sigma(gamma)) - mu(gamma)', sep='')
result = est.comp.in.reg.model.fn(n, p, pop.distr, out.distr, shape.0=shape0, rate.0 =rate0, BB=BB, method.ind=1) # only AdaReg


# --------------------------#
#      regression model     #
#  with lognorm-norm noise  #
# --------------------------#

n = 2000
p = 2
pop.distr = 'lognorm'
out.distr = 'norm'
logmean0 = 0
logsd0 = 0.25 # 0.3
BB = 100
noise.model.nm = paste('(1 - pi_1) lognorm(', logmean0, ',', logsd0,') + pi_1 N(mu_1, sigma(lognorm)) - mu(lognorm)', sep='')
result = est.comp.in.reg.model.fn(n, p, pop.distr, out.distr, meanlog.0=logmean0, sdlog.0=logsd0, BB=BB, method.ind=1) # only AdaReg

