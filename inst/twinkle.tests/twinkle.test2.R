#################################################################################
##
##   R package twinkle by Alexios Ghalanos Copyright (C) 2014.
##   This file is part of the R package twinkle.
##
##   The R package twinkle is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package twinkle is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################


twinkle.test2a = function(cluster = NULL){
	tic = Sys.time()
	require(rugarch)
	# 1-regime
	set.seed(12)
	sim = fitted(arfimapath(arfimaspec(mean.model=list(armaOrder=c(2,0)), fixed.pars=list(mu=0.04, ar1 = 0.8, ar2=-0.1, sigma=0.1)), n.sim=2000))
	dat = xts(sim, as.Date(1:2000))
	dat[500] = -0.7
	dat[501] = 0.3
	dat[502] = 0.3
	dat[503] = -0.7	
	fp = rep(1, 2000)
	fp[500:503] = 0
	
	fp = xts(fp, index(dat))
	spec = starspec(
			mean.model = list(regimes = 1, include.intercept = c(1), arOrder = c(2)),
			variance.model = list(dynamic = FALSE, model = "sGARCH", 
					garchOrder = c(1, 1), submodel = NULL, external.regressors = NULL, 
					variance.targeting = FALSE),
			distribution.model = "norm", fixed.prob=fp, intmodel=1)
	fit0 = starfit(spec, data = dat, out.sample = 0, solver = "optim", solver.control = solver.control)
	
	fp = xts(rep(1, 2000), index(dat))
	spec = starspec(
			mean.model = list(regimes = 1, include.intercept = c(1), arOrder = c(2)),
			variance.model = list(dynamic = FALSE, model = "sGARCH", 
					garchOrder = c(1, 1), submodel = NULL, external.regressors = NULL, 
					variance.targeting = FALSE),
			distribution.model = "norm", fixed.prob=fp, intmodel=1)
	solver.control=list(maxit=150000, alpha=1, beta=0.4, gamma=2, reltol=1e-12, trace=1,method="BFGS")
	fit1 = starfit(spec, data = dat, out.sample = 0, solver = "optim", solver.control = solver.control)
	fit2 = arfimafit(arfimaspec(mean.model=list(armaOrder=c(2,0))), data = dat)
	
	residuals(fit0)[500]
	residuals(fit1)[500]
	residuals(fit2)[500]
	
	sum(residuals(fit0)^2)
	sum(residuals(fit1)^2)
	sum(residuals(fit2)^2)
	
	rbind(coef(fit0), coef(fit1), coef(fit2))
	c(likelihood(fit0), likelihood(fit1), likelihood(fit2))
}