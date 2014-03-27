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
	data(forex)
	# to replicate the results of van Dijk and Franses use:
	# next observation carried backward (i.e. fromLast=TRUE)
	# for missing values
	fx = na.locf(forex, fromLast = TRUE)
	fxw = fx[which(weekdays(index(forex))=="Wednesday"),4]
	dx = ROC(fxw, na.pad=FALSE)*100
	# threshold variable a running mean of the absolute returns
	dxv = runMean(abs(dx), n=4)
	dxv[1:3] = c(as.numeric(abs(dx[1,1])), as.numeric(mean(abs(dx[1:2,1]))), 
			as.numeric(mean(abs(dx[1:3,1]))))
	# period considered
	idx = 1:521
	spec1 = starspec(
			mean.model = list(regimes = 2, include.intercept = c(1,1), arOrder = c(0,2), 
					type = c("y","x")[2], x = dxv[idx], dstar = F, lags = 1, 
					external.regressors = NULL, fun = NULL),
			variance.model = list(dynamic = FALSE, model = "sGARCH", 
					garchOrder = c(1, 1), submodel = NULL, external.regressors = NULL, 
					variance.targeting = FALSE), 
			distribution.model = "std")
	solver.control=list(maxit=150000, alpha=1, beta=0.4, gamma=1.2, reltol=1e-12, trace=1,method="BFGS",
			n.restarts = 4)
	set.seed(5)
	fit1 = starfit(spec1, data = dx[idx], out.sample = 0, solver = "msoptim", solver.control = solver.control)
	
	
	spec2 = starspec(
			mean.model = list(regimes = 2, include.intercept = c(1,1), arOrder = c(0,2), 
					type = c("y","x")[2], x = dxv[idx], dstar = F, lags = 1, 
					external.regressors = NULL, fun = NULL),
			variance.model = list(dynamic = T, model = "sGARCH", 
					garchOrder = c(1, 1), submodel = NULL, external.regressors = NULL, 
					variance.targeting = T), 
			distribution.model = "std")
	solver.control=list(maxit=150000, alpha=1, beta=0.4, gamma=1.2, reltol=1e-12, trace=1,method="BFGS",
			n.restarts = 4)
	set.seed(5)
	fit2 = starfit(spec2, data = dx[idx], out.sample = 0, solver = "msoptim", solver.control = solver.control)
	
	spec3 = ugarchspec(mean.model=list(armaOrder=c(2,0)), variance.model=list(variance.targeting=TRUE), distribution = "std")
	fit3 = ugarchfit(spec3, dx[idx], solver="nlminb")
	
	plot(dx[idx])
	lines(quantile(fit1, 0.05), col=2)
	lines(quantile(fit3, 0.05), col=4)
	lines(quantile(fit2, 0.05), col=3)
	
	
}