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

# Simulation Tests:
# -- path versus sim
# -- unconditional
# -- methods and plots
# -- GIRF

# static variance
twinkle.test2a = function(cluster = NULL){
	tic = Sys.time()
	require(twinkle)
	data(forex)	
	
	# 1. Equivalence of sim and path methods for self-exciting model with yfun
	
	# load("./data/forex.rda")
	fx = na.locf(forex, fromLast = TRUE)
	fxw = fx[which(weekdays(index(forex))=="Wednesday"),4]
	dx = ROC(fxw, na.pad=FALSE)*100
	fun = function(x){
		x = as.numeric(x)
		N = length(x)
		if(N<4){
			y = abs(x)	
		} else{
			y = runMean(abs(x), n=4)
			y[1:3] = c(abs(x[1]), mean(abs(x[1:2])), mean(abs(x[1:3])))
		}
		return(y)
	}
	
	arorder=list(c(0,0),c(0,1), c(1,0), c(1,1), c(2,1))
	int = list(c(0,0), c(1,0))
	test = matrix(NA, 2, 5)
	for(i in 1:2){
		for(j in 1:5){
			spec = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[1], ylags = 1, yfun = fun))
			solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
			fit = starfit(spec, data = dx[1:521], solver = "strategy", solver.control = solver.control, n=5)
			sim  = starsim(fit, n.sim=5000)
			specx = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]],
							statevar = c("y","s")[1], ylags = 1, yfun = fun))
			setfixed(specx)<-as.list(coef(fit))
			path = starpath(specx, n.sim=5000, prereturns = tail(dx[1:521], 6), 
					custom.dist=list(name="sample", distfit=sim@simulation$residSim/coef(fit)["sigma"]))			
			test[i,j] = as.integer(all.equal(fitted(sim), fitted(path)))
			#print(test[i,j])
		}
	}
	zz <- file("test2a-1.txt", open="wt")
	sink(zz)
	print(test)
	sink(type="message")
	sink()
	close(zz)
	
	
	# 2. Equivalence of sim and path methods for self-exciting model without yfun
	arorder=list(c(0,0),c(0,1), c(1,0), c(1,1), c(2,1))
	int = list(c(0,0), c(1,0))
	test = matrix(NA, 2, 5)
	for(i in 1:2){
		for(j in 1:5){
			spec = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[1], ylags = 1))
			solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
			fit = starfit(spec, data = dx[1:521], solver = "strategy", solver.control = solver.control, n=5)
			sim  = starsim(fit, n.sim=5000)
			specx = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]],
							statevar = c("y","s")[1], ylags = 1))
			setfixed(specx)<-as.list(coef(fit))
			path = starpath(specx, n.sim=5000, prereturns = tail(dx[1:521], 6), 
					custom.dist=list(name="sample", distfit=sim@simulation$residSim/coef(fit)["sigma"]))
			test[i,j] = as.integer(all.equal(fitted(sim), fitted(path)))
			#print(test[i,j])
		}
	}
	zz <- file("test2a-2.txt", open="wt")
	sink(zz)
	print(test)
	sink(type="message")
	sink()
	close(zz)
	
	
	# 3. Equivalence of sim and path methods for s-var model
	dxv = runMean(abs(dx), n=4)
	dxv[1:3] = c(as.numeric(abs(dx[1,1])), as.numeric(mean(abs(dx[1:2,1]))), 
			as.numeric(mean(abs(dx[1:3,1]))))
	dxv = lag(dxv, 1)
	dxv[1] = 0
	arorder=list(c(0,0),c(0,1), c(1,0), c(1,1), c(2,1))
	int = list(c(0,0), c(1,0))
	test = matrix(NA, 2, 5)
	
	for(i in 1:2){
		for(j in 1:5){
			spec = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[2], s = dxv[1:521]))
			solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
			fit = starfit(spec, data = dx[1:521], solver = "strategy", solver.control = solver.control, n=5)
			ssim = sample(as.numeric(dxv), 5000, replace=TRUE)
			sim  = starsim(fit, n.sim=5000, ssim = list(matrix(ssim, ncol=1)))
			specx = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[2], s = dxv[500:521]))
			setfixed(specx)<-as.list(coef(fit))
			path = starpath(specx, n.sim=5000, prereturns = tail(dx[1:521], 6), 
					custom.dist=list(name="sample", distfit=sim@simulation$residSim/coef(fit)["sigma"]),
					ssim= list(matrix(ssim, ncol=1)))
			test[i,j] = as.integer(all.equal(fitted(sim), fitted(path)))
		}
	}
	zz <- file("test2a-3.txt", open="wt")
	sink(zz)
	print(test)
	sink(type="message")
	sink()
	close(zz)
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

# mixture variance
twinkle.test2b = function(cluster = NULL){
	tic = Sys.time()
	require(twinkle)
	data(forex)	
	
	# 1. Equivalence of sim and path methods for self-exciting model with yfun
	
	# load("./data/forex.rda")
	fx = na.locf(forex, fromLast = TRUE)
	fxw = fx[which(weekdays(index(forex))=="Wednesday"),4]
	dx = ROC(fxw, na.pad=FALSE)*100
	fun = function(x){
		x = as.numeric(x)
		N = length(x)
		if(N<4){
			y = abs(x)	
		} else{
			y = runMean(abs(x), n=4)
			y[1:3] = c(abs(x[1]), mean(abs(x[1:2])), mean(abs(x[1:3])))
		}
		return(y)
	}
	
	arorder=list(c(0,0),c(0,1), c(1,0), c(1,1), c(2,1))
	int = list(c(0,0), c(1,0))
	test = matrix(NA, 2, 5)
	for(i in 1:2){
		for(j in 1:5){
			spec = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[1], ylags = 1, yfun = fun),
					variance.model=list(dynamic=TRUE, model="mixture"))
			solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
			fit = starfit(spec, data = dx[1:521], solver = "strategy", solver.control = solver.control, n=5)
			sim  = starsim(fit, n.sim=5000)
			specx = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]],
							statevar = c("y","s")[1], ylags = 1, yfun = fun),
					variance.model=list(dynamic=TRUE, model="mixture"))
			setfixed(specx)<-as.list(coef(fit))
			path = starpath(specx, n.sim=5000, prereturns = tail(dx[1:521], 12), preresiduals = tail(residuals(fit),2),
					custom.dist=list(name="sample", distfit=sim@simulation$residSim/sigma(sim)))
			test[i,j] = as.integer(all.equal(fitted(sim), fitted(path)))
			print(test[i,j])
		}
	}
	zz <- file("test2b-1.txt", open="wt")
	sink(zz)
	print(test)
	sink(type="message")
	sink()
	close(zz)
	
	
	# 2. Equivalence of sim and path methods for self-exciting model without yfun
	arorder=list(c(0,0),c(0,1), c(1,0), c(1,1), c(2,1))
	int = list(c(0,0), c(1,0))
	test = matrix(NA, 2, 5)
	for(i in 1:2){
		for(j in 1:5){
			spec = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[1], ylags = 1),
					variance.model=list(dynamic=TRUE, model="mixture"))
			solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
			fit = starfit(spec, data = dx[1:521], solver = "strategy", solver.control = solver.control, n=5)
			sim  = starsim(fit, n.sim=5000)
			specx = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]],
							statevar = c("y","s")[1], ylags = 1),
					variance.model=list(dynamic=TRUE, model="mixture"))
			setfixed(specx)<-as.list(coef(fit))
			path = starpath(specx, n.sim=5000, prereturns = tail(dx[1:521], 6), preresiduals = tail(residuals(fit),2),
					custom.dist=list(name="sample", distfit=sim@simulation$residSim/sigma(sim)))
			test[i,j] = as.integer(all.equal(fitted(sim), fitted(path)))
			print(test[i,j])
		}
	}
	zz <- file("test2b-2.txt", open="wt")
	sink(zz)
	print(test)
	sink(type="message")
	sink()
	close(zz)
	
	
	# 3. Equivalence of sim and path methods for s-var model
	dxv = runMean(abs(dx), n=4)
	dxv[1:3] = c(as.numeric(abs(dx[1,1])), as.numeric(mean(abs(dx[1:2,1]))), 
			as.numeric(mean(abs(dx[1:3,1]))))
	dxv = lag(dxv, 1)
	dxv[1] = 0
	arorder=list(c(0,0),c(0,1), c(1,0), c(1,1), c(2,1))
	int = list(c(0,0), c(1,0))
	test = matrix(NA, 2, 5)
	
	for(i in 1:2){
		for(j in 1:5){
			spec = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[2], s = dxv[1:521]),
					variance.model=list(dynamic=TRUE, model="mixture"))
			solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
			fit = starfit(spec, data = dx[1:521], solver = "strategy", solver.control = solver.control, n=5)
			ssim = sample(as.numeric(dxv), 5000, replace=TRUE)
			sim  = starsim(fit, n.sim=5000, ssim = list(matrix(ssim, ncol=1)))
			specx = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[2], s = dxv[500:521]),
					variance.model=list(dynamic=TRUE, model="mixture"))
			setfixed(specx)<-as.list(coef(fit))
			path = starpath(specx, n.sim=5000, prereturns = tail(dx[1:521], 6), 
					custom.dist=list(name="sample", distfit=sim@simulation$residSim/sigma(sim)),
					ssim=list(matrix(ssim, ncol=1)))
			test[i,j] = as.integer(all.equal(fitted(sim), fitted(path)))
		}
	}
	zz <- file("test2b-3.txt", open="wt")
	sink(zz)
	print(test)
	sink(type="message")
	sink()
	close(zz)
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

# GARCH variance
twinkle.test2c = function(cluster = NULL){
	tic = Sys.time()
	require(twinkle)
	data(forex)	
	
	# 1. Equivalence of sim and path methods for self-exciting model with yfun
	
	# load("./data/forex.rda")
	fx = na.locf(forex, fromLast = TRUE)
	fxw = fx[which(weekdays(index(forex))=="Wednesday"),4]
	dx = ROC(fxw, na.pad=FALSE)*100
	fun = function(x){
		x = as.numeric(x)
		N = length(x)
		if(N<4){
			y = abs(x)	
		} else{
			y = runMean(abs(x), n=4)
			y[1:3] = c(abs(x[1]), mean(abs(x[1:2])), mean(abs(x[1:3])))
		}
		return(y)
	}
	
	arorder=list(c(0,0),c(0,1), c(1,0), c(1,1), c(2,1))
	int = list(c(0,0), c(1,0))
	test = matrix(NA, 2, 5)
	for(i in 1:2){
		for(j in 1:5){
			spec = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[1], ylags = 1, yfun = fun),
					variance.model=list(dynamic=TRUE, model="sGARCH"))
			solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
			fit = starfit(spec, data = dx[1:521], solver = "strategy", solver.control = solver.control, n=10)
			sim  = starsim(fit, n.sim=5000)
			specx = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]],
							statevar = c("y","s")[1], ylags = 1, yfun = fun),
					variance.model=list(dynamic=TRUE, model="sGARCH"))
			setfixed(specx)<-as.list(coef(fit))
			path = starpath(specx, n.sim=5000, prereturns = tail(dx[1:521], 12), preresiduals = tail(residuals(fit),2),
					presigma = tail(sigma(fit), 2),
					custom.dist=list(name="sample", distfit=sim@simulation$residSim/sigma(sim)))
			test[i,j] = as.integer(all.equal(fitted(sim), fitted(path)))
			print(test[i,j])
		}
	}
	zz <- file("test2c-1.txt", open="wt")
	sink(zz)
	print(test)
	sink(type="message")
	sink()
	close(zz)
	
	
	# 2. Equivalence of sim and path methods for self-exciting model without yfun
	arorder=list(c(0,0),c(0,1), c(1,0), c(1,1), c(2,1))
	int = list(c(0,0), c(1,0))
	test = matrix(NA, 2, 5)
	for(i in 1:2){
		for(j in 1:5){
			spec = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[1], ylags = 1),
					variance.model=list(dynamic=TRUE, model="sGARCH"))
			# set reltol to a high number in order to converge quickly since
			solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-6, method="BFGS", parsearch=TRUE)
			fit = starfit(spec, data = dx[1:521], solver = "optim", solver.control = solver.control, n=5)
			sim  = starsim(fit, n.sim=5000)
			specx = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]],
							statevar = c("y","s")[1], ylags = 1),
					variance.model=list(dynamic=TRUE, model="sGARCH"))
			setfixed(specx)<-as.list(coef(fit))
			path = starpath(specx, n.sim=5000, prereturns = tail(dx[1:521], 6), preresiduals = tail(residuals(fit),2),
					presigma = tail(sigma(fit), 2),
					custom.dist=list(name="sample", distfit=sim@simulation$residSim/sigma(sim)))
			test[i,j] = as.integer(all.equal(fitted(sim), fitted(path)))
			print(test[i,j])
		}
	}
	zz <- file("test2c-2.txt", open="wt")
	sink(zz)
	print(test)
	sink(type="message")
	sink()
	close(zz)
	
	
	# 3. Equivalence of sim and path methods for s-var model
	dxv = runMean(abs(dx), n=4)
	dxv[1:3] = c(as.numeric(abs(dx[1,1])), as.numeric(mean(abs(dx[1:2,1]))), 
			as.numeric(mean(abs(dx[1:3,1]))))
	dxv = lag(dxv, 1)
	dxv[1] = 0
	arorder=list(c(0,0),c(0,1), c(1,0), c(1,1), c(2,1))
	int = list(c(0,0), c(1,0))
	test = matrix(NA, 2, 5)
	
	for(i in 1:2){
		for(j in 1:5){
			spec = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[2], s = dxv[1:521]),
					variance.model=list(dynamic=TRUE, model="sGARCH"))
			solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-8, method="BFGS", parsearch=TRUE)
			fit = starfit(spec, data = dx[1:521], solver = "strategy", solver.control = solver.control, n=5)
			ssim = sample(as.numeric(dxv), 5000, replace=TRUE)
			sim  = starsim(fit, n.sim=5000, ssim = list(matrix(ssim, ncol=1)))
			specx = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[2], s = dxv[500:521]),
					variance.model=list(dynamic=TRUE, model="sGARCH"))
			setfixed(specx)<-as.list(coef(fit))
			path = starpath(specx, n.sim=5000, prereturns = tail(dx[1:521], 6), 
					preresiduals = tail(residuals(fit),2),
					presigma = tail(sigma(fit), 2),
					custom.dist=list(name="sample", distfit=sim@simulation$residSim/sigma(sim)),
					ssim=list(matrix(ssim, ncol=1)))
			test[i,j] = as.integer(all.equal(fitted(sim), fitted(path)))
			print(test[i,j])
		}
	}
	zz <- file("test2c-3.txt", open="wt")
	sink(zz)
	print(test)
	sink(type="message")
	sink()
	close(zz)
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

twinkle.test4c = function(cluster = NULL){
	
	tic = Sys.time()
	require(twinkle)
	require(quantmod)
	data(forex)	
	
	fx = na.locf(forex, fromLast = TRUE)
	fxw = fx[which(weekdays(index(forex))=="Wednesday"),4]
	dx = ROC(fxw, na.pad=FALSE)*100
	fun = function(x){
		x = as.numeric(x)
		N = length(x)
		if(N<4){
			y = abs(x)	
		} else{
			y = runMean(abs(x), n=4)
			y[1:3] = c(abs(x[1]), mean(abs(x[1:2])), mean(abs(x[1:3])))
		}
		return(y)
	}
	
	spec = starspec(mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(1,1), 
					statevar = c("y","s")[1], ylags = 1, yfun = fun))
	solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
	fit = starfit(spec, data = dx[1:521], solver = "strategy", solver.control = solver.control, n=5)
	save(fit, file="/Users/alexiosg/Development/Twinkle/simulation/fit.rda")
	
	cf = coef(fit)
	
	sim  = starsim(fit, n.sim=300, m.sim=100)
	
	S300 = fitted(sim)
	# N = 300
	S300 = xts(S300, as.Date(1:300))
	save(S300, file="/Users/alexiosg/Development/Twinkle/simulation/S300.rda")
	if(!is.null(cluster)){
		clusterEvalQ(cluster, library(twinkle))
		clusterEvalQ(cluster, library(quantmod))
		clusterExport(cluster, c("S300","fun"), envir = environment())
		X300 = parLapply(cluster, 1:100, function(i){
					spec = starspec(mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(1,1), 
									statevar = c("y","s")[1], ylags = 1, yfun = fun))
					setbounds(spec)<-list(s1.phi0=c(-5, 0), s2.phi0=c(0, 5))
					ctrl1=list(maxit=5000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
					ctrl2=list(maxit=5000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE, n.restarts=3)
					fit1 = try(starfit(spec, data = S300[,i], solver = "strategy", solver.control = ctrl1, n=6), silent=TRUE)
					fit2 = try(starfit(spec, data = S300[,i], solver = "msoptim", solver.control = ctrl2), silent=TRUE)
					if(inherits(fit1, 'try-error') || convergence(fit1)!=0) cf1 = rep(NA, 8) else cf1 = c(coef(fit1), likelihood(fit1))
					if(inherits(fit2, 'try-error') || convergence(fit2)!=0) cf2 = rep(NA, 8) else cf2 = c(coef(fit2), likelihood(fit2))
					cfx = list(strategy = cf1, msoptim = cf2)
					return(cfx)
				})
	} else{
		X300 = lapply(1:100, function(i){
					spec = starspec(mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(1,1), 
									statevar = c("y","s")[1], ylags = 1, yfun = fun))
					setbounds(spec)<-list(s1.phi0=c(-5, 0), s2.phi0=c(0, 5))
					ctrl1=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
					ctrl2=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE, n.restarts=3)
					fit1 = starfit(spec, data = S300[,i], solver = "strategy", solver.control = ctrl1, n=6)
					fit2 = starfit(spec, data = S300[,i], solver = "msoptim", solver.control = ctrl2)
					cfx = list(strategy = c(coef(fit1), likelihood(fit1)), msoptim = c(coef(fit2), likelihood(fit2)))
					return(cfx)
				})
	}
	save(X300, file="/Users/alexiosg/Development/Twinkle/simulation/X300.rda")
	
	cf1 = t(sapply(X300, function(x) x$strategy))
	cf2 = t(sapply(X300, function(x) x$msoptim))
	colSd = function(x, na.rm=FALSE) apply(x, 2, function(x) sd(x, na.rm=na.rm))
	cbind(colMeans(cf1, na.rm=TRUE)[1:7], colMeans(cf2, na.rm=TRUE)[1:7], cf)
	cbind(colSd(cf1, na.rm=TRUE)[1:7], colSd(cf2, na.rm=TRUE)[1:7], fit@fit$matcoef[,2])
	
	sim  = starsim(fit, n.sim=600, m.sim=100)
	S600 = fitted(sim)
	# N = 600
	S600 = xts(S600, as.Date(1:600))
	save(S600, file="/Users/alexiosg/Development/Twinkle/simulation/S600.rda")
	if(!is.null(cluster)){
		clusterEvalQ(cluster, library(twinkle))
		clusterEvalQ(cluster, library(quantmod))
		clusterExport(cluster, c("S600","fun"), envir = environment())
		X600 = parLapply(cluster, 1:100, function(i){
					spec = starspec(mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(1,1), 
									statevar = c("y","s")[1], ylags = 1, yfun = fun))
					setbounds(spec)<-list(s1.phi0=c(-5, 0), s2.phi0=c(0, 5))
					ctrl1=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
					ctrl2=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE, n.restarts=3)
					fit1 = try(starfit(spec, data = S600[,i], solver = "strategy", solver.control = ctrl1, n=6), silent=TRUE)
					fit2 = try(starfit(spec, data = S600[,i], solver = "msoptim", solver.control = ctrl2), silent=TRUE)
					if(inherits(fit1, 'try-error') || convergence(fit1)!=0) cf1 = rep(NA, 8) else cf1 = c(coef(fit1), likelihood(fit1))
					if(inherits(fit2, 'try-error') || convergence(fit2)!=0) cf2 = rep(NA, 8) else cf2 = c(coef(fit2), likelihood(fit2))
					cfx = list(strategy = cf1, msoptim = cf2)
					return(cfx)
				})
	} else{
		X600 = lapply(1:100, function(i){
					spec = starspec(mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(1,1), 
									statevar = c("y","s")[1], ylags = 1, yfun = fun))
					setbounds(spec)<-list(s1.phi0=c(-5, 0), s2.phi0=c(0, 5))
					ctrl1=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
					ctrl2=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE, n.restarts=3)
					fit1 = starfit(spec, data = S600[,i], solver = "strategy", solver.control = ctrl1, n=6)
					fit2 = starfit(spec, data = S600[,i], solver = "msoptim", solver.control = ctrl2)
					cfx = list(strategy = c(coef(fit1), likelihood(fit1)), msoptim = c(coef(fit2), likelihood(fit2)))
					return(cfx)
				})
	}
	
	save(X600, file="/Users/alexiosg/Development/Twinkle/simulation/X600.rda")
	
	sim  = starsim(fit, n.sim=900, m.sim=100)
	S900 = fitted(sim)
	# N = 600
	S900 = xts(S900, as.Date(1:900))
	save(S900, file="/Users/alexiosg/Development/Twinkle/simulation/S900.rda")
	if(!is.null(cluster)){
		clusterEvalQ(cluster, library(twinkle))
		clusterEvalQ(cluster, library(quantmod))
		clusterExport(cluster, c("S900","fun"), envir = environment())
		X900 = parLapply(cluster, 1:100, function(i){
					spec = starspec(mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(1,1), 
									statevar = c("y","s")[1], ylags = 1, yfun = fun))
					setbounds(spec)<-list(s1.phi0=c(-5, 0), s2.phi0=c(0, 5))
					ctrl1=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
					ctrl2=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE, n.restarts=3)
					fit1 = try(starfit(spec, data = S900[,i], solver = "strategy", solver.control = ctrl1, n=6), silent=TRUE)
					fit2 = try(starfit(spec, data = S900[,i], solver = "msoptim", solver.control = ctrl2), silent=TRUE)
					if(inherits(fit1, 'try-error') || convergence(fit1)!=0) cf1 = rep(NA, 8) else cf1 = c(coef(fit1), likelihood(fit1))
					if(inherits(fit2, 'try-error') || convergence(fit2)!=0) cf2 = rep(NA, 8) else cf2 = c(coef(fit2), likelihood(fit2))
					cfx = list(strategy = cf1, msoptim = cf2)
					return(cfx)
				})
	} else{
		X900 = lapply(1:100, function(i){
					spec = starspec(mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(1,1), 
									statevar = c("y","s")[1], ylags = 1, yfun = fun))
					setbounds(spec)<-list(s1.phi0=c(-5, 0), s2.phi0=c(0, 5))
					ctrl1=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
					ctrl2=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE, n.restarts=3)
					fit1 = starfit(spec, data = S900[,i], solver = "strategy", solver.control = ctrl1, n=6)
					fit2 = starfit(spec, data = S900[,i], solver = "msoptim", solver.control = ctrl2)
					cfx = list(strategy = c(coef(fit1), likelihood(fit1)), msoptim = c(coef(fit2), likelihood(fit2)))
					return(cfx)
				})
	}
	
	save(X900, file="/Users/alexiosg/Development/Twinkle/simulation/X900.rda")
	
	
	sim  = starsim(fit, n.sim=1200, m.sim=100)
	S1200 = fitted(sim)
	# N = 1200
	S1200 = xts(S1200, as.Date(1:1200))
	save(S1200, file="/Users/alexiosg/Development/Twinkle/simulation/S1200.rda")
	if(!is.null(cluster)){
		clusterEvalQ(cluster, library(twinkle))
		clusterEvalQ(cluster, library(quantmod))
		clusterExport(cluster, c("S1200","fun"), envir = environment())
		X1200 = parLapply(cluster, 1:100, function(i){
					spec = starspec(mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(1,1), 
									statevar = c("y","s")[1], ylags = 1, yfun = fun))
					setbounds(spec)<-list(s1.phi0=c(-5, 0), s2.phi0=c(0, 5))
					ctrl1=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
					ctrl2=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE, n.restarts=3)
					fit1 = try(starfit(spec, data = S1200[,i], solver = "strategy", solver.control = ctrl1, n=6), silent=TRUE)
					fit2 = try(starfit(spec, data = S1200[,i], solver = "msoptim", solver.control = ctrl2), silent=TRUE)
					if(inherits(fit1, 'try-error') || convergence(fit1)!=0) cf1 = rep(NA, 8) else cf1 = c(coef(fit1), likelihood(fit1))
					if(inherits(fit2, 'try-error') || convergence(fit2)!=0) cf2 = rep(NA, 8) else cf2 = c(coef(fit2), likelihood(fit2))
					cfx = list(strategy = cf1, msoptim = cf2)
					return(cfx)
				})
	} else{
		X1200 = lapply(1:100, function(i){
					spec = starspec(mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(1,1), 
									statevar = c("y","s")[1], ylags = 1, yfun = fun))
					setbounds(spec)<-list(s1.phi0=c(-5, 0), s2.phi0=c(0, 5))
					ctrl1=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
					ctrl2=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE, n.restarts=3)
					fit1 = starfit(spec, data = S1200[,i], solver = "strategy", solver.control = ctrl1, n=6)
					fit2 = starfit(spec, data = S1200[,i], solver = "msoptim", solver.control = ctrl2)
					cfx = list(strategy = c(coef(fit1), likelihood(fit1)), msoptim = c(coef(fit2), likelihood(fit2)))
					return(cfx)
				})
	}
	
	save(X1200, file="/Users/alexiosg/Development/Twinkle/simulation/X1200.rda")
	
	
	
	sim  = starsim(fit, n.sim=1500, m.sim=100)
	S1500 = fitted(sim)
	# N = 1500
	S1500 = xts(S1500, as.Date(1:1500))
	save(S1500, file="/Users/alexiosg/Development/Twinkle/simulation/S1500.rda")
	if(!is.null(cluster)){
		clusterEvalQ(cluster, library(twinkle))
		clusterEvalQ(cluster, library(quantmod))
		clusterExport(cluster, c("S1500","fun"), envir = environment())
		X1500 = parLapply(cluster, 1:100, function(i){
					spec = starspec(mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(1,1), 
									statevar = c("y","s")[1], ylags = 1, yfun = fun))
					setbounds(spec)<-list(s1.phi0=c(-5, 0), s2.phi0=c(0, 5))
					ctrl1=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
					ctrl2=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE, n.restarts=3)
					fit1 = try(starfit(spec, data = S1500[,i], solver = "strategy", solver.control = ctrl1, n=6), silent=TRUE)
					fit2 = try(starfit(spec, data = S1500[,i], solver = "msoptim", solver.control = ctrl2), silent=TRUE)
					if(inherits(fit1, 'try-error') || convergence(fit1)!=0) cf1 = rep(NA, 8) else cf1 = c(coef(fit1), likelihood(fit1))
					if(inherits(fit2, 'try-error') || convergence(fit2)!=0) cf2 = rep(NA, 8) else cf2 = c(coef(fit2), likelihood(fit2))
					cfx = list(strategy = cf1, msoptim = cf2)
					return(cfx)
				})
	} else{
		X1500 = lapply(1:100, function(i){
					spec = starspec(mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(1,1), 
									statevar = c("y","s")[1], ylags = 1, yfun = fun))
					setbounds(spec)<-list(s1.phi0=c(-5, 0), s2.phi0=c(0, 5))
					ctrl1=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
					ctrl2=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE, n.restarts=3)
					fit1 = starfit(spec, data = S1500[,i], solver = "strategy", solver.control = ctrl1, n=6)
					fit2 = starfit(spec, data = S1500[,i], solver = "msoptim", solver.control = ctrl2)
					cfx = list(strategy = c(coef(fit1), likelihood(fit1)), msoptim = c(coef(fit2), likelihood(fit2)))
					return(cfx)
				})
	}
	
	save(X1500, file="/Users/alexiosg/Development/Twinkle/simulation/X1500.rda")

	load(file="/Users/alexiosg/Development/Twinkle/simulation/X300.rda")
	load(file="/Users/alexiosg/Development/Twinkle/simulation/X600.rda")
	load(file="/Users/alexiosg/Development/Twinkle/simulation/X900.rda")
	load(file="/Users/alexiosg/Development/Twinkle/simulation/X1200.rda")
	load(file="/Users/alexiosg/Development/Twinkle/simulation/X1500.rda")
	
	
	.rmse = function(est, act)
	{
		exc = unique(which(is.na(est), arr.ind = TRUE)[,1])
		if(length(exc)>0) est = est[-exc, , drop = FALSE]
		n = dim(est)[1]
		diff = est - repmat(t(act), n, 1)
		ans = apply(diff, 2, FUN = function(x) sum((x)^2)/(n-1) )
		return(sqrt(ans))
	}
	
	cnames = names(cf)
	
	cf["s1.c"]=cf["s1.c"]/cf["s1.alpha1"]
	for(i in 1:100){
		if(!is.na(X300[[i]]$strategy["s1.c"])) X300[[i]]$strategy["s1.c"] = X300[[i]]$strategy["s1.c"]/X300[[i]]$strategy["s1.alpha1"]
		if(!is.na(X300[[i]]$msoptim["s1.c"])) X300[[i]]$msoptim["s1.c"] = X300[[i]]$msoptim["s1.c"]/X300[[i]]$msoptim["s1.alpha1"]
		if(!is.na(X600[[i]]$strategy["s1.c"])) X600[[i]]$strategy["s1.c"] = X600[[i]]$strategy["s1.c"]/X600[[i]]$strategy["s1.alpha1"]
		if(!is.na(X600[[i]]$msoptim["s1.c"])) X600[[i]]$msoptim["s1.c"] = X600[[i]]$msoptim["s1.c"]/X600[[i]]$msoptim["s1.alpha1"]
		if(!is.na(X900[[i]]$strategy["s1.c"])) X900[[i]]$strategy["s1.c"] = X900[[i]]$strategy["s1.c"]/X900[[i]]$strategy["s1.alpha1"]
		if(!is.na(X900[[i]]$msoptim["s1.c"])) X900[[i]]$msoptim["s1.c"] = X900[[i]]$msoptim["s1.c"]/X900[[i]]$msoptim["s1.alpha1"]
		if(!is.na(X1200[[i]]$strategy["s1.c"])) X1200[[i]]$strategy["s1.c"] = X1200[[i]]$strategy["s1.c"]/X1200[[i]]$strategy["s1.alpha1"]
		if(!is.na(X1200[[i]]$msoptim["s1.c"])) X1200[[i]]$msoptim["s1.c"] = X1200[[i]]$msoptim["s1.c"]/X1200[[i]]$msoptim["s1.alpha1"]
		if(!is.na(X1500[[i]]$strategy["s1.c"])) X1500[[i]]$strategy["s1.c"] = X1500[[i]]$strategy["s1.c"]/X1500[[i]]$strategy["s1.alpha1"]
		if(!is.na(X1500[[i]]$msoptim["s1.c"])) X1500[[i]]$msoptim["s1.c"] = X1500[[i]]$msoptim["s1.c"]/X1500[[i]]$msoptim["s1.alpha1"]
	}
	
	strategyrmse = rbind(
			rugarch:::.rmse(t(sapply(X300, function(x) x$strategy))[,-8], cf),
			rugarch:::.rmse(t(sapply(X600, function(x) x$strategy))[,-8], cf),
			rugarch:::.rmse(t(sapply(X900, function(x) x$strategy))[,-8], cf),
			rugarch:::.rmse(t(sapply(X1200, function(x) x$strategy))[,-8], cf),
			rugarch:::.rmse(t(sapply(X1500, function(x) x$strategy))[,-8], cf))
	
	multistartrmse = rbind(
			rugarch:::.rmse(t(sapply(X300, function(x) x$msoptim))[,-8], cf),
			rugarch:::.rmse(t(sapply(X600, function(x) x$msoptim))[,-8], cf),
			rugarch:::.rmse(t(sapply(X900, function(x) x$msoptim))[,-8], cf),
			rugarch:::.rmse(t(sapply(X1200, function(x) x$msoptim))[,-8], cf),
			rugarch:::.rmse(t(sapply(X1500, function(x) x$msoptim))[,-8], cf))
	
	s1.phi0 = 
	
	s1.phi1 = cbind(t(sapply(X300, function(x) x$strategy))[,1],
			t(sapply(X600, function(x) x$strategy))[,1],
			t(sapply(X900, function(x) x$strategy))[,1],
			t(sapply(X1200, function(x) x$strategy))[,1],
			t(sapply(X1500, function(x) x$strategy))[,1])
	
	mx = list()
	mx[[1]] = expression(phi['1,0'])
	mx[[2]] = expression(phi['1,1'])
	mx[[3]] = expression(phi['2,0'])
	mx[[4]] = expression(phi['2,1'])
	mx[[5]] = expression(c)
	mx[[7]] = expression(sigma)
	
	par(mfrow=c(2,3))	
	for(i in c(1,2,3,4,5,7)){
		tmp = cbind(t(sapply(X300, function(x) x$strategy))[,i],
				t(sapply(X600, function(x) x$strategy))[,i],
				t(sapply(X900, function(x) x$strategy))[,i],
				t(sapply(X1200, function(x) x$strategy))[,i],
				t(sapply(X1500, function(x) x$strategy))[,i])
		boxplot(tmp, names=c("300","600","900","1200","1500"), col = "bisque")
		abline(h=cf[i], col=2, lty=2)
		title(main = mx[[i]])
	}
	
	
	
		
}