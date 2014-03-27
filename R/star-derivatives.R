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


.stars2sLLHgrad1 = function(pars, arglist)
{
	modelinc = arglist$model$modelinc
	data = arglist$data
	returnType = arglist$returnType
	starenv = arglist$starenv
	# rejoin fixed and pars
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars	
	ipars[estidx, 1] = pars
	if(arglist$transform==1){
		if(modelinc[1]>0 && ipars["s1.phi0",4]==1) ipars["s1.phi0",1] = logtransform(ipars["s1.phi0",1], ipars["s1.phi0","LB"], ipars["s1.phi0","UB"])
		if(modelinc[4]>0 && ipars["s2.phi0",4]==1) ipars["s2.phi0",1] = logtransform(ipars["s2.phi0",1], ipars["s2.phi0","LB"], ipars["s2.phi0","UB"])
		if(modelinc[2]>0 && ipars["s1.phi1",4]==1) ipars["s1.phi1",1] = logtransform(ipars["s1.phi1",1], ipars["s1.phi1","LB"], ipars["s1.phi1","UB"])
		if(modelinc[5]>0 && ipars["s2.phi1",4]==1) ipars["s2.phi1",1] = logtransform(ipars["s2.phi1",1], ipars["s2.phi1","LB"], ipars["s2.phi1","UB"])
		if(modelinc[22]>0 && ipars["sigma",4]==1) ipars["sigma",1] = logtransform(ipars["sigma",1], ipars["sigma","LB"], ipars["sigma","UB"])
		if(modelinc[33]>0 && ipars["skew",4]==1) ipars["skew",1] = logtransform(ipars["skew",1], ipars["skew","LB"], ipars["skew","UB"])
		if(modelinc[34]>0 && ipars["shape",4]==1) ipars["shape",1] = logtransform(ipars["shape",1], ipars["shape","LB"], ipars["shape","UB"])
		if(modelinc[35]>0 && ipars["ghlambda",4]==1) ipars["ghlambda",1] = logtransform(ipars["ghlambda",1], ipars["ghlambda","LB"], ipars["ghlambda","UB"])
		if(modelinc[15]>0 && ipars["s1.beta",4]==1) ipars[idx["s1.beta",1],1] = logtransform(ipars[idx["s1.beta",1],1], -0.99, 0.99)
	}
	if(modelinc[39]==1){
		probs = arglist$probs
	} else{
		xprobs = dstar2(ipars, arglist)
		probs = xprobs$probs
		pmu = xprobs$pmu
		initp = xprobs$initp
	}
	trace = arglist$trace
	T = length(data)
	mexdata = arglist$mexdata
	fit.control = arglist$fit.control	
	m = model$maxOrder
	N = c(m,T)
	distribution = model$modeldesc$distribution
	# distribution number
	# 1 = norm, 2=snorm, 3=std, 4=sstd, 5=ged, 6=sged, 7=nig
	dist = model$modeldesc$distno
	if(modelinc[3]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)
	
	if(modelinc[1]>0) dm1 = rep(0, T) else dm1 = 0
	if(modelinc[2]>0) dphi1 = rep(0, T*modelinc[2]) else dphi1 = 0
	if(modelinc[3]>0) dxi1 = rep(0, T*modelinc[3]) else dxi1 = 0
	if(modelinc[4]>0) dm2 = rep(0, T) else dm2 = 0
	if(modelinc[5]>0) dphi2 = rep(0, T*modelinc[5]) else dphi2 = 0
	if(modelinc[6]>0) dxi2 = rep(0, T*modelinc[6]) else dxi2 = 0
	if(modelinc[13]>0) dc = rep(0, T) else dc = 0
	if(modelinc[14]>0) dalpha = rep(0, T*modelinc[14]) else dalpha = 0
	if(modelinc[15]>0) dbeta = rep(0, T) else dbeta = 0
	dsigma = rep(0, T)
	tmp = try( .C("stars2normalderiv1", 
					model = as.integer(modelinc[1:50]), 
					pars = as.double(ipars[,1]), 
					dpars = as.double(ipars[,1]*0),
					idx = as.integer(idx[,1]-1), 
					x = as.double(data), 
					res = double(T),  
					mexdata = mexdata, 
					prob = as.double(probs), 
					y = as.double(arglist$XL),
					mpt = as.double(pmu),
					mpinit = as.double(initp),
					constm = double(2*T), 
					condm = double(2*T), 
					m = as.integer(m), 
					T = as.integer(T), 
					dm1 = as.double(dm1), 
					dm2 = as.double(dm2), 
					dxi1 = as.double(dxi1),
					dxi2 = as.double(dxi2),
					dphi1 = as.double(dphi1),
					dphi2 = as.double(dphi2),
					dc = as.double(dc),
					dalpha = as.double(dalpha),
					dbeta = as.double(dbeta),
					dsigma = as.double(dsigma),
					PACKAGE = "twinkle"), silent = TRUE )
	if(arglist$derivarg == 1){
		if(inherits(tmp, 'try-error')){
			ans = rep(1e5, length(pnames))
		} else{
			g = NULL
			if(modelinc[1]>0) g = c(g, tmp$dpars[idx["s1.mu",1]])
			if(modelinc[2]>0) g = c(g, tmp$dpars[idx["s1.phi",1]:idx["s1.phi",2]])
			if(modelinc[3]>0) g = c(g, tmp$dpars[idx["s1.xi",1]:idx["s1.xi",2]])
			if(modelinc[4]>0) g = c(g, tmp$dpars[idx["s2.mu",1]])
			if(modelinc[5]>0) g = c(g, tmp$dpars[idx["s2.phi",1]:idx["s2.phi",2]])
			if(modelinc[6]>0) g = c(g, tmp$dpars[idx["s2.xi",1]:idx["s2.xi",2]])
			if(modelinc[13]>0) g = c(g, tmp$dpars[idx["s1.c",1]])
			if(modelinc[14]>0) g = c(g, tmp$dpars[idx["s1.alpha",1]:idx["s1.alpha",2]])
			if(modelinc[15]>0) g = c(g, tmp$dpars[idx["s1.beta",1]])
			g = c(g, tmp$dpars[idx["sigma",1]])
			g = -1*g
			names(g)<-pnames
		}
	} else{
		if(!inherits(tmp, 'try-error')){
			g = matrix(1, nrow=T)
			if(modelinc[1]>0) g = cbind(g, tmp$dm1)
			if(modelinc[2]>0) g = cbind(g, matrix(tmp$dphi1, ncol = modelinc[2]))
			if(modelinc[3]>0) g = cbind(g, matrix(tmp$dxi1, ncol = modelinc[3]))
			if(modelinc[4]>0) g = cbind(g, tmp$dm2)
			if(modelinc[5]>0) g = cbind(g, matrix(tmp$dphi2, ncol = modelinc[5]))
			if(modelinc[6]>0) g = cbind(g, matrix(tmp$dxi2, ncol = modelinc[6]))
			if(modelinc[13]>0) g = cbind(g, tmp$dc)
			if(modelinc[14]>0) g = cbind(g, matrix(tmp$dalpha, ncol = modelinc[14]))
			if(modelinc[15]>0) g = cbind(g, tmp$dbeta)
			g = cbind(g, tmp$dsigma)
			g = g[,-1]
			g = -1*g
			colnames(g)<-pnames
		} else{
			g = tmp
		}
	}
	return(g)
}
.stars2sLLHgrad2 = function(pars, arglist)
{
	modelinc = arglist$model$modelinc
	data = arglist$data
	returnType = arglist$returnType
	starenv = arglist$starenv
	# rejoin fixed and pars
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars	
	ipars[estidx, 1] = pars
	if(arglist$transform==1){
		if(modelinc[1]>0 && ipars["s1.phi0",4]==1) ipars["s1.phi0",1] = logtransform(ipars["s1.phi0",1], ipars["s1.phi0","LB"], ipars["s1.phi0","UB"])
		if(modelinc[4]>0 && ipars["s2.phi0",4]==1) ipars["s2.phi0",1] = logtransform(ipars["s2.phi0",1], ipars["s2.phi0","LB"], ipars["s2.phi0","UB"])
		if(modelinc[2]>0 && ipars["s1.phi1",4]==1) ipars["s1.phi1",1] = logtransform(ipars["s1.phi1",1], ipars["s1.phi1","LB"], ipars["s1.phi1","UB"])
		if(modelinc[5]>0 && ipars["s2.phi1",4]==1) ipars["s2.phi1",1] = logtransform(ipars["s2.phi1",1], ipars["s2.phi1","LB"], ipars["s2.phi1","UB"])
		if(modelinc[22]>0 && ipars["sigma",4]==1) ipars["sigma",1] = logtransform(ipars["sigma",1], ipars["sigma","LB"], ipars["sigma","UB"])
		if(modelinc[33]>0 && ipars["skew",4]==1) ipars["skew",1] = logtransform(ipars["skew",1], ipars["skew","LB"], ipars["skew","UB"])
		if(modelinc[34]>0 && ipars["shape",4]==1) ipars["shape",1] = logtransform(ipars["shape",1], ipars["shape","LB"], ipars["shape","UB"])
		if(modelinc[35]>0 && ipars["ghlambda",4]==1) ipars["ghlambda",1] = logtransform(ipars["ghlambda",1], ipars["ghlambda","LB"], ipars["ghlambda","UB"])
		if(modelinc[15]>0 && ipars["s1.beta",4]==1) ipars[idx["s1.beta",1],1] = logtransform(ipars[idx["s1.beta",1],1], -0.99, 0.99)
	}
	if(modelinc[39]==1){
		probs = arglist$probs
		pi0=0
		XL = NULL
		meanXL = 0
	} else{
		xprobs = dstar2(ipars, arglist)
		probs = xprobs$probs
		pmu = xprobs$pmu
		pi0 = xprobs$initp
		XL = arglist$XL
		meanXL = colMeans(XL)
	}
	trace = arglist$trace
	T = length(data)
	mexdata = arglist$mexdata
	fit.control = arglist$fit.control	
	m = model$maxOrder
	N = c(m,T)
	distribution = model$modeldesc$distribution
	# distribution number
	# 1 = norm, 2=snorm, 3=std, 4=sstd, 5=ged, 6=sged, 7=nig
	dist = model$modeldesc$distno
	if(modelinc[3]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)
	
	alpha = ipars[idx["s1.alpha",1]:idx["s1.alpha",2],1]
	beta = ipars[idx["s1.beta",1]:idx["s1.beta",2],1]
	
	dpidalpha = matrix(0, ncol = ncol(XL), nrow = nrow(XL))
	dpi0dalpha = rep(0, ncol(XL))
	for(i in 1:ncol(XL)) dpi0alpha[i] = meanXL[i]/(1-beta)
	dpidbeta = matrix(0, ncol = 1, nrow = nrow(XL))
	dpi0dbeta = (ipars["s1.c",1]+sum(alpha*meanXL))/(1-beta)^2
	dpidc = matrix(0, ncol = 1, nrow = nrow(XL))
	dpi0dc = 1/(1-beta)
	
	dLdalpha = matrix(0, ncol = ncol(XL), nrow=T)
	dLdc = dLdbeta = rep(0, T)
	dLdphi20 = dLdphi10 = rep(0, T)
	dLdphi11 = dLdphi22 = matrix(0, ncol = modelinc[2])
	dLdksi11 = dLdksi22 = matrix(0, ncol = modelinc[3])
	dLdsigma = rep(0, T)
	dL0dalpha = 0
	dL0dbeta  = 0
	dL0dc = 0
	dL0dphi10 =  dL0dphi20 = 0
	dL0dphi11 =  dL0dphi22 = rep(0, modelinc[2])
	dL0dksi11 = dL0dksi22 = rep(0, modelinc[3])
	dL0dsigma = 0
	for(i in 1:T){
		
		part1 = -e[i]/sigma^2
		part2 = exp(-pi[i])
		part3 = (1+part2)^2
		s1mu = sum(ipars[idx["s1.xi",1]:idx["s1.xi",2],1]*mexdata[i,]) + ipars["s1.phi0",1]
		s2mu = sum(ipars[idx["s2.xi",1]:idx["s2.xi",2],1]*mexdata[i,]) + ipars["s2.phi0",1]
		if(modelinc[2]>0){
			for(j in 1:modelinc[2]){
				if(j>i){
					s1mu = s1mu +  ipars[idx["s1.phi",1]+j-1,1]*x[i-j]
					s2mu = s2mu +  ipars[idx["s2.phi",1]+j-1,1]*x[i-j]
				}
			}
		}
		for(j in 1:ncol(XL)){
			leftpart = 
			dLdalpha[i,j] = ipars[idx["s1.alpha",1]+j-1,1]
		}
	}
	
}


