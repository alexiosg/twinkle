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

# McLeod and Li Test of Heteroscedasticity:
# asymptotically chisquared distribution, assuming that 
# lags/n is small and lags is moderately large

mclitest = function(x, lag.max = 10, p = 0, q = 0)
{
	n = length(x)
	y = acf(x^2, lag.max = lag.max, plot  = FALSE)$acf^2
	STATISTIC = n*(n+2)*cumsum( (1/(n-(1:lag.max))) * y )
	pvalues = 1-pchisq(STATISTIC, 1:lag.max - p - q)
	return(list(STATISTIC = STATISTIC, pvalues = pvalues))
}

bgstartest = function(object, method = "Chisq", lag.max = 10)
{
	z = score(object)
	r = as.numeric(residuals(object))
	n = nrow(z)
	e = matrix(0, ncol = lag.max, nrow = n)
	for(i in 1:lag.max){e[(i+1):n,i] = r[1:(n-i)]}
	z = cbind(z, e)
	z = z[-c(1:lag.max),,drop=FALSE]
	r = r[-c(1:lag.max)]
	Ru = residuals(lm(r~z-1))
	nq = ncol(z)
	T = nrow(z)
	if(method=="F"){
		# F-Statistic ~ F(df1=lag.max, df2=T-nq-lag.max) 
		stat = ((T-nq-lag.max)/lag.max)*(sum(r^2) - sum(Ru^2))/(sum(Ru^2))
		dof = c(lag.max, T-nq-lag.max)
		pvalue = pf(stat, df1 = dof[1], df2 = dof[2], lower.tail=FALSE)
		names(dof)=c("df1","df2")
	} else if(method=="Chisq"){
		# LM Statistic ~ chisq(lag.max)
		stat = T*(sum(r^2) - sum(Ru^2))/(sum(Ru^2))
		dof = lag.max
		pvalue = pchisq(stat, dof, lower.tail = FALSE)
		names(dof)="df"
	} else{
		stop("\nunrecognized method.\n")
	}
	names(stat)<-"LM test"
	names(pvalue)<-"p.value"
	test = list(statistic = stat, parameter = dof,
			method = paste("LM Test (",method,") for residual serial correlation in STAR model",sep=""), 
			p.value = pvalue,
			data.name = "STARfit object", coeffficients = NULL)
	class(test)<-"htest"
	return(test)
}


# LM linearity testing against 2 regime STAR
#
#   Performs an 3rd order Taylor expansion LM test
#
#   str: an nlar.struct object
#   rob
#   sig
linearityTest.star <- function(str, thVar, externThVar=FALSE,
		rob=FALSE, sig=0.05, trace=TRUE, ...)
{
	n.used <- NROW(str$xx);  # The number of lagged samples
	
	# Build the regressand vector
	y_t <- str$yy;
	
	# Regressors under the null
	xH0 <- cbind(1, str$xx)
	
	# Get the transition variable
	s_t <- thVar
	
	# Linear Model (null hypothesis)
	linearModel <- lm(y_t ~ . , data=data.frame(xH0))
	
	u_t <- linearModel$residuals;
	SSE0 <- sum(u_t^2)
	
	# Regressors under the alternative
	if (externThVar) {
#    if (rob) {} else {
		tmp <- rep(s_t, NCOL(str$xx) + 1)
		dim(tmp) <- c(length(s_t), NCOL(str$xx) + 1)
		xH1 <- cbind(cbind(1, str$xx) * tmp, cbind(1, str$xx) * (tmp^2),
				cbind(1, str$xx) * (tmp^3))
	} else {
#    if (rob) {} else {
		tmp <- rep(s_t, NCOL(str$xx))
		dim(tmp) <- c(length(s_t), NCOL(str$xx))
		xH1 <- cbind(str$xx * tmp, str$xx * (tmp^2), str$xx * (tmp^3))
	}
	
	# Standarize the regressors
	Z <- cbind(xH0, xH1);
	nZ <- NCOL(Z);
	sdZ <- apply(Z,2,sd)
	dim(sdZ) <- c(1, nZ)
	sdZ <- kronecker(matrix(1,n.used,1), sdZ) # repeat sdZ n.used rows
	Z[,2:nZ] <- Z[,2:nZ] / sdZ[,2:nZ]
	
	# Nonlinear model (alternative hypothesis)
	nonlinearModel <- lm(u_t ~ ., data=data.frame(Z));
	e_t <- nonlinearModel$residuals;
	SSE1 <- sum(e_t^2)
	
	# Compute the test statistic
	n <- NCOL(xH0);
	m <- NCOL(xH1);
	
	F = ((SSE0 - SSE1) / m) / (SSE1 / (n.used - m - n));
	
	# Look up the statistic in the table, get the p-value
	pValue <- pf(F, m, n.used - m - n, lower.tail = FALSE);
	
	if (pValue >= sig) {
		return(list(isLinear = TRUE, pValue = pValue));
	}
	else {
		return(list(isLinear = FALSE, pValue = pValue));
	}
	
}