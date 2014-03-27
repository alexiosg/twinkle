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
# asymptotically chisquared (m ��������� p ��������� q) distribution, assuming that 
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
