/* 
Estimation of Nonparametric instrumental variable (NPIV) models
Version 0.2.0 30th May 2016
This program estimates the function g(x) in
Y = g(X) + e with E(e|Z)=0
where Y is a scalar dependent variable ("depvar"), 
X is a scalar endogenous variable ("expvar"), and 
Z a scalar instrument ("inst").
Syntax:
npivreg depvar expvar inst
*/

program define npivreg
		version 12
		
		// initializations
		syntax varlist(numeric)
		display "varlist is `varlist'"
		
		// defining temporary names (disappear after running command)
		tempname b p Yhat
		tempvar beta P
		
		capture ssc install bspline
		capture drop basisexpvar* basisinst* npest*
		tempname 
		global mylist `varlist'
		global depvar   : word 1 of $mylist
		global expvar   : word 2 of $mylist
		global inst     : word 3 of $mylist
		
		// generate B-spline bases for X and Z
		quietly bspline, xvar($expvar) gen(basisexpvar) power(2)
		quietly bspline, xvar($inst) gen(basisinst) power(3)
		
		// compute NPIV fitted value
		mata : P 		= st_data(., "basisexpvar*")
		mata : Q 		= st_data(., "basisinst*")
		mata : Y 		= st_data(., "$depvar")
		mata : MQ 		= Q*invsym(Q'*Q)*Q'
		mata : b  		= invsym(P'*MQ*P)*P'*MQ*Y
		mata : Yhat 	= P*b
		
		mata : st_matrix("`b'", b)
		mata : st_matrix("`p'", P)
		mata : st_matrix("`Yhat'", Yhat)

		// plot NPIV fitted value
		svmat `Yhat', name(npest) // compute and return the fitted value with name "npest"
		svmat `p', name(`P')
		svmat `b', name(`beta')
		label variable npest "NPIV fitted value"
		scatter $depvar $expvar, msym(circle_hollow) || line npest $expvar, sort	
		drop basisexpvar* basisinst*
end
