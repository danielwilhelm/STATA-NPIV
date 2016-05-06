/* 
Estimation of Nonparametric instrumental variable (NPIV) models
Version 0.1.0 6th May 2016

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
		
		capture ssc install bspline
		capture drop px* qx* npest* P*
		global mylist `varlist'
		global depvar   : word 1 of $mylist
		global expvar   : word 2 of $mylist
		global inst     : word 3 of $mylist
		
		// generate B-spline bases for X and Z
		quietly bspline, xvar($expvar) gen(px) power(2)
		quietly bspline, xvar($inst) gen(qx) power(3)
		
		// compute NPIV fitted value
		mata : P 		= st_data(., "px*")
		mata : Q 		= st_data(., "qx*")
		mata : Y 		= st_data(., "$depvar")
		mata : MQ 		= Q*invsym(Q'*Q)*Q'
		mata : b  		= invsym(P'*MQ*P)*P'*MQ*Y
		mata : Yhat 	= P*b
		
		mata : st_matrix("beta", b)
		mata : st_matrix("P", P)
		mata : st_matrix("npest", Yhat)

		// plot NPIV fitted value
		svmat npest, name(npest)
		svmat P, name(P)
		svmat beta, name(beta)
		label variable npest "NPIV fitted value"
		label variable P1 "NPIV bspline base 1"
		label variable P2 "NPIV bspline base 2"
		label variable P3 "NPIV bspline base 3"
		label variable beta "NPIV coefficients"
		scatter $depvar $expvar, msym(circle_hollow) || line npest $expvar, sort	

end
		
