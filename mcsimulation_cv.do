clear
capture program drop npivregcv
capture program drop npivreg

// set the seed
set seed 1234

// number of observations = sample size
set obs 100

// instrument z is generated from standard normal dist.
generate double z = rnormal(0, 1)

/* errors u and v are generated from Joint normal dist
    u       ~     N( 0,  1.0   0.5 )
	v              ( 0,  0.5   1.0 )
   correlation btw u and v = 0.5
*/
matrix C = (1, .5 \ .5, 1)
drawnorm u v, corr(C)

// DGP for x
generate x = z + v

// DGP for y : g(x) is increasing in x
//generate true_y = exp(x/2)/(1+exp(x/2))
generate true_y = -exp(x/2)
generate y =  true_y + u

// provide npivreg estimator 'npest' (fitted value)
// and coefficients of series estimation
npivregcv y x z, pctile(5) decreasing

// show optimal knots by cross validation
display opt_knot

quietly line true_y x, sort || line npest grid
