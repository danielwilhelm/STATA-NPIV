clear
capture program drop npivreg

// number of observations = sample size
set obs 10000

// instrument z is generated from standard normal dist.
generate double z = rnormal(0, 1)

/* errors u and v are generated from Joint normal dist
    u       ~     N( 0,  1.0   0.5 )
	v              ( 0,  0.5   1.0 )
   correlation btw u and v = 0.5
*/
matrix C = (1, .5 \ .5, 1)
drawnorm u v, n(10000) corr(C)

// DGP for x
generate x = 2*z + v
// DGP for y
generate y = exp(0.5*x) + u
generate true_y = exp(0.5*x)

// NPIV regression
npivreg y x z, num_exp(3) num_inst(4) power_exp(3) power_inst(4) bspline

// Comparison of true y and fitted value
line true_y x, sort || line npest $expvar, sort
