clear
capture program drop npivreg

// set the seed
set seed 1234

// number of observations = sample size
set obs 100

// instrument z is generated from standard normal dist.
generate double z = rnormal(0, 1)
generate double w1 = rnormal(0, 1)
generate double w2 = rnormal(0, 1)

/* errors u and v are generated from Joint normal dist
    u       ~     N( 0,  1.0   0.5 )
	v              ( 0,  0.5   1.0 )
   correlation btw u and v = 0.5
*/
matrix C = (1, .5 \ .5, 1)
drawnorm u v, corr(C)

// DGP for x
generate x = 2*z + v
// DGP for y
// generate true_y = exp(1.5*x)
generate true_y = exp(0.5*x) / (1 + exp(0.5*x))
generate y =  true_y + 0.5*w1 + 0.3*w2 + u
//generate y =  true_y + 0.5*w1 + u


// NPIV regression 
npivreg y x z w1 w2, power_exp(2) power_inst(3) num_exp(3) num_inst(4) pctile(2) increasing
// Comparison of true y and fitted value (drawing a chart)
quietly line true_y x, sort || line npest grid
