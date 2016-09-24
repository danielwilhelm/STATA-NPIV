clear
capture program drop npivreg

// set the seed
set seed 1234

// number of observations = sample size
set obs 1000

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
generate true_y = exp(x)
generate y =  true_y + u

// NPIV regression with decreasing shape restriction
npivreg y x z, power_exp(2) power_inst(3) num_exp(4) num_inst(4) pctile(1) increasing
generate ygrid = exp(grid)
quietly line ygrid grid, sort || line npest grid, title("increasing restriction") name(increasing, replace) 

capture drop ygrid true_y y
// DGP for y : g(x) is decreasing
generate true_y = -exp(x)
generate y =  true_y + u

// NPIV regression with decreasing shape restriction
npivreg y x z, power_exp(2) power_inst(3) num_exp(4) num_inst(4) pctile(1) decreasing
generate ygrid = -exp(grid)
quietly line ygrid grid, sort || line npest grid, title("decreasing restriction") name(decreasing, replace) 

graph close _all

// Combined graph 
gr combine increasing decreasing, cols(2) title("increasing(left) and decreasing(right)") ycom
