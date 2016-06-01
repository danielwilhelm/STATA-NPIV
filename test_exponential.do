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
generate x = 2*z + v
// DGP for y
// generate true_y = exp(1.5*x)
generate true_y = exp(0.5*x) / (1 + exp(0.5*x))
generate y =  true_y + u


// NPIV regression 
npivreg y x z, num_exp(2) num_inst(2) power_exp(2) power_inst(4)
// Comparison of true y and fitted value (drawing a chart)
quietly line true_y x, sort || line npest $expvar, sort title("power= (2,2), knots = (2,4)") name(setting_poly, replace)

npivreg y x z, num_exp(2) num_inst(2) power_exp(2) power_inst(4) bspline
quietly line true_y x, sort || line npest $expvar, sort title("power= (2,2), knots = (2,4)") name(setting_bspl, replace)

npivreg y x z, num_exp(2) num_inst(2) power_exp(4) power_inst(8)
quietly line true_y x, sort || line npest $expvar, sort title("power= (2,2), knots = (4,8)") name(setting_poly2, replace)

npivreg y x z, num_exp(5) num_inst(5) power_exp(4) power_inst(8) bspline
quietly line true_y x, sort || line npest $expvar, sort title("power= (5,5), knots = (4,8)") name(setting_bspl2, replace)

// npivreg y x z, num_exp(9) num_inst(10) power_exp(9) power_inst(10)
// quietly line true_y x, sort || line npest $expvar, sort title("power= (9,10), knots = (9,10)") name(setting_poly3, replace)

// npivreg y x z, num_exp(9) num_inst(10) power_exp(9) power_inst(10) bspline
// quietly line true_y x, sort || line npest $expvar, sort title("power= (9,10), knots = (9,10)") name(setting_bspl3, replace)

graph close _all

// Combined graph for all the settings
// gr combine setting_poly setting_bspl setting_poly2 setting_bspl2 setting_poly3 setting_bspl3, cols(2) title("Polyspline(left) vs bspline(right)")
gr combine setting_poly setting_bspl setting_poly2 setting_bspl2, cols(2) title("Polyspline(left) vs bspline(right)") ycom
