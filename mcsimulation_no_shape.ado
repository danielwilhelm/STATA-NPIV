// This command produces a random sample from a data generating process
// and compute npiv estimate without shape restriction

program define mcsimulation_no_shape, eclass
capture program drop npiv
clear 
// sample size for mc simulation is only possible upto 800 in Stata IC.
// In MP or SE, it can be upto 11,000.
set obs 800
set matsize 800 

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
generate true_y = exp(0.5*x) / (1 + exp(0.5*x))
generate y =  true_y + u

// NPIV regression 
npiv y x z, power_exp(2) power_inst(3) num_exp(3) num_inst(4) pctile(1)
mkmat npest1, matrix(A)
matrix B = A'

// return the estimated function values at grid points
ereturn post B

end
