clear
capture program drop npivregcv
capture program drop npivreg

// set the seed
set seed 1234

// number of observations = sample size
set obs 10000

// instrument z is generated from standard normal dist.
generate double z = rnormal(0, 1)
generate double w = rnormal(0, 1)

/* errors u and v are generated from Joint normal dist
    u       ~     N( 0,  1.0   0.5 )
	v              ( 0,  0.5   1.0 )
   correlation btw u and v = 0.5
*/
matrix C = (1, .5 \ .5, 1)
drawnorm u v, corr(C)

// DGP for x
generate x = 2*z + v

// DGP for y : g(x) is increasing in x
generate true_y = exp(0.5*x)/(1 + exp(0.5*x))
//generate true_y = sqrt(exp(x))
//generate true_y = exp(x/2)

generate y =  true_y + 0.5*w + u

// provide npivreg estimator 'npest' (fitted value) with increasing shape restriction
// and coefficients of series estimation
npivcv y x z w, pctile(1)
// save this one-step estimate in onestepest
generate onestepest = npest

// get Ey = W'γ
mata
bw = st_matrix("e(b)")
W  = st_data(., "w", 0)
n  = rows(bw)
np = n - cols(W)
Ey = bw[(np+1)..n]*W
end

getmata Ey, force

// define Y = y - Ey and run npivreg of Y on x only 
quietly generate Y = y - Ey
npivcv Y x z, pctile(5)

quietly line true_y x, sort || line npest grid || line onestepest grid
