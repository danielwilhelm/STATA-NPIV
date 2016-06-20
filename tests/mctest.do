// This do file provides Monte Carlo simulation results with mctest_npiv.ado file.
// Modify the number of reps and run

clear
program drop mctest_npiv

// This stata implemented command does MC simulation
// The estimated function in each trial is stored in _b_r# vector
simulate _b, rep(5000) : mctest_npiv

// Mata codes for simple matrix algebra
mata
est      = st_data(.,.) // load simulation results
meanest  = mean(est)'   // compute average estimate
std      = sqrt(diagonal(variance(est))) // compute standard deviation at each grid point
ci_left  = meanest - std*1.96 // 2.5 percentile of reps  
ci_right = meanest + std*1.96 // 97.5 percentile of reps
end

// create grid and true value
mctest_npiv 
generate true = exp(0.5*grid) / (1 + exp(0.5*grid))

// Mata matrices to Stata matrices
mata
st_matrix("meanest", meanest)
st_matrix("ci_left", ci_left)
st_matrix("ci_right", ci_right)
end

// Store Stata matrices in variable space
svmat meanest, name(meanest)
svmat ci_left, name(ci_left)
svmat ci_right, name(ci_right)

// Line graph
#delimit ;
line true grid, sort || line meanest grid || line ci_left grid, lcolor(red) || line ci_right grid, 
    lcolor(red)	xtitle("x") ytitle("y")	title("Monte Carlo Simulation") 
	legend(label(1 "True function") label(2 "Averaged estimate") label(3 "-2SD") label(4 "+2SD"));
#delimit cr
	
