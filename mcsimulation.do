// This do file provides Monte Carlo simulation results with mctest_npiv.ado file.
// Modify the number of reps and run

clear
capture program drop mcsimulation_no_shape
capture program drop mcsimulation_shape

// This stata implemented command does MC simulation
// simulation for npiv
// # of knots, powers can be modified in mctest_npiv
// The estimated function in each trial is stored in _b_r# vector
simulate _b, rep(10) : mcsimulation_no_shape

// Mata codes for simple matrix algebra
mata
est      = st_data(.,.) // load simulation results
meanest  = mean(est)'   // compute average estimate
std      = sqrt(diagonal(variance(est))) // compute standard deviation at each grid point
ci_left  = meanest - std*1.96 // 2.5 percentile of reps  
ci_right = meanest + std*1.96 // 97.5 percentile of reps
end

// create grid and true value
mcsimulation_no_shape
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
    lcolor(red)	xtitle("x") ytitle("y")	title("npiv") 
	legend(label(1 "True fn") label(2 "Aver. est.") label(3 "-2SD") label(4 "+2SD")) name(monotone, replace);
#delimit cr

// simulation for npivmonotne
// # of knots, powers can be modified in mctest_npivmono
// The estimated function in each trial is stored in _b_r# vector
simulate _b, rep(10) : mcsimulation_shape

// Mata codes for simple matrix algebra
mata
est      = st_data(.,.) // load simulation results
meanest  = mean(est)'   // compute average estimate
std      = sqrt(diagonal(variance(est))) // compute standard deviation at each grid point
ci_left  = meanest - std*1.96 // 2.5 percentile of reps  
ci_right = meanest + std*1.96 // 97.5 percentile of reps
end

// create grid and true value
mcsimulation_no_shape 
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
    lcolor(red)	xtitle("x") ytitle("y")	title("npivmonotone") 
	legend(label(1 "True fn") label(2 "Aver. est.") label(3 "-2SD") label(4 "+2SD")) name(unconst, replace);
#delimit cr

graph close _all

gr combine monotone unconst, title("npiv(left) vs npivmonotone(right)") ycommon xcommon
