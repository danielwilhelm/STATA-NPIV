clear
program drop npivreg

// Load a data set "auto"
sysuse auto

// Compute NPIV estimate with bspline bases and draw a plot 
npivreg price mpg trunk, power_exp(4) power_inst(5) num_exp(4) num_inst(5) bspline
scatter price mpg, msym(circle_hollow) || line npest mpg, sort title("power(4,5), knots(4,5)") name(bspline_result, replace)

// Compute NPIV estimate with polynomial bases and draw a plot 
npivreg price mpg trunk, power_exp(4) power_inst(5) num_exp(4) num_inst(5) 
scatter price mpg, msym(circle_hollow) || line npest mpg, sort title("power(4,5), knots(4,5)") name(poly_result, replace)

graph close _all
// draw a combined graph
gr combine poly_result bspline_result, title("Polyspline(left) vs Bspline(right)")

