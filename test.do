clear
capture program drop npivreg
capture program drop npivoptim

// Load a data set "auto"
sysuse auto

// Compute NPIV estimate with bspline bases and draw a plot (using the closed form solution)
npivreg price mpg trunk, power_exp(2) power_inst(2) num_exp(2) num_inst(5) pctile(5) bspline
scatter price mpg, msym(circle_hollow) || line npest grid, title("power(2,2), knots(2,5)") name(bspline_result, replace)

// Compute NPIV estimate with polynomial bases and draw a plot (using the closed form solution)
npivreg price mpg trunk, power_exp(2) power_inst(2) num_exp(2) num_inst(5) pctile(5)
scatter price mpg, msym(circle_hollow) || line npest grid, title("power(2,2), knots(2,5)") name(poly_result, replace)

// Compute NPIV estimate with bspline bases and draw a plot (using optimization routine)
npivoptim price mpg trunk, power_exp(2) power_inst(2) num_exp(2) num_inst(5) pctile(5) bspline 
scatter price mpg, msym(circle_hollow) || line npest grid, title("power(2,2), knots(2,5)") name(bspline_result_op, replace)

// Compute NPIV estimate with polynomial bases and draw a plot (using optimization routine)
npivoptim price mpg trunk, power_exp(2) power_inst(2) num_exp(2) num_inst(5) pctile(5)
scatter price mpg, msym(circle_hollow) || line npest grid, title("power(2,2), knots(2,5)") name(poly_result_op, replace)

graph close _all

//draw a combined graph
gr combine poly_result bspline_result poly_result_op bspline_result_op, title("Polyspline(left) vs Bspline(right)")

