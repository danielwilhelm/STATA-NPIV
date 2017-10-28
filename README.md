# STATA-NPIV 
Authors : Dongwoo Kim and Daniel Wilhelm

This project provides two Stata commands for nonparametric estimation of instrumental variable (NPIV) models with or without imposing monotonicity restrictions on the function of interest. The command `npivreg` implements the estimators with user-chosen tuning parameters and `npivregcv` with tuning parameters chosen by cross-validation.

Files contained in this package:

- The file `npivreg.ado` contains the `npivreg` command.
- The file `npivregcv.ado` contains the `npivregcv` command.
- The file `example_no_shape_restriction.do` contains an example that simulates data and then estimates the NPIV model without imposing the monotonicity restriction.
- The file `example_shape_restriction.do` contains an example that simulates data and then estimates the NPIV model imposing the monotonicity restriction.
- The file `mcsimulation_cv.do` contains an example that simulates data and then estimates the NPIV model imposing the monotonicity restriction and using cross-validated tuning parameters.


## Installation
1. Download the package.
2. Install the `bspline` and `polyspline` functions from ssc by opening STATA and typing
	
	```
	ssc install bspline
	ssc install polyspline
	```

3. Change into the directory containing this package.
4. Use the commands `npivreg` and `npivregcv` as described below.

## Syntax
The commands `npivreg` and `npivregcv` estimate the function g(x) in the NPIV model

```
Y = g(X) + Z'γ+ e      E[e | Z, W] = 0
```

where
- Y is a scalar dependent variable (`depvar`) 
- X is a scalar endogenous variable (`expvar`)
- W is a scalar instrument (`inst`)
- Z is a vector of exogeneous covariats (`exovar`) - optional (need not be included), and 

Syntax:

```
npivreg depvar expvar inst [exovar] [if] [in] [, power_exp(#) power_inst(#) num_exp(#) num_inst(#) polynomial increasing decreasing]
npivregcv depvar expvar inst [exovar] [if] [in] [, power_exp(#) power_inst(#) maxknot(#) polynomial increasing decreasing]
```

where
- `power_exp` is the power of basis functions for X (default = 2).
- `power_inst` is the power of basis functions for W (default = 3).
- `num_exp` is the number of knots for X (default = 2) when using a B-spline basis.
- `num_inst` is the number of knots for W (default = 3) when using a B-spline basis.
- `maxknot` is the maximum number of knots upto which cross-validation is executed. (default = 5)
- adding the `polynomial` option makes the estimator use a polynomial spline basis (default is B-spline basis).
- adding the `increasing` option makes the estimator impose the constraint that the estimator of g(x) is increasing in x.
- adding the `decreasing` option makes the estimator impose the constraint that the estimator of g(x) is decreasing in x.

In the current version of the code, shape restrictions can only be imposed with the B-spline basis. The code will produce an error when attempting to impose shape restrictions with a polynomial basis. When imposing a shape constraint, the current code also imposes the power of the B-spline for X to be 2.

Users can freely modify the power and the type of basis functions and the number of knots
when shape restrictions are not imposed.

If options are left unspecified, the command runs on the default settings.

The command `npivregcv` estimates the function g(x) in the NPIV model using cross-validation to find the optimal number of knots (setting `num_exp=num_inst`).


## Output

The commands `npivreg` and `npivregcv` save their estimates of g(x) over a grid of values for x in the variable `npest`. In addition, the commands are of the e-class so other results such as the coefficient vector are stored in e(). 


## Examples

NPIV estimation with default options:
```
npivreg y x w z
```

NPIV estimation with B-spline bases of powers 2 and 3, and 3 and 4 knots (for X and W, respectively):
```
npivreg y x w z, power_exp(2) power_inst(3) num_exp(3) num_inst(4)
```

NPIV estimation with B-spline bases of powers 2 and 3, and 3 and 4 knots (for X and W, respectively), imposing that the estimator of g(x) is increasing:
```
npivreg y x w z, power_exp(2) power_inst(3) num_exp(3) num_inst(4) increasing
```
NPIV estimation using cross-validation to determine the optimal number of knots, imposing that the estimator of g(x) is increasing:
```
npivregcv y x w z, power_exp(2) power_inst(3) increasing
```


# Reference
[Chetverikov, D. and Wilhelm, D. (2017), "Nonparametric Instrumental Variable Estimation Under Monotonicity", Econometrica, 85(4), p. 1303-1320](http://onlinelibrary.wiley.com/doi/10.3982/ECTA13639/full)

[Chetverikov, D. and Wilhelm, D. (2017), "Nonparametric Instrumental Variable Estimation Under Monotonicity", cemmap working paper](http://www.ucl.ac.uk/~uctpdwi/papers/cwp141717.pdf)
