# STATA-NPIV (Author : Dongwoo Kim and Daniel Wilhelm)
This project provides a Stata command `npivreg` for nonparametric estimation of instrumental variable (NPIV) models with or without monotonicity restrictions and 'npivregcv' for cross validation.

Files contained in this package:

- The file 'npivreg.ado' contains the `npivreg` command.
- The file 'npivregcv.ado' contains the `npivregcv` command.
- The file 'example_no_shape_restriction.do' contains an example that simulates data and then estimates the NPIV model without the monotonicity restriction.
- The file 'example_shape_restriction.do' contains an example that simulates data and then estimates the NPIV model with the monotonicity restriction.
- The file 'mcsimulation_cv.do' contains an example that simulates data and then estimates the NPIV model with the monotonicity restriction and the optimal number of knots found by cross validation.


## Installation
1. Download the package.
2. Install the `bspline` and `polyspline` functions from ssc by opening STATA and typing
	
	```
	ssc install bspline
	ssc install polyspline
	```

3. Change into the directory containing this package.
4. Use the command `npivreg` as described below.

## Syntax
The command `npivreg` estimates the function `g(x)` in the NPIV model

```
Y = g(X) + e      E[e | Z] = 0
```

where
- Y is a scalar dependent variable (`depvar`) 
- X is a scalar endogenous variable (`expvar`)
- Z a scalar instrument (`inst`)

Syntax:

```
npivreg depvar expvar inst [if] [in] [, power_exp(#) power_inst(#) num_exp(#) num_inst(#) polynomial increasing decreasing]
```

where
- `power_exp` is the power of basis functions for X (defalut = 2)
- `power_inst` is the power of basis functions for Z (defalut = 3)
- `num_exp` is the number of knots for Z (defalut = 2) when using a B-spline basis
- `num_inst` is the number of knots for Z (defalut = 3) when using a B-spline basis
- `polonomial` option gives the basis functions for polynomial spline (default is B-spline).
- `increasing` option imposes an increasing shape restriction on function g(x).
- `decreasing` option imposes an decreasing shape restriction on function g(x).

In the current version of the code, shape restrictions can only be imposed with the B-spline basis. The code will produce an error when attempting to impose shape restrictions with a polynomial basis. When imposing a shape constraint, the current code also imposes the power of the B-spline for X to be 2.

Users can freely modify the power and the type of basis functions and the number of knots
when shape restrictions are not imposed.

If options are left unspecified, the command runs on the default settings.

The command `npivregcv` estimates the function `g(x)` in the NPIV model with cross validation. It generates the npiv estimate for `g(x)` and also shows the optimal number of knots when a user types 'display opt_knot'.

## Examples

NPIV estimation with default options:
```
npivreg y x z
```

NPIV estimation with B-spline bases of powers 2 and 3, and 3 and 4 knots (for X and Z, respectively):
```
npivreg y x z, power_exp(2) power_inst(3) num_exp(3) num_inst(4)
```

NPIV estimation with B-spline bases of powers 2 and 3, and 3 and 4 knots (for X and Z, respectively), imposing that the estimator of g(x) is increasing:
```
npivreg y x z, power_exp(2) power_inst(3) num_exp(3) num_inst(4) increasing
```
NPIV estimation with cross validation, imposing that the estimator of g(x) is increasing
```
npivregcv y x z, power_exp(2) power_inst(3) increasing
```


# Reference
[Chetverikov, D. and Wilhelm, D. (2016), "Nonparametric Instrumental Variable Estimation Under Monotonicity", cemmap working paper](http://www.ucl.ac.uk/~uctpdwi/papers/cwp481616.pdf)
