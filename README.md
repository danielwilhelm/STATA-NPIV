# STATA-NPIV
This project provides a Stata command `npivreg` for nonparametric estimation of instrumental variable models with or without monotonicity restrictions.

Files contained in this package:

- The file 'npivreg.ado' contains the `npivreg` command.
- The file 'example_no_shape_restriction.do' contains an example that simulates data and then estimates the NPIV model without the monotonicity restriction.
- The file 'example_shape_restriction.do' contains an example that simulates data and then estimates the NPIV model with the monotonicity restriction.

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