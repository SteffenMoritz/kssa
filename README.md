
# kssa

<!-- badges: start -->
[![R-CMD-check](https://github.com/SteffenMoritz/kssa/workflows/R-CMD-check/badge.svg)](https://github.com/SteffenMoritz/kssa/actions)
[![Codecov test coverage](https://codecov.io/gh/SteffenMoritz/kssa/branch/master/graph/badge.svg)](https://app.codecov.io/gh/SteffenMoritz/kssa?branch=master)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The goal of kssa is to ...

## Installation

You can install the development version of kssa like so:

``` r
library(devtools)
install_github("SteffenMoritz/kssa")
```

## Example

# create a numeric vector with 20% missing data
#' x = c(1, 5, 6, 8, 4, NA, 5, 4, NA, NA)
#'
#' # convert x to a time series object
#' x_ts = ts(x)
#'
#' # apply the kssa algorithm with 2 segments,
#' # 10 iterations and 20% of missing data.
#' # Remember that percentmd must match with
#' # the real percentaje of missing data in the
#' # input time series
#'
#' results_kssa = kssa(x_ts,
#'                start_method = "all",
#'                methods = "all",
#'                segments = 2,
#'                iterations = 10,
#'                percentmd = 0.2)
#'
#' # print results
#' results_kssa
#'
#' # plot complete results with Root Mean Squared Error for easy
#' # interpretation
#' kssa_plot(results_kssa, type = "complete", metric = "rmse")
#'
#' # Conclusion: Since the kssa_plot is ordered from lower to
#' # higher (left to right) average error, the method
#' # exponential_ma (exponential moving average) is
#' # the best to impute missing data in x_tx.
#'


``` r
library(kssa)
## basic example code
```
git config --global user.email "you@example.com"
