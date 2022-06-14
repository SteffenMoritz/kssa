#' @title get_imputations function
#'
#' @description Function to get imputations from methods compared by kssa
#' @param x_ts A ts object with missing data to be imputed
#' @param methods A string or string vector indicating the method or methods
#' @param seed Numeric. Any number
#'
#' @return A list of imputed time series with the selected methods
#'
#' @examples
#' \donttest{
#' # Get imputed values for airgap_na_ts with the methods of
#' # Create 20% random missing data in tsAirgapComplete time series from imputeTS
#' set.seed(1234)
#' library("imputeTS")
#' library("kssa")
#' airgap_na <- missMethods::delete_MCAR(as.data.frame(tsAirgapComplete), 0.2)
#'
#' # Convert co2_na to time series object
#' airgap_na_ts <- ts(airgap_na, start = c(1959, 1), end = c(1997, 12), frequency = 12)
#'
#' my_imputations <- get_imputations(airgap_na_ts, methods = "all")
#'
#' # my_imputations contains the imputed time series with all methods.
#' # Access it and choose the one from the best method for your purposes
#'
#' my_imputations$seadec
#' plot.ts(my_imputations$seadec)
#' }
#'
#' @importFrom imputeTS na_kalman na_interpolation na_ma na_seadec na_locf
#' @importFrom forecast na.interp
#' @importFrom stats sd cor ts window
#' @importFrom Metrics rmse mase mape smape
#' @importFrom zoo coredata index
#' @importFrom rlang .data
#' @export

get_imputations <- function(x_ts, # Time-series
                            methods = "all", # Can select various or all
                            seed = 1234) {
  # Generate df of imputation methods and formulas
  df_of_methods <- data.frame(
    "methods" = c(
      "auto.arima", "StructTS", "linear_i",
      "spline_i", "stine_i", "simple_ma", "linear_ma",
      "exponential_ma", "seadec", "locf", "stl"
    ),
    "formulas_x_ts" = c(
      "na_kalman(x_ts,model='auto.arima',smooth = TRUE,nit = -1)",
      "na_kalman(x_ts,model='StructTS',smooth = TRUE,nit = -1)",
      "na_interpolation(x_ts,option='linear')",
      "na_interpolation(x_ts,option='spline')",
      "na_interpolation(x_ts,option='stine')",
      "na_ma(x_ts,k=3,weighting = 'simple')",
      "na_ma(x_ts,k=3,weighting = 'linear')",
      "na_ma(x_ts,k=3,weighting = 'exponential')",
      "na_seadec(x_ts,algorithm = 'kalman')",
      "na_locf(x_ts,option = 'locf',na_remaining = 'rev')",
      "na.interp(x_ts)"
    ),
    "formulas_actual_ts" = c(
      "na_kalman(newmdsimulation,model='auto.arima',smooth = TRUE,nit = -1)",
      "na_kalman(newmdsimulation,model='StructTS',smooth = TRUE,nit = -1)",
      "na_interpolation(newmdsimulation,option='linear')",
      "na_interpolation(newmdsimulation,option='spline')",
      "na_interpolation(newmdsimulation,option='stine')",
      "na_ma(newmdsimulation,k=3,weighting = 'simple')",
      "na_ma(newmdsimulation,k=3,weighting = 'linear')",
      "na_ma(newmdsimulation,k=3,weighting = 'exponential')",
      "na_seadec(newmdsimulation,algorithm = 'kalman')",
      "na_locf(newmdsimulation,option = 'locf',na_remaining = 'rev')",
      "na.interp(newmdsimulation)"
    )
  )

  if (length(methods) == 1 && methods == "all") { # Define 'all' statement
    methods <- c(
      "auto.arima", "StructTS", "linear_i",
      "spline_i", "stine_i", "simple_ma", "linear_ma",
      "exponential_ma", "seadec", "locf", "stl"
    )
  } else {
    methods <- methods
  }

  check <- methods %in% df_of_methods$methods

  if (all(check) == TRUE) {
    results <- c()
    for (i in 1:length(methods)) {
      # First imputation
      set.seed(seed)
      imputation <- eval(parse(text = df_of_methods$formulas_x_ts[df_of_methods$methods == methods[i]]))
      results <- append(results, list(imputation))
    }
    names(results) <- methods
    return(results)
  } else {
    print(paste0(
      "The methods '",
      paste(as.character(methods[which(!check)]),
        collapse = ", "
      ),
      "' in actual_method parameter, are not in the list of available options"
    ))
  }
}
