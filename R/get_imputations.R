#' @title get_imputations function
#'
#' @description Function to get imputation
#' @param x_ts A ts object to be imputed
#' @param method A character with method or methods types for imputation
#'
#' @return The imputated \code{x_ts}
#'
#' @importFrom imputeTS na_kalman na_interpolation na_ma na_seadec na_locf
#' @importFrom forecast na.interp
#' @importFrom stats sd cor ts window
#' @importFrom Metrics rmse mase mape smape
#' @importFrom zoo coredata index
#' @importFrom rlang .data
#' @export

get_imputations <- function(
    x_ts, # Time-series
    methods = 'all', # Can select various or all
    seed = 1234
){
  # Generate df of imputation methods and formulas
  df_of_methods <- data.frame(
    "methods" =c("auto.arima", "StructTS", "linear",
                 "spline", "stine", "simple", "malinear",
                 "exponential", "kalman", "nalocf", "decomp"),
    "formulas_x_ts" = c("na_kalman(x_ts,model='auto.arima',smooth = TRUE,nit = -1)",
                        "na_kalman(x_ts,model='StructTS',smooth = TRUE,nit = -1)",
                        "na_interpolation(x_ts,option='linear')",
                        "na_interpolation(x_ts,option='spline')",
                        "na_interpolation(x_ts,option='stine')",
                        "na_ma(x_ts,k=3,weighting = 'simple')",
                        "na_ma(x_ts,k=3,weighting = 'linear')",
                        "na_ma(x_ts,k=3,weighting = 'exponential')",
                        "na_seadec(x_ts,algorithm = 'kalman')",
                        "na_locf(x_ts,option = 'locf',na_remaining = 'rev')",
                        "na.interp(x_ts)"),
    "formulas_actual_ts" = c("na_kalman(newmdsimulation,model='auto.arima',smooth = TRUE,nit = -1)",
                             "na_kalman(newmdsimulation,model='StructTS',smooth = TRUE,nit = -1)",
                             "na_interpolation(newmdsimulation,option='linear')",
                             "na_interpolation(newmdsimulation,option='spline')",
                             "na_interpolation(newmdsimulation,option='stine')",
                             "na_ma(newmdsimulation,k=3,weighting = 'simple')",
                             "na_ma(newmdsimulation,k=3,weighting = 'linear')",
                             "na_ma(newmdsimulation,k=3,weighting = 'exponential')",
                             "na_seadec(newmdsimulation,algorithm = 'kalman')",
                             "na_locf(newmdsimulation,option = 'locf',na_remaining = 'rev')",
                             "na.interp(newmdsimulation)")
  )

  if(length(methods) == 1 && methods == 'all'){ # Define 'all' statement
    methods <- c("auto.arima", "StructTS", "linear",
                 "spline", "stine", "simple", "malinear",
                 "exponential", "kalman", "nalocf", "decomp")
  } else {
    methods <- methods
  }

  check <- methods %in% df_of_methods$methods

  if (all(check) == TRUE){
    results = c()
    for (i in 1:length(methods)) {
      # First imputation
      set.seed(seed); imputation <- eval(parse(text = df_of_methods$formulas_x_ts[df_of_methods$methods == methods[i]]))
      results <- append(results, list(imputation))
    }
    names(results) <- methods
    return(results)
  } else {
    print(paste0("The methods '",
                 paste(as.character(methods[which(!check)]),
                       collapse = ", "),
                 "' in actual_method parameter, are not in the list of available options"))
  }
}
