#' @title kssa_plot function
#'
#' @description Function to plot the results of kssa for easy interpretation
#' @param results An object with results produced with function \code{\link{kssa}}
#' @param type A character value with the type of plot to show.
#' It can be "summary" or "complete".
#'
#' @param metric A character with the performance metric to be plotted.
#' It can be  "rmse", "mase," "cor", or "smape"
#'
#' \itemize{
#'    \item{"rmse" -  Root Mean Squared Error} (default choice)
#'    \item{"mase" - Mean Absolute Scaled Error}
#'    \item{"smape" - Symmetric Mean Absolute Percentage Error}
#'    \item{"cor" - Pearson correlation coefficient}
#'    }
#'
#' For further details on these metrics please check package Metrics
#'
#' @return A plot of kssa results in which imputation methods are ordered
#' from lower to higher (left to right) error.
#'
#' @examples
#' \donttest{
#' # Plot the results obtained in the example from function kssa
#'
#' # Create 20% random missing data in tsAirgapComplete time series from imputeTS
#' set.seed(1234)
#' library("kssa")
#' library("imputeTS")
#' airgap_na <- missMethods::delete_MCAR(as.data.frame(tsAirgapComplete), 0.2)
#'
#' # Convert co2_na to time series object
#' airgap_na_ts <- ts(airgap_na, start = c(1959, 1), end = c(1997, 12), frequency = 12)
#'
#' # Apply the kssa algorithm with 5 segments,
#' # 10 iterations, 20% of missing data, and
#' # compare among all available methods in the package.
#' # Remember that percentmd must match with
#' # the real percentage of missing data in the
#' # input co2_na_ts time series
#'
#' results_kssa <- kssa(airgap_na_ts,
#'   start_methods = "all",
#'   actual_methods = "all",
#'   segments = 5,
#'   iterations = 10,
#'   percentmd = 0.2
#' )
#'
#' kssa_plot(results_kssa, type = "complete", metric = "rmse")
#'
#' # Conclusion: Since kssa_plot is ordered from lower to
#' # higher error (left to right), method 'linear_i' is the best to
#' # impute missing data in airgap_na_ts. Notice that method 'locf' is the worst
#'
#' # To obtain imputations with the best method, or any method of preference
#' # please use function get_imputations
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom stats reorder
#' @importFrom ggplot2 .data
#' @importFrom methods is
#' @export

kssa_plot <- function(results, # results object from apply kssa function
                      type, # complete or summary
                      metric # if type == complete put any metric like rmse, cor, mase, smape
) {
  # First condition with kssa.object
  if (is(results[[1]],"kssa.table") & is(results[[2]],"kssa.table")) {

    # Extract complete data from kssa.tables object
    plot1 <- data.frame(
      "start_methods" = results[[1]]$start_methods,
      "actual_methods" = results[[1]]$actual_methods,
      "rmse" = results[[1]]$rmse,
      "cor" = results[[1]]$cor,
      "mase" = results[[1]]$mase,
      "smape" = results[[1]]$smape
    )

    # Extract summary information from kssa.table objects
    plot2 <- data.frame(
      "start_methods" = results[[2]]$start_methods,
      "actual_methods" = results[[2]]$actual_methods,
      "mean_rmse" = results[[2]]$mean_rmse,
      "std_rmse" = results[[2]]$std_rmse,
      "mean_cor" = results[[2]]$mean_cor,
      "std_cor" = results[[2]]$std_cor,
      "mean_mase" = results[[2]]$mean_mase,
      "std_mase" = results[[2]]$std_mase,
      "mean_smape" = results[[2]]$mean_smape,
      "std_smape" = results[[2]]$std_smape
    )

    # Generates list of metrics
    list_of_metrics <- c("rmse", "cor", "mase", "smape")

    # Get check for validation process
    check <- metric %in% list_of_metrics


    # Process validation and if all clear, ---> Plot
    if (type == "complete") {
      if (all(check) == TRUE) {
        plot <- ggplot(data = plot1, aes(x = stats::reorder(.data$actual_methods, eval(parse(text = metric))), y = eval(parse(text = metric)), fill = .data$start_methods)) +
          geom_boxplot() +
          labs(title = "Complete plot", x = "Actual methods", y = casefold(metric), fill = "Start methods") +
          theme_classic()
      }
      else {
        print('metric must be "rmse","cor", "mase" or "smape".
             Otherwise please check the defined metrics when run kssa')
      }
    }
    else if (type == "summary") {
      if (all(check) == TRUE) {
        plot <- ggplot(data = plot2, aes(x = stats::reorder(.data$actual_methods, eval(parse(text = paste0("mean_", metric)))), y = eval(parse(text = paste0("mean_", metric))), fill = .data$start_methods)) +
          geom_bar(
            stat = "identity", color = "black",
            position = position_dodge()
          ) +
          geom_errorbar(aes(
            ymin = eval(parse(text = paste0("mean_", metric))) - eval(parse(text = paste0("std_", metric))),
            ymax = eval(parse(text = paste0("mean_", metric))) + eval(parse(text = paste0("std_", metric)))
          ),
          width = .2,
          position = position_dodge(.9)
          ) +
          labs(title = "Summary plot", x = "Actual methods", y = casefold(metric)) +
          theme_classic()
      }
      else {
        print('metric must be "rmse","cor", "mase" or "smape".
             Otherwise please check the defined metrics when run kssa')
      }
    }
    else {
      print('type must be "complete" or "summary"')
    }
  }
  else {
    print("Error: data must be a kssa.table object")
  }
  return(plot)
} # End of code
