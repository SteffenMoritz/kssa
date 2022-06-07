#' @title kssa_plot function
#'
#' @description Function to plot the results of kssa for easy interpretation
#' @param results A list object with results produced with function \code{\link{kssa}}
#' @param type A character value with the type of plot to show.
#' It can be "summary" or "complete".
#'
#' @param metric A character with the performance metric to be plotted
#' It can be  "rmse", "mase," "cor", or "smape"
#'
#' For further details on these metrics please check package \code{\link{Metrics}}
#'
#' @return A plot of kssa results ordered from lower to higher (left to right)
#' averaged error
#'
#' @examples # create a numeric vector with 20% missing data
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
#' @import ggplot2
#' @import dplyr
#' @importFrom stats reorder
#' @importFrom ggplot2 .data
#' @export

kssa_plot <- function(
  results, # results object from apply kssa function
  type, # complete or summary
  metric # if type == complete put any metric like rmse, cor, mase, smape
){
  # First condition with kssa.object
  if (class(results[[1]]) == 'kssa.table' & class(results[[2]]) == 'kssa.table'){

    # Extract complete data from kssa.tables object
    plot1 <- data.frame("start_method" = results[[1]]$start_method,
                              "actual_method" = results[[1]]$actual_method,
                              "rmse" = results[[1]]$rmse,
                              "cor" = results[[1]]$cor,
                              "mase" = results[[1]]$mase,
                              "smape"= results[[1]]$smape)

    # Extract summary information from kssa.table object
    plot2 <- data.frame("start_method" = results[[2]]$start_method,
                        "actual_method" = results[[2]]$actual_method,
                        "mean_rmse" = results[[2]]$mean_rmse,
                        "std_rmse" = results[[2]]$std_rmse,
                        "mean_cor" = results[[2]]$mean_cor,
                        "std_cor" = results[[2]]$std_cor,
                        "mean_mase" = results[[2]]$mean_mase,
                        "std_mase" = results[[2]]$std_mase,
                        "mean_smape"= results[[2]]$mean_smape,
                        "std_smape"= results[[2]]$std_smape)

    # Generates list of metrics
    list_of_metrics <- c("rmse", "cor", "mase", "smape")

    # Get check for validation process
    check <- metric %in% list_of_metrics


    # Process validation and if all clear, ---> Plot
    if (type == "complete"){
      if (all(check) == TRUE){
        plot <- ggplot(data = plot1, aes(x = stats::reorder(.data$actual_method, eval(parse(text = metric))), y = eval(parse(text = metric)), fill = .data$start_method)) +
          geom_boxplot()+
          labs(title = "Complete plot", x = "Actual method", y = casefold(metric), fill = "Start method") +
          theme_classic()
      }
      else {
        print('metric must be "rmse","cor", "mase" or "smape".
             Otherwise please check the defined metrics when run kssa')
      }
    }
    else if (type == "summary"){
      if (all(check) == TRUE){
        plot <- ggplot(data = plot2, aes(x = stats::reorder(.data$actual_method, eval(parse(text = paste0('mean_',metric)))), y = eval(parse(text = paste0('mean_',metric))), fill = .data$start_method)) +
          geom_bar(stat="identity", color="black",
                   position=position_dodge()) +
          geom_errorbar(aes(ymin=eval(parse(text = paste0('mean_',metric)))-eval(parse(text = paste0('std_',metric))),
                            ymax=eval(parse(text = paste0('mean_',metric)))+eval(parse(text = paste0('std_',metric)))),
                        width=.2,
                        position=position_dodge(.9)) +
          labs(title = "Summary plot", x = "Actual method", y = casefold(metric)) +
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
