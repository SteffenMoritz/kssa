#' @title kssa_plot function
#'
#' @description Function to plot kssa results
#' @param results A list object with kssa.table objects
#' @param type A character with type of plot complete or summary
#' @param metric A character with metric choosed to be plot
#'
#' @return The plot of \code{results}
#'
#' @import ggplot2
#' @import kssa
#' @import dplyr
#' @import zoo
#' @export

kssa_plot <- function(
  results, #' results object from apply kssa function
  type, #' complete or summary
  metric #' if type == complete put any metric like rmse, cor, mase, smape
){
  #' Fisrt condition with kssa.object
  if (class(results[[1]]) == 'kssa.table' & class(results[[2]]) == 'kssa.table'){

    #' Extract complete data from kssa.tables object
    plot1 <- data.frame("start.method" = results[[1]]$start.method,
                              "actual.method" = results[[1]]$actual.method,
                              "rmse" = results[[1]]$rmse,
                              "cor" = results[[1]]$cor,
                              "mase" = results[[1]]$mase,
                              "smape"= results[[1]]$smape)

    #' Extract summary information from kssa.table object
    plot2 <- data.frame("start.method" = results[[2]]$start.method,
                        "actual.method" = results[[2]]$actual.method,
                        "mean_rmse" = results[[2]]$mean_rmse,
                        "std_rmse" = results[[2]]$std_rmse,
                        "mean_cor" = results[[2]]$mean_cor,
                        "std_cor" = results[[2]]$std_cor,
                        "mean_mase" = results[[2]]$mean_mase,
                        "std_mase" = results[[2]]$std_mase,
                        "mean_smape"= results[[2]]$mean_smape,
                        "std_smape"= results[[2]]$std_smape)

    #' Generates list of metrics
    list_of_metrics <- c("rmse", "cor", "mase", "smape")

    #' Get check for validation process
    check <- metric %in% list_of_metrics


    #' Process validation and if all clear, ---> Plot
    if (type == "complete"){
      if (all(check) == TRUE){
        plot <- ggplot(data = plot1, aes(x = reorder(actual.method, eval(parse(text = metric))), y = eval(parse(text = metric)), fill = start.method)) +
          geom_boxplot()+
          labs(title = "Complete plot", x = "Actual method", y = casefold(metric), fill = "Start method") +
          theme_classic()
      }
      else {
        prin('metric must be "rmse","cor", "mase" or "smape".
             Otherwise please check the defined metrics when run kssa')
      }
    }
    else if (type == "summary"){
      if (all(check) == TRUE){
        plot <- ggplot(data = plot2, aes(x = reorder(actual.method, eval(parse(text = paste0('mean_',metric)))), y = eval(parse(text = paste0('mean_',metric))), fill = start.method)) +
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
        prin('metric must be "rmse","cor", "mase" or "smape".
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
} #' End
