#' @title kssa Algorithm
#'
#' @description Function to add two numbers
#' @param x A number
#' @param y A number
#'
#' @return The sum of \code{x} and \code{y}
#'
#' @import imputeTS
#' @import imputeFin
#' @import sjmisc
#' @import Metrics
#' @import dplyr
#' @import snpar
#' @import zoo
#' @import forecast
#' @export

kssa <- function(x_ts, #Time-series
                 start.method, # Can select various
                 methods, # Can select various
                 segments, #Number of segments to ts be divided
                 iterations, #Number of replicates
                 percentmd) { #Percentage of new MD in simulations

  results <- data.frame( #Create data frame where put the final results
    "start.method" = as.character(),
    "actual.method" = as.character(),
    "rmse" = numeric(),
    "cor" = numeric(),
    "mase" = numeric(),
    "smape"= numeric()
    )

#Function to get original positions of MD
  get_mdoriginal <- function(x, y){
    na_original <- which(is.na(x)) #Allows get NA index
    weight_index <- rep(1, length(y)) #Asign weigth of 1
    weight_index[na_original] <- 0 #Asign weigth of 0 to original MD
    return(weight_index) #Return index
  }

#Function to split TS
    split_arbitrary <- function(x, percentmd, segments, mdoriginal){ #Arguments
    size_window_B <- seq(#generate B point of time window
      from = round(length(x)/segments), #from length / nparts
                         to = length(x), #to max length of TS
                         by = round(length(x)/segments)+1 #step by step
      )
    size_window_A <- size_window_B + 1 #Generate A point of time window
    size_window_B <- c(size_window_B, length(x)) #Put last one length
    #size_window_A <- c(size_window_A, size_window_B[length(size_window_B)-1]+1)
    size_window_A <- c(1, size_window_A) #Put firs point on side A

    index_time <- index(x) #get indexes of window

    chunks <- c() #Vector to contain segments

    for (i in 1:segments) { #Run loop over segments
      chunk <- window(x = x, start = index_time[size_window_A[i]], #create chunk
                      end = index_time[size_window_B[i]])

      m_a2 <- sample(x = index(chunk), #Take new sample for simlate new MD
                     size = (round(length(chunk)*percentmd)), replace = F,
                     prob = mdoriginal[size_window_A[i]:size_window_B[i]])

      temp <- index(chunk) %in% m_a2

      chunk[temp == TRUE] <- NA

      chunks <- append(chunks, chunk)
    }
    return(chunks)
  }

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

  #Check if start methods are in list of avaliable methods
  check <- start.method %in% df_of_methods$methods

  if (all(check) == TRUE){
    for (i in 1:length(start.method)) {
      # First imputation
      first_imputed <- eval(parse(text = df_of_methods$formulas_x_ts[df_of_methods$methods == start.method[i]]))

      # Get MD original
      mdoriginal <- get_mdoriginal(x = x_ts, y = first_imputed)

      for (k in 1:iterations) {
        # Put new simulated MD
        newmdsimulation <- split_arbitrary(x = first_imputed, percentmd = percentmd, segments = segments, mdoriginal = mdoriginal)

        for (j in 1:length(methods)) {
          #Check if selected methods are in list of avaliable methods
          check2 <- methods %in% df_of_methods$methods

          if (all(check) == TRUE){ # if all works correctly
            # Actual imputation
            actual_imputation <- eval(parse(text = df_of_methods$formulas_actual_ts[df_of_methods$methods == methods[j]]))#Here it goes good

            # Get scores
            rmse <- rmse(coredata(first_imputed),coredata(actual_imputation))
            cor <- cor(coredata(first_imputed),coredata(actual_imputation))^2
            mase <- mase(coredata(first_imputed),coredata(actual_imputation))
            smape <- smape(coredata(first_imputed),coredata(actual_imputation))

            #Storage scores temporarly
            tempresults <- data.frame("start.method" = start.method[i],
                                      "actual.method" = methods[j],
                                      "rmse" = rmse,
                                      "cor" = cor,
                                      "mase" = mase,
                                      "smape"= smape)

            # Append to final results
            results <- bind_rows(results, tempresults)
          }
          else {
            print(paste0("The methods '",
                         paste(as.character(methods[which(!check2)]),
                               collapse = ", "),
                         "' in actual.method parameter, are not in the list of available options"))
            }
        }
      }
    }
  }
  else {
    print(paste0("The methods '",
                 paste(as.character(start.method[which(!check)]),
                       collapse = ", "),
                 "' in start.method parameter, are not in the list of available options"))
  }
  summary_results <- results %>%
    group_by(start.method, actual.method) %>%
    summarise(mean_rmse = mean(rmse),
              std_rmse = sd(rmse),
              mean_cor = mean(cor),
              std_cor = sd(cor),
              mean_mase = mean(mase),
              std_mase = sd(mase),
              mean_smape = mean(smape),
              std_smape = sd(smape))
  list_results <- list(results, summary_results)
  return(list_results)
}
