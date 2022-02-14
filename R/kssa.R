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
#' @export

kssa <- function(x_ts, #Time-series
                 start.method, # Can select various
                 methods, # Can select various
                 segments = 6, #Number of segments to ts be divided
                 iterations = 30, #Number of replicates
                 percentmd = .3, #Percentage of new MD in simulations
                 metric) { #Metrics to evaluate results

  # #' @import forecast

  get_mdoriginal <- function(x, y){
    na_original <- which(is.na(x))
    weight_index <- rep(1, length(y))
    weight_index[na_original] <- 0
    return(weight_index)
  }

  split_arbitrary <- function(x, percentmd, n_parts, mdoriginal){
    size_window_B <- seq(from = round(length(x)/n_parts),
                         to = length(x),
                         by = round(length(x)/n_parts)+1)
    size_window_A <- size_window_B - (n_parts + 1)
    size_window_B <- c(size_window_B, length(x))
    size_window_A <- c(size_window_A, size_window_B[length(size_window_B)-1]+1)
    size_window_A[1] <- 1

    index_time <- index(x)

    chunks <- c()

    for (i in 1:n_parts) {
      chunk <- window(x = x, start = index_time[size_window_A[i]],
                      end = index_time[size_window_B[i]])

      m_a2 <- sample(x = chunk, size = (round(length(chunk)*percentmd)), replace = F,
                     prob = mdoriginal[size_window_A[i]:size_window_B[i]])

      chunk[index(m_a2)] <- NA

      chunks <- append(chunks, chunk)
    }
    return(chunks)
  }


# 1. IMPUT DATA FOR DIFFERENT SPECIES SEPEC ####

  df_of_methods <- data.frame("methods" =c("auto.arima", "StructTS", "linear",
                                           "spline", "stine", "simple", "malinear",
                                           "exponential", "kalman", "nalocf", "decomp"),
                              "formulas" = c("na_kalman(x_ts,model='auto.arima',smooth = TRUE,nit = -1)",
                                             "na_kalman(x_ts,model='StructTS',smooth = TRUE,nit = -1)",
                                             "na_interpolation(x_ts,option='linear')",
                                             "na_interpolation(x_ts,option='spline')",
                                             "na_interpolation(x_ts,option='stine')",
                                             "na_ma(x_ts,k=3,weighting = 'simple')",
                                             "na_ma(x_ts,k=3,weighting = 'linear')",
                                             "na_ma(x_ts,k=3,weighting = 'exponential')",
                                             "na_seadec(x_ts,algorithm = 'kalman')",
                                             "na_locf(x_ts,option = 'locf',na_remaining = 'rev')",
                                             "na.interp(x_ts)")
                              )

  for (i in 1:length(df_of_methods$methods)) {
    check <- start.method %in% df_of_methods$methods
    if (all(check) == TRUE){
      for (j in 1:length(start.method)) {
        first_imputed <- eval(parse(text = df_of_methods$formulas[df_of_methods$methods == start.method[j]]))
#Continue... [Here] <- El siguiente paso es introducir mi funcion, luego hay que ver el resultado contemplando una salida multiple de metodo inicial y de metodo secundario
      }
    }
    else {
      print(paste0("The methods '",
                   paste(as.character(start.method[which(!check)]),
                         collapse = ", "),
                   "' in start.method parameter, are not in the list of available options"))

    }
  }

# This is for the starting imputations
if (start.method == "interpolation") {
  start_data <- imputeTS::na.interpolation(x_ts)
}
else if (start.method == "auto.arima") {
  start_data <- imputeTS::na.interpolation(x_ts)
}
  else if (start.method == "auto.arima") {
    start_data <- imputeTS::na.interpolation(x_ts)
  }

#x_imput_aa=na_kalman(x_ts,model="auto.arima",smooth = TRUE,nit = -1) #m?todo Structured Model
#x_imput_st=na_kalman(x_ts,model="StructTS",smooth = TRUE,nit = -1) #m?todo auto.arima
#x_imput_itl=na_interpolation(x_ts,option="linear") #m?todo interpolacion lineal
#x_imput_its=na_interpolation(x_ts,option="spline") #m?todo interpolaci?n spline
#x_imput_itst=na_interpolation(x_ts,option="stine") #m?todo interpolaci?n spline
#x_imput_mas=na_ma(x_ts,k=3,weighting = "simple") #m?todo moving average simple
#x_imput_mal=na_ma(x_ts,k=3,weighting = "linear") #m?todo moving average linear
#x_imput_mae=na_ma(x_ts,k=3,weighting = "exponential") #m?todo moving average exponential
#x_imput_seadec=na_seadec(x_ts,algorithm = "kalman") #m?todo de imputacion kalman en serie estacionalmente descompuesta
#x_imput_locf=na_locf(x_ts,option = "locf",na_remaining = "rev") #metodo de ultima observacion hacia delante
#x_imput_stl=na.interp(x_ts) #metodo de descomposici?n STL con interpolaci?n lineal




#2. VALIDATION WHEN SIMULATING THE NA'S STRUCTURE OF DATA IN CHUNKSSS ####


#here we need a code line to split the time series and specify where to put simulation windows

for (i in 1:iterations) {
  # 1. Simulate new MD
  brayans_function()

  # 2. Make imputations on simulated MD
  #y_impute <- na_kalman(y_na,model="StructTS",nit = -1)

  for (i in methods) {

    if (i == "intpoltaion")
  }
   y_impute <- na_interpolation(y_na)

  # 3. Calculate performance metrics
  nas <-sum(is.na(y_na))
  rmse <- rmse(serie.temporal,y_impute)
  cor <- cor(serie.temporal,y_impute)^2
  mase <- mase(serie.temporal,y_impute)
  smape <- smape(serie.temporal,y_impute)


  # 4. Store metrics in a temporal table
  temporal1 <- data.frame(nas,rmse, cor,mase,smape)
  colnames(temporal1)=c("nas","rmse","cor","mase","smape")

  # 5. Assign metrics in the temporal table to a final table
  results_table <- bind_rows(results_table, temporal1)
}

#Create a data.frame with final table
return(results_table)

}

