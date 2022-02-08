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
#'
#' @export

kssa <- function(x_ts, start.method = "interpolation", methods = c("auto.arima", "interpolation"), segments = 6, iterations = 30, metric) {

  # #' @import forecast


############### 1. IMPUT DATA FOR DIFFERENT SPECIES SEPEC #############

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




################################################################################
#################### 2. VALIDATION WHEN SIMULATING THE NA'S STRUCTURE OF DATA IN CHUNKSSS ###############


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

