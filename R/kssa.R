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

kssa <- function(x_ts, start.method, segments = 6, iterations, metric) {

  # #' @import forecast

## Long code KSSA


############### 1. IMPUT DATA FOR DIFFERENT SPECIES SEPEC #############



# This is for the starting imputations

x_imput_aa=na_kalman(x_ts,model="auto.arima",smooth = TRUE,nit = -1) #m?todo Structured Model
x_imput_st=na_kalman(x_ts,model="StructTS",smooth = TRUE,nit = -1) #m?todo auto.arima
x_imput_itl=na_interpolation(x_ts,option="linear") #m?todo interpolacion lineal
x_imput_its=na_interpolation(x_ts,option="spline") #m?todo interpolaci?n spline
x_imput_itst=na_interpolation(x_ts,option="stine") #m?todo interpolaci?n spline
x_imput_mas=na_ma(x_ts,k=3,weighting = "simple") #m?todo moving average simple
x_imput_mal=na_ma(x_ts,k=3,weighting = "linear") #m?todo moving average linear
x_imput_mae=na_ma(x_ts,k=3,weighting = "exponential") #m?todo moving average exponential
x_imput_seadec=na_seadec(x_ts,algorithm = "kalman") #m?todo de imputacion kalman en serie estacionalmente descompuesta
x_imput_locf=na_locf(x_ts,option = "locf",na_remaining = "rev") #metodo de ultima observacion hacia delante
#x_imput_stl=na.interp(x_ts) #metodo de descomposici?n STL con interpolaci?n lineal




################################################################################
#################### 2. VALIDATION WHEN SIMULATING THE NA'S STRUCTURE OF DATA IN CHUNKSSS ###############
#Con los mismos datos pero los imputados
#Hacerlo para cada m?todo por separado
#HAcerlo para cada especie por separado
#cambiar especie y m?todo en <-serie.temporal



#times to repeat the loop (same for all)
iterations=1000

#here we need a code line to split the time series and specify where to put simulation windows

gap1=seq(1,12,1);gap2=seq(13,25,1);gap3=seq(26,38,1);gap4=seq(39,51,1);gap5=seq(52,64,1);gap6=seq(65,77,1);gap7=seq(78,90,1);gap8=seq(91,101,1)




#### FOR STRUCTS

results_table <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(results_table) <- c('nas','rmse', 'cor','mase','smape')
serie.temporal <- x_imput  # CAMBIAR AQUI la serie con la que voy a trabajar
# Eliminar la serie de datos inventada
results_table <- results_table[-1,]

for (i in 1:iterations) {
  # 1. Simulate new MD
  y_na <- sjmisc::set_na(x_imput, na=  c(x_imput[sample(gap1,1):sample(gap1,1)],x_imput[sample(gap2,1):sample(gap2,1)],x_imput[sample(gap3,1):sample(gap3,1)],x_imput[sample(gap4,1):sample(gap4,1)],x_imput[sample(gap5,1):sample(gap5,1)],x_imput[sample(gap6,1):sample(gap6,1)],x_imput[sample(gap7,1):sample(gap7,1)],x_imput[sample(gap8,1):sample(gap8,1)]))

  # 2. Make imputations on simulated MD
  #y_impute <- na_kalman(y_na,model="StructTS",nit = -1)
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

