#' @title Add two numbers
#'
#' @description Function to add two numbers
#' @param x A number
#' @param y A number
#' @return The sum of \code{x} and \code{y}
#' @examples
#' add(1, 1)
#' add(10, 1)
#'
#' @export

add <- function(x, y) {

  z <- x + y
  print("first function")

  return(z)
}

## Long code KSSA

library(sjmisc)
library(imputeTS)
library(Metrics)
library(dplyr)
library(snpar)
library(BaylorEdPsych)
library(mvnmle)
library(forecast)
setwd("C:/Users/USUARIO/Dropbox/Mis Art?culos/Art?culo Marlon/Datos")#PC1


############### 1. IMPUT DATA FOR DIFFERENT SPECIES SEPEC #############


dfull=read.table("sepec_full_series.txt",h=T)
#extract data for the full pacifico
dfullpac=dfull[dfull$site=="pacifico",]
#replace zeroes by NA because we do not know if zeroes are really no landings or abscent data
dfullpac[dfullpac==0]=NA
#check table
View(dfullpac)
#check NA's for each species
as.data.frame(summary(dfullpac))


#convert to time series ts
bagre_ts=ts(dfullpac$bagre,start=c(2012,1),end=c(2020,5),frequency = 12)
camaron_blanco_ts=ts(dfullpac$camaron_blanco,start=c(2012,1),end=c(2020,5),frequency = 12)
camaron_pomadilla_ts=ts(dfullpac$camaron_pomadilla,start=c(2012,1),end=c(2020,5),frequency = 12)
cherna_rosada_ts=ts(dfullpac$cherna_rosada,start=c(2012,1),end=c(2020,5),frequency = 12)
corvina_ts=ts(dfullpac$corvina,start=c(2012,1),end=c(2020,5),frequency = 12)
dorado_ts=ts(dfullpac$dorado,start=c(2012,1),end=c(2020,5),frequency = 12)
jaiba_ts=ts(dfullpac$jaiba,start=c(2012,1),end=c(2020,5),frequency = 12)
merluza_ts=ts(dfullpac$merluza,start=c(2012,1),end=c(2020,5),frequency = 12)
pargo_lunarejo_ts=ts(dfullpac$pargo_lunarejo,start=c(2012,1),end=c(2020,5),frequency = 12)
pargo_rojo_ts=ts(dfullpac$pargo_rojo,start=c(2012,1),end=c(2020,5),frequency = 12)
pelada_ts=ts(dfullpac$pelada,start=c(2012,1),end=c(2020,5),frequency = 12)
picua_ts=ts(dfullpac$picua,start=c(2012,1),end=c(2020,5),frequency = 12)
sierra_ts=ts(dfullpac$sierra,start=c(2012,1),end=c(2020,5),frequency = 12)
camaron_titi_ts=ts(dfullpac$camaron_titi,start=c(2012,1),end=c(2020,5),frequency = 12)
yft_ts=ts(dfullpac$yft,start=c(2012,1),end=c(2020,5),frequency = 12)
skj_ts=ts(dfullpac$skp,start=c(2012,1),end=c(2020,5),frequency = 12)

#distribution graph of NA's for each species

ggplot_na_distribution(bagre_ts)
ggplot_na_distribution(camaron_blanco_ts)
ggplot_na_distribution(camaron_pomadilla_ts)
ggplot_na_distribution(corvina_ts)
ggplot_na_distribution(dorado_ts)
ggplot_na_distribution(sierra_ts)
ggplot_na_distribution(yft_ts)
ggplot_na_distribution(skj_ts)



#imput missing data to bagre and create time series with NA's deleted
ggplot_na_distribution(cherna_rosada_ts)
statsNA(corvina_ts) #24 NA, 23.8%
bagre.runs=ifelse(bagre_ts >= 0, 1, 0);bagre.runs[is.na(bagre.runs)] <- 0 #convert time series to a binary vector where numbers = 1 and NA = 0
runs.test(bagre.runs) # runs test to check if NAs appear at random
bagre_no_na=na_remove(bagre_ts) #remove na's from bagre_ts
bagre_imput_aa=na_kalman(bagre_ts,model="auto.arima",smooth = TRUE,nit = -1) #m?todo Structured Model
bagre_imput_st=na_kalman(bagre_ts,model="StructTS",smooth = TRUE,nit = -1) #m?todo auto.arima
bagre_imput_itl=na_interpolation(bagre_ts,option="linear") #m?todo interpolacion lineal
bagre_imput_its=na_interpolation(bagre_ts,option="spline") #m?todo interpolaci?n spline
bagre_imput_itst=na_interpolation(bagre_ts,option="stine") #m?todo interpolaci?n spline
bagre_imput_mas=na_ma(bagre_ts,k=3,weighting = "simple") #m?todo moving average simple
bagre_imput_mal=na_ma(bagre_ts,k=3,weighting = "linear") #m?todo moving average linear
bagre_imput_mae=na_ma(bagre_ts,k=3,weighting = "exponential") #m?todo moving average exponential
bagre_imput_seadec=na_seadec(bagre_ts,algorithm = "kalman") #m?todo de imputacion kalman en serie estacionalmente descompuesta
bagre_imput_locf=na_locf(bagre_ts,option = "locf",na_remaining = "rev") #metodo de ultima observacion hacia delante
bagre_imput_stl=na.interp(bagre_ts) #metodo de descomposici?n STL con interpolaci?n lineal



#graficar
plot.ts(bagre_ts,lwd=16,col="gray")
lines(bagre_imput_aa,lwd=2,col="red")
lines(bagre_imput_st,lwd=2,col="blue")
lines(bagre_imput_itl,lwd=2,col="orange")
lines(bagre_imput_its,lwd=2,col="purple")
lines(bagre_imput_itst,lwd=2,col="green")
lines(bagre_imput_mas,lwd=2,col="seashell3")
lines(bagre_imput_mal,lwd=2,col="salmon1")
lines(bagre_imput_mae,lwd=2,col="blue4")
lines(bagre_imput_seadec,lwd=2,col="orange")

#imput missing data to camaron_blanco
camaron_blanco_no_na=na_remove(camaron_blanco_ts)
ggplot_na_distribution(camaron_blanco_ts)
statsNA(camaron_blanco_ts) #24 NA, 23.8%
camaron_blanco.runs=ifelse(camaron_blanco_ts >= 0, 1, 0);camaron_blanco.runs[is.na(camaron_blanco.runs)] <- 0 #convert time series to a binary vector where numbers = 1 and NA = 0
runs.test(camaron_blanco.runs) # runs test to check if NAs appear at random
camaron_blanco_imput_aa=na_kalman(camaron_blanco_ts,model="auto.arima",smooth = TRUE,nit = -1) #m?todo Structured Model
camaron_blanco_imput_st=na_kalman(camaron_blanco_ts,model="StructTS",smooth = TRUE,nit = -1) #m?todo auto.arima
camaron_blanco_imput_itl=na_interpolation(camaron_blanco_ts,option="linear") #m?todo interpolacion lineal
camaron_blanco_imput_its=na_interpolation(camaron_blanco_ts,option="spline") #m?todo interpolaci?n spline
camaron_blanco_imput_itst=na_interpolation(camaron_blanco_ts,option="stine") #m?todo interpolaci?n spline
camaron_blanco_imput_mas=na_ma(camaron_blanco_ts,k=3,weighting = "simple") #m?todo moving average simple
camaron_blanco_imput_mal=na_ma(camaron_blanco_ts,k=3,weighting = "linear") #m?todo moving average linear
camaron_blanco_imput_mae=na_ma(camaron_blanco_ts,k=3,weighting = "exponential") #m?todo moving average exponential
camaron_blanco_imput_seadec=na_seadec(camaron_blanco_ts,algorithm = "kalman") #m?todo de imputacion kalman en serie estacionalmente descompuesta
camaron_blanco_imput_locf=na_locf(camaron_blanco_ts,option = "locf",na_remaining = "rev") #metodo de ultima observacion hacia delante
camaron_blanco_imput_stl=na.interp(camaron_blanco_ts) #metodo de ultima observacion hacia delante


#graficar
plot.ts(camaron_blanco_ts,lwd=16,col="gray")
lines(camaron_blanco_imput_aa,lwd=2,col="red")
lines(camaron_blanco_imput_st,lwd=2,col="blue")
lines(camaron_blanco_imput_itl,lwd=2,col="orange")
lines(camaron_blanco_imput_its,lwd=2,col="purple")
lines(camaron_blanco_imput_itst,lwd=2,col="green")
lines(camaron_blanco_imput_mas,lwd=2,col="seashell3")
lines(camaron_blanco_imput_mal,lwd=2,col="salmon1")
lines(camaron_blanco_imput_mae,lwd=2,col="blue4")

#imput missing data to camaron_pomadilla
camaron_pomadilla_no_na=na_remove(camaron_pomadilla_ts)
ggplot_na_distribution(camaron_pomadilla_ts)
statsNA(camaron_pomadilla_ts) #24 NA, 23.8%
camaron_pomadilla.runs=ifelse(camaron_pomadilla_ts >= 0, 1, 0);camaron_pomadilla.runs[is.na(camaron_pomadilla.runs)] <- 0 #convert time series to a binary vector where numbers = 1 and NA = 0
runs.test(camaron_pomadilla.runs) # runs test to check if NAs appear at random
camaron_pomadilla_imput_aa=na_kalman(camaron_pomadilla_ts,model="auto.arima",smooth = TRUE,nit = -1) #m?todo Structured Model
camaron_pomadilla_imput_st=na_kalman(camaron_pomadilla_ts,model="StructTS",smooth = TRUE,nit = -1) #m?todo auto.arima
camaron_pomadilla_imput_itl=na_interpolation(camaron_pomadilla_ts,option="linear") #m?todo interpolacion lineal
camaron_pomadilla_imput_its=na_interpolation(camaron_pomadilla_ts,option="spline") #m?todo interpolaci?n spline
camaron_pomadilla_imput_itst=na_interpolation(camaron_pomadilla_ts,option="stine") #m?todo interpolaci?n spline
camaron_pomadilla_imput_mas=na_ma(camaron_pomadilla_ts,k=3,weighting = "simple") #m?todo moving average simple
camaron_pomadilla_imput_mal=na_ma(camaron_pomadilla_ts,k=3,weighting = "linear") #m?todo moving average linear
camaron_pomadilla_imput_mae=na_ma(camaron_pomadilla_ts,k=3,weighting = "exponential") #m?todo moving average exponential
camaron_pomadilla_imput_seadec=na_seadec(camaron_pomadilla_ts,algorithm = "kalman") #m?todo de imputacion kalman en serie estacionalmente descompuesta
camaron_pomadilla_imput_locf=na_locf(camaron_pomadilla_ts,option = "locf",na_remaining = "rev") #metodo de ultima observacion hacia delante
camaron_pomadilla_imput_stl=na.interp(camaron_pomadilla_ts) #metodo de descomposici?n STL con interpolaci?n lineal


#graficar
plot.ts(camaron_pomadilla_ts,lwd=16,col="gray")
lines(camaron_pomadilla_imput_aa,lwd=2,col="red")
lines(camaron_pomadilla_imput_st,lwd=2,col="blue")
lines(camaron_pomadilla_imput_itl,lwd=2,col="orange")
lines(camaron_pomadilla_imput_its,lwd=2,col="purple")
lines(camaron_pomadilla_imput_itst,lwd=2,col="green")
lines(camaron_pomadilla_imput_mas,lwd=2,col="seashell3")
lines(camaron_pomadilla_imput_mal,lwd=2,col="salmon1")
lines(camaron_pomadilla_imput_mae,lwd=2,col="blue4")

#imput missing data to cherna_rosada
ggplot_na_distribution(cherna_rosada_ts)
statsNA(cherna_rosada_ts) #24 NA, 23.8%
cherna_rosada.runs=ifelse(cherna_rosada_ts >= 0, 1, 0);cherna_rosada.runs[is.na(cherna_rosada.runs)] <- 0 #convert time series to a binary vector where numbers = 1 and NA = 0
runs.test(cherna_rosada.runs) # runs test to check if NAs appear at random
cherna_rosada_imput_aa=na_kalman(cherna_rosada_ts,model="auto.arima",smooth = TRUE,nit = -1) #m?todo Structured Model
cherna_rosada_imput_st=na_kalman(cherna_rosada_ts,model="StructTS",smooth = TRUE,nit = -1) #m?todo auto.arima
cherna_rosada_imput_itl=na_interpolation(cherna_rosada_ts,option="linear") #m?todo interpolacion lineal
cherna_rosada_imput_its=na_interpolation(cherna_rosada_ts,option="spline") #m?todo interpolaci?n spline
cherna_rosada_imput_itst=na_interpolation(cherna_rosada_ts,option="stine") #m?todo interpolaci?n spline
cherna_rosada_imput_mas=na_ma(cherna_rosada_ts,k=3,weighting = "simple") #m?todo moving average simple
cherna_rosada_imput_mal=na_ma(cherna_rosada_ts,k=3,weighting = "linear") #m?todo moving average linear
cherna_rosada_imput_mae=na_ma(cherna_rosada_ts,k=3,weighting = "exponential") #m?todo moving average exponential
cherna_rosada_imput_seadec=na_seadec(cherna_rosada_ts,algorithm = "kalman") #m?todo de imputacion kalman en serie estacionalmente descompuesta
cherna_rosada_imput_locf=na_locf(cherna_rosada_ts,option = "locf",na_remaining = "rev") #metodo de ultima observacion hacia delante
cherna_rosada_imput_stl=na.interp(cherna_rosada_ts) #metodo de descomposici?n STL con interpolaci?n lineal


#graficar
plot.ts(cherna_rosada_ts,lwd=16,col="gray")
lines(cherna_rosada_imput_aa,lwd=2,col="red")
lines(cherna_rosada_imput_st,lwd=2,col="blue")
lines(cherna_rosada_imput_itl,lwd=2,col="orange")
lines(cherna_rosada_imput_its,lwd=2,col="purple")
lines(cherna_rosada_imput_itst,lwd=2,col="green")
lines(cherna_rosada_imput_mas,lwd=2,col="seashell3")
lines(cherna_rosada_imput_mal,lwd=2,col="salmon1")
lines(cherna_rosada_imput_mae,lwd=2,col="blue4")

#imput missing data to corvina
ggplot_na_distribution(corvina_ts)
statsNA(corvina_ts) #24 NA, 23.8%
corvina.runs=ifelse(corvina_ts >= 0, 1, 0);corvina.runs[is.na(corvina.runs)] <- 0 #convert time series to a binary vector where numbers = 1 and NA = 0
runs.test(corvina.runs) # runs test to check if NAs appear at random
corvina_imput_aa=na_kalman(corvina_ts,model="auto.arima",smooth = TRUE,nit = -1) #m?todo Structured Model
corvina_imput_st=na_kalman(corvina_ts,model="StructTS",smooth = TRUE,nit = -1) #m?todo auto.arima
corvina_imput_itl=na_interpolation(corvina_ts,option="linear") #m?todo interpolacion lineal
corvina_imput_its=na_interpolation(corvina_ts,option="spline") #m?todo interpolaci?n spline
corvina_imput_itst=na_interpolation(corvina_ts,option="stine") #m?todo interpolaci?n spline
corvina_imput_mas=na_ma(corvina_ts,k=3,weighting = "simple") #m?todo moving average simple
corvina_imput_mal=na_ma(corvina_ts,k=3,weighting = "linear") #m?todo moving average linear
corvina_imput_mae=na_ma(corvina_ts,k=3,weighting = "exponential") #m?todo moving average exponential
corvina_imput_seadec=na_seadec(corvina_ts,algorithm = "kalman") #m?todo de imputacion kalman en serie estacionalmente descompuesta
corvina_imput_locf=na_locf(corvina_ts,option = "locf",na_remaining = "rev") #metodo de ultima observacion hacia delante
corvina_imput_stl=na.interp(corvina_ts) #metodo de descomposici?n STL con interpolaci?n lineal


#graficar
plot.ts(corvina_ts,lwd=16,col="gray")
lines(corvina_imput_aa,lwd=2,col="red")
lines(corvina_imput_st,lwd=2,col="blue")
lines(corvina_imput_itl,lwd=2,col="orange")
lines(corvina_imput_its,lwd=2,col="purple")
lines(corvina_imput_itst,lwd=2,col="green")
lines(corvina_imput_mas,lwd=2,col="seashell3")
lines(corvina_imput_mal,lwd=2,col="salmon1")
lines(corvina_imput_mae,lwd=2,col="blue4")


#imput missing data to dorado
statsNA(dorado_ts) #24 NA, 23.8%
dorado_imput_aa=na_kalman(dorado_ts,model="auto.arima",smooth = TRUE,nit = -1) #m?todo Structured Model
dorado_imput_st=na_kalman(dorado_ts,model="StructTS",smooth = TRUE,nit = -1) #m?todo auto.arima
dorado_imput_itl=na_interpolation(dorado_ts,option="linear") #m?todo interpolacion lineal
dorado_imput_its=na_interpolation(dorado_ts,option="spline") #m?todo interpolaci?n spline
dorado_imput_itst=na_interpolation(dorado_ts,option="stine") #m?todo interpolaci?n spline
dorado_imput_mas=na_ma(dorado_ts,k=3,weighting = "simple") #m?todo moving average simple
dorado_imput_mal=na_ma(dorado_ts,k=3,weighting = "linear") #m?todo moving average linear
dorado_imput_mae=na_ma(dorado_ts,k=3,weighting = "exponential") #m?todo moving average exponential
dorado_imput_seadec=na_seadec(dorado_ts,algorithm = "kalman") #m?todo de imputacion kalman en serie estacionalmente descompuesta
dorado_imput_locf=na_locf(dorado_ts,option = "locf",na_remaining = "rev") #metodo de ultima observacion hacia delante

#graficar
plot.ts(dorado_ts,lwd=16,col="gray")
lines(dorado_imput_aa,lwd=2,col="red")
lines(dorado_imput_st,lwd=2,col="blue")
lines(dorado_imput_itl,lwd=2,col="orange")
lines(dorado_imput_its,lwd=2,col="purple")
lines(dorado_imput_itst,lwd=2,col="green")
lines(dorado_imput_mas,lwd=2,col="seashell3")
lines(dorado_imput_mal,lwd=2,col="salmon1")
lines(dorado_imput_mae,lwd=2,col="blue4")################# VALIDACION DE LA IMPUTACION #################################

#imput missing data to jaiba
statsNA(jaiba_ts) #24 NA, 23.8%
jaiba_imput_aa=na_kalman(jaiba_ts,model="auto.arima",smooth = TRUE,nit = -1) #m?todo Structured Model
jaiba_imput_st=na_kalman(jaiba_ts,model="StructTS",smooth = TRUE,nit = -1) #m?todo auto.arima
jaiba_imput_itl=na_interpolation(jaiba_ts,option="linear") #m?todo interpolacion lineal
jaiba_imput_its=na_interpolation(jaiba_ts,option="spline") #m?todo interpolaci?n spline
jaiba_imput_itst=na_interpolation(jaiba_ts,option="stine") #m?todo interpolaci?n spline
jaiba_imput_mas=na_ma(jaiba_ts,k=3,weighting = "simple") #m?todo moving average simple
jaiba_imput_mal=na_ma(jaiba_ts,k=3,weighting = "linear") #m?todo moving average linear
jaiba_imput_mae=na_ma(jaiba_ts,k=3,weighting = "exponential") #m?todo moving average exponential
jaiba_imput_seadec=na_seadec(jaiba_ts,algorithm = "kalman") #m?todo de imputacion kalman en serie estacionalmente descompuesta
jaiba_imput_locf=na_locf(jaiba_ts,option = "locf",na_remaining = "rev") #metodo de ultima observacion hacia delante


#graficar
plot.ts(jaiba_ts,lwd=16,col="gray")
lines(jaiba_imput_aa,lwd=2,col="red")
lines(jaiba_imput_st,lwd=2,col="blue")
lines(jaiba_imput_itl,lwd=2,col="orange")
lines(jaiba_imput_its,lwd=2,col="purple")
lines(jaiba_imput_itst,lwd=2,col="green")
lines(jaiba_imput_mas,lwd=2,col="seashell3")
lines(jaiba_imput_mal,lwd=2,col="salmon1")
lines(jaiba_imput_mae,lwd=2,col="blue4")

#imput missing data to merluza
statsNA(merluza_ts) #24 NA, 23.8%
merluza_imput_aa=na_kalman(merluza_ts,model="auto.arima",smooth = TRUE,nit = -1) #m?todo Structured Model
merluza_imput_st=na_kalman(merluza_ts,model="StructTS",smooth = TRUE,nit = -1) #m?todo auto.arima
merluza_imput_itl=na_interpolation(merluza_ts,option="linear") #m?todo interpolacion lineal
merluza_imput_its=na_interpolation(merluza_ts,option="spline") #m?todo interpolaci?n spline
merluza_imput_itst=na_interpolation(merluza_ts,option="stine") #m?todo interpolaci?n spline
merluza_imput_mas=na_ma(merluza_ts,k=3,weighting = "simple") #m?todo moving average simple
merluza_imput_mal=na_ma(merluza_ts,k=3,weighting = "linear") #m?todo moving average linear
merluza_imput_mae=na_ma(merluza_ts,k=3,weighting = "exponential") #m?todo moving average exponential
merluza_imput_seadec=na_seadec(merluza_ts,algorithm = "kalman") #m?todo de imputacion kalman en serie estacionalmente descompuesta
merluza_imput_locf=na_locf(merluza_ts,option = "locf",na_remaining = "rev") #metodo de ultima observacion hacia delante


#graficar
plot.ts(merluza_ts,lwd=16,col="gray")
lines(merluza_imput_aa,lwd=2,col="red")
lines(merluza_imput_st,lwd=2,col="blue")
lines(merluza_imput_itl,lwd=2,col="orange")
lines(merluza_imput_its,lwd=2,col="purple")
lines(merluza_imput_itst,lwd=2,col="green")
lines(merluza_imput_mas,lwd=2,col="seashell3")
lines(merluza_imput_mal,lwd=2,col="salmon1")
lines(merluza_imput_mae,lwd=2,col="blue4")

#imput missing data to pargo_lunarejo
statsNA(pargo_lunarejo_ts) #24 NA, 23.8%
pargo_lunarejo_imput_aa=na_kalman(pargo_lunarejo_ts,model="auto.arima",smooth = TRUE,nit = -1) #m?todo Structured Model
pargo_lunarejo_imput_st=na_kalman(pargo_lunarejo_ts,model="StructTS",smooth = TRUE,nit = -1) #m?todo auto.arima
pargo_lunarejo_imput_itl=na_interpolation(pargo_lunarejo_ts,option="linear") #m?todo interpolacion lineal
pargo_lunarejo_imput_its=na_interpolation(pargo_lunarejo_ts,option="spline") #m?todo interpolaci?n spline
pargo_lunarejo_imput_itst=na_interpolation(pargo_lunarejo_ts,option="stine") #m?todo interpolaci?n spline
pargo_lunarejo_imput_mas=na_ma(pargo_lunarejo_ts,k=3,weighting = "simple") #m?todo moving average simple
pargo_lunarejo_imput_mal=na_ma(pargo_lunarejo_ts,k=3,weighting = "linear") #m?todo moving average linear
pargo_lunarejo_imput_mae=na_ma(pargo_lunarejo_ts,k=3,weighting = "exponential") #m?todo moving average exponential
pargo_lunarejo_imput_seadec=na_seadec(pargo_lunarejo_ts,algorithm = "kalman") #m?todo de imputacion kalman en serie estacionalmente descompuesta
pargo_lunarejo_imput_locf=na_locf(pargo_lunarejo_ts,option = "locf",na_remaining = "rev") #metodo de ultima observacion hacia delante


#graficar
plot.ts(pargo_lunarejo_ts,lwd=16,col="gray")
lines(pargo_lunarejo_imput_aa,lwd=2,col="red")
lines(pargo_lunarejo_imput_st,lwd=2,col="blue")
lines(pargo_lunarejo_imput_itl,lwd=2,col="orange")
lines(pargo_lunarejo_imput_its,lwd=2,col="purple")
lines(pargo_lunarejo_imput_itst,lwd=2,col="green")
lines(pargo_lunarejo_imput_mas,lwd=2,col="seashell3")
lines(pargo_lunarejo_imput_mal,lwd=2,col="salmon1")
lines(pargo_lunarejo_imput_mae,lwd=2,col="blue4")

#imput missing data to pargo_rojo
ggplot_na_distribution(pargo_rojo_ts)
statsNA(pargo_rojo_ts) #24 NA, 23.8%

pargo_rojo_imput_aa=na_kalman(pargo_rojo_ts,model="auto.arima",smooth = TRUE,nit = -1) #m?todo Structured Model
pargo_rojo_imput_st=na_kalman(pargo_rojo_ts,model="StructTS",smooth = TRUE,nit = -1) #m?todo auto.arima
pargo_rojo_imput_itl=na_interpolation(pargo_rojo_ts,option="linear") #m?todo interpolacion lineal
pargo_rojo_imput_its=na_interpolation(pargo_rojo_ts,option="spline") #m?todo interpolaci?n spline
pargo_rojo_imput_itst=na_interpolation(pargo_rojo_ts,option="stine") #m?todo interpolaci?n spline
pargo_rojo_imput_mas=na_ma(pargo_rojo_ts,k=3,weighting = "simple") #m?todo moving average simple
pargo_rojo_imput_mal=na_ma(pargo_rojo_ts,k=3,weighting = "linear") #m?todo moving average linear
pargo_rojo_imput_mae=na_ma(pargo_rojo_ts,k=3,weighting = "exponential") #m?todo moving average exponential
pargo_rojo_imput_seadec=na_seadec(pargo_rojo_ts,algorithm = "kalman") #m?todo de imputacion kalman en serie estacionalmente descompuesta
pargo_rojo_imput_locf=na_locf(pargo_rojo_ts,option = "locf",na_remaining = "rev") #metodo de ultima observacion hacia delante
pargo_rojo_imput_stl=na.interp(pargo_rojo_ts) #metodo de descomposici?n STL con interpolaci?n lineal


#graficar
plot.ts(pargo_rojo_ts,lwd=16,col="gray")
lines(pargo_rojo_imput_aa,lwd=2,col="red")
lines(pargo_rojo_imput_st,lwd=2,col="blue")
lines(pargo_rojo_imput_itl,lwd=2,col="orange")
lines(pargo_rojo_imput_its,lwd=2,col="purple")
lines(pargo_rojo_imput_itst,lwd=2,col="green")
lines(pargo_rojo_imput_mas,lwd=2,col="seashell3")
lines(pargo_rojo_imput_mal,lwd=2,col="salmon1")
lines(pargo_rojo_imput_mae,lwd=2,col="blue4")#EN LOS MISMOS DATOS DE SEPEC

#imput missing data to pelada
ggplot_na_distribution(pelada_ts)
statsNA(pelada_ts) #24 NA, 23.8%
pargo_rojo.runs=ifelse(pargo_rojo_ts >= 0, 1, 0);pargo_rojo.runs[is.na(pargo_rojo.runs)] <- 0 #convert time series to a binary vector where numbers = 1 and NA = 0
runs.test(pargo_rojo.runs) # runs test to check if NAs appear at random
pelada_imput_aa=na_kalman(pelada_ts,model="auto.arima",smooth = TRUE,nit = -1) #m?todo Structured Model
pelada_imput_st=na_kalman(pelada_ts,model="StructTS",smooth = TRUE,nit = -1) #m?todo auto.arima
pelada_imput_itl=na_interpolation(pelada_ts,option="linear") #m?todo interpolacion lineal
pelada_imput_its=na_interpolation(pelada_ts,option="spline") #m?todo interpolaci?n spline
pelada_imput_itst=na_interpolation(pelada_ts,option="stine") #m?todo interpolaci?n spline
pelada_imput_mas=na_ma(pelada_ts,k=3,weighting = "simple") #m?todo moving average simple
pelada_imput_mal=na_ma(pelada_ts,k=3,weighting = "linear") #m?todo moving average linear
pelada_imput_mae=na_ma(pelada_ts,k=3,weighting = "exponential") #m?todo moving average exponential
pelada_imput_seadec=na_seadec(pelada_ts,algorithm = "kalman") #m?todo de imputacion kalman en serie estacionalmente descompuesta
pelada_imput_locf=na_locf(pelada_ts,option = "locf",na_remaining = "rev") #metodo de ultima observacion hacia delante
pelada_imput_stl=na.interp(pelada_ts) #metodo de descomposici?n STL con interpolaci?n lineal


#graficar
plot.ts(pelada_ts,lwd=16,col="gray")
lines(pelada_imput_aa,lwd=2,col="red")
lines(pelada_imput_st,lwd=2,col="blue")
lines(pelada_imput_itl,lwd=2,col="orange")
lines(pelada_imput_its,lwd=2,col="purple")
lines(pelada_imput_itst,lwd=2,col="green")
lines(pelada_imput_mas,lwd=2,col="seashell3")
lines(pelada_imput_mal,lwd=2,col="salmon1")
lines(pelada_imput_mae,lwd=2,col="blue4")#LOOP INCLUYENDO N?MERO DE NA'S Y NUMERO DE VECES QUE SE REPITE EL CICLO (LOOP ANIDADO)

#imput missing data to picua
ggplot_na_distribution(picua_ts)
statsNA(picua_ts) #24 NA, 23.8%
picua_imput_aa=na_kalman(picua_ts,model="auto.arima",smooth = TRUE,nit = -1) #m?todo Structured Model
picua_imput_st=na_kalman(picua_ts,model="StructTS",smooth = TRUE,nit = -1) #m?todo auto.arima
picua_imput_itl=na_interpolation(picua_ts,option="linear") #m?todo interpolacion lineal
picua_imput_its=na_interpolation(picua_ts,option="spline") #m?todo interpolaci?n spline
picua_imput_itst=na_interpolation(picua_ts,option="stine") #m?todo interpolaci?n spline
picua_imput_mas=na_ma(picua_ts,k=3,weighting = "simple") #m?todo moving average simple
picua_imput_mal=na_ma(picua_ts,k=3,weighting = "linear") #m?todo moving average linear
picua_imput_mae=na_ma(picua_ts,k=3,weighting = "exponential") #m?todo moving average exponential
picua_imput_seadec=na_seadec(picua_ts,algorithm = "kalman") #m?todo de imputacion kalman en serie estacionalmente descompuesta
picua_imput_locf=na_locf(picua_ts,option = "locf",na_remaining = "rev") #metodo de ultima observacion hacia delante
picua_imput_stl=na.interp(picua_ts) #metodo de descomposici?n STL con interpolaci?n lineal


#graficar
plot.ts(picua_ts,lwd=16,col="gray")
lines(picua_imput_aa,lwd=2,col="red")
lines(picua_imput_st,lwd=2,col="blue")
lines(picua_imput_itl,lwd=2,col="orange")
lines(picua_imput_its,lwd=2,col="purple")
lines(picua_imput_itst,lwd=2,col="green")
lines(picua_imput_mas,lwd=2,col="seashell3")
lines(picua_imput_mal,lwd=2,col="salmon1")
lines(picua_imput_mae,lwd=2,col="blue4")#ESTE HAY QUE HACERLO PARA CADA M?TODO DE IMPUTACION POR SEPARADO (no de una todo)

#imput missing data to sierra
ggplot_na_distribution(sierra_ts)
statsNA(sierra_ts) #24 NA, 23.8%
sierra.runs=ifelse(sierra_ts >= 0, 1, 0);sierra.runs[is.na(sierra.runs)] <- 0 #convert time series to a binary vector where numbers = 1 and NA = 0
runs.test(sierra.runs) # runs test to check if NAs appear at random
sierra_imput_aa=na_kalman(sierra_ts,model="auto.arima",smooth = TRUE,nit = -1) #m?todo Structured Model
sierra_imput_st=na_kalman(sierra_ts,model="StructTS",smooth = TRUE,nit = -1) #m?todo auto.arima
sierra_imput_itl=na_interpolation(sierra_ts,option="linear") #m?todo interpolacion lineal
sierra_imput_its=na_interpolation(sierra_ts,option="spline") #m?todo interpolaci?n spline
sierra_imput_itst=na_interpolation(sierra_ts,option="stine") #m?todo interpolaci?n spline
sierra_imput_mas=na_ma(sierra_ts,k=3,weighting = "simple") #m?todo moving average simple
sierra_imput_mal=na_ma(sierra_ts,k=3,weighting = "linear") #m?todo moving average linear
sierra_imput_mae=na_ma(sierra_ts,k=3,weighting = "exponential") #m?todo moving average exponential
sierra_imput_seadec=na_seadec(sierra_ts,algorithm = "kalman") #m?todo de imputacion kalman en serie estacionalmente descompuesta
sierra_imput_locf=na_locf(sierra_ts,option = "locf",na_remaining = "rev") #metodo de ultima observacion hacia delante
sierra_imput_seadec=na_seadec(sierra_ts,algorithm = "kalman") #m?todo de imputacion kalman en serie estacionalmente descompuesta
sierra_imput_locf=na_locf(sierra_ts,option = "locf",na_remaining = "rev") #metodo de ultima observacion hacia delante
sierra_imput_stl=na.interp(sierra_ts) #metodo de descomposici?n STL con interpolaci?n lineal


#graficar
plot.ts(sierra_ts,lwd=16,col="gray")
lines(sierra_imput_aa,lwd=2,col="red")
lines(sierra_imput_st,lwd=2,col="blue")
lines(sierra_imput_itl,lwd=2,col="orange")
lines(sierra_imput_its,lwd=2,col="purple")
lines(sierra_imput_itst,lwd=2,col="green")
lines(sierra_imput_mas,lwd=2,col="seashell3")
lines(sierra_imput_mal,lwd=2,col="salmon1")
lines(sierra_imput_mae,lwd=2,col="blue4")#CAMBIAR LA ESPECIE Y EL M?TODO EN >serie.temporal



#imput missing data to camaron_titi
statsNA(camaron_titi_ts) #24 NA, 23.8%
camaron_titi_imput_aa=na_kalman(camaron_titi_ts,model="auto.arima",smooth = TRUE,nit = -1) #m?todo Structured Model
camaron_titi_imput_st=na_kalman(camaron_titi_ts,model="StructTS",smooth = TRUE,nit = -1) #m?todo auto.arima
camaron_titi_imput_itl=na_interpolation(camaron_titi_ts,option="linear") #m?todo interpolacion lineal
camaron_titi_imput_its=na_interpolation(camaron_titi_ts,option="spline") #m?todo interpolaci?n spline
camaron_titi_imput_itst=na_interpolation(camaron_titi_ts,option="stine") #m?todo interpolaci?n spline
camaron_titi_imput_mas=na_ma(camaron_titi_ts,k=3,weighting = "simple") #m?todo moving average simple
camaron_titi_imput_mal=na_ma(camaron_titi_ts,k=3,weighting = "linear") #m?todo moving average linear
camaron_titi_imput_mae=na_ma(camaron_titi_ts,k=3,weighting = "exponential") #m?todo moving average exponential
camaron_titi_imput_seadec=na_seadec(camaron_titi_ts,algorithm = "kalman") #m?todo de imputacion kalman en serie estacionalmente descompuesta
camaron_titi_imput_locf=na_locf(camaron_titi_ts,option = "locf",na_remaining = "rev") #metodo de ultima observacion hacia delante


#graficar
plot.ts(camaron_titi_ts,lwd=16,col="gray")
lines(camaron_titi_imput_aa,lwd=2,col="red")
lines(camaron_titi_imput_st,lwd=2,col="blue")
lines(camaron_titi_imput_itl,lwd=2,col="orange")
lines(camaron_titi_imput_its,lwd=2,col="purple")
lines(camaron_titi_imput_itst,lwd=2,col="green")
lines(camaron_titi_imput_mas,lwd=2,col="seashell3")
lines(camaron_titi_imput_mal,lwd=2,col="salmon1")
lines(camaron_titi_imput_mae,lwd=2,col="blue4")



#imput missing data to yft
ggplot_na_distribution(yft_ts)
statsNA(yft_ts) #24 NA, 23.8%
yft.runs=ifelse(yft_ts >= 0, 1, 0);yft.runs[is.na(yft.runs)] <- 0 #convert time series to a binary vector where numbers = 1 and NA = 0
runs.test(yft.runs) # runs test to check if NAs appear at random
yft_imput_aa=na_kalman(yft_ts,model="auto.arima",smooth = TRUE,nit = -1) #m?todo Structured Model
yft_imput_st=na_kalman(yft_ts,model="StructTS",smooth = TRUE,nit = -1) #m?todo auto.arima
yft_imput_itl=na_interpolation(yft_ts,option="linear") #m?todo interpolacion lineal
yft_imput_its=na_interpolation(yft_ts,option="spline") #m?todo interpolaci?n spline
yft_imput_itst=na_interpolation(yft_ts,option="stine") #m?todo interpolaci?n spline
yft_imput_mas=na_ma(yft_ts,k=3,weighting = "simple") #m?todo moving average simple
yft_imput_mal=na_ma(yft_ts,k=3,weighting = "linear") #m?todo moving average linear
yft_imput_mae=na_ma(yft_ts,k=3,weighting = "exponential") #m?todo moving average exponential
yft_imput_seadec=na_seadec(yft_ts,algorithm = "kalman") #m?todo de imputacion kalman en serie estacionalmente descompuesta
yft_imput_locf=na_locf(yft_ts,option = "locf",na_remaining = "rev") #metodo de ultima observacion hacia delante
yft_imput_stl=na.interp(yft_ts) #metodo de descomposici?n STL con interpolaci?n lineal


#graficar
plot.ts(yft_ts,lwd=16,col="gray")
lines(yft_imput_aa,lwd=2,col="red")
lines(yft_imput_st,lwd=2,col="blue")
lines(yft_imput_itl,lwd=2,col="orange")
lines(yft_imput_its,lwd=2,col="purple")
lines(yft_imput_itst,lwd=2,col="green")
lines(yft_imput_mas,lwd=2,col="seashell3")
lines(yft_imput_mal,lwd=2,col="salmon1")
lines(yft_imput_mae,lwd=2,col="blue4")


#imput missing data to skj
statsNA(skj_ts) #24 NA, 23.8%
skj_imput_aa=na_kalman(skj_ts,model="auto.arima",smooth = TRUE,nit = -1) #m?todo Structured Model
skj_imput_st=na_kalman(skj_ts,model="StructTS",smooth = TRUE,nit = -1) #m?todo auto.arima
skj_imput_itl=na_interpolation(skj_ts,option="linear") #m?todo interpolacion lineal
skj_imput_its=na_interpolation(skj_ts,option="spline") #m?todo interpolaci?n spline
skj_imput_itst=na_interpolation(skj_ts,option="stine") #m?todo interpolaci?n spline
skj_imput_mas=na_ma(skj_ts,k=3,weighting = "simple") #m?todo moving average simple
skj_imput_mal=na_ma(skj_ts,k=3,weighting = "linear") #m?todo moving average linear
skj_imput_mae=na_ma(skj_ts,k=3,weighting = "exponential") #m?todo moving average exponential
skj_imput_seadec=na_seadec(skj_ts,algorithm = "kalman") #m?todo de imputacion kalman en serie estacionalmente descompuesta
skj_imput_locf=na_locf(skj_ts,option = "locf",na_remaining = "rev") #metodo de ultima observacion hacia delante


#graficar
plot.ts(skj_ts,lwd=16,col="gray")
lines(skj_imput_aa,lwd=2,col="red")
lines(skj_imput_st,lwd=2,col="blue")
lines(skj_imput_itl,lwd=2,col="orange")
lines(skj_imput_its,lwd=2,col="purple")
lines(skj_imput_itst,lwd=2,col="green")
lines(skj_imput_mas,lwd=2,col="seashell3")
lines(skj_imput_mal,lwd=2,col="salmon1")
lines(skj_imput_mae,lwd=2,col="blue4")








################################################################################
#################### 2. VALIDATION WHEN ASSUMING ONLY RANDOM NA'S###############
#Con los mismos datos pero los imputados
#Hacerlo para cada m?todo por separado
#HAcerlo para cada especie por separado
#cambiar especie y m?todo en <-serie.temporal


#para StructTS----------------------------------------------------------------

cantidad <- 5 #Cantidad de veces que voy a repetir mi loop | ciclo
cuantos_na <- seq(5,50,5) # Crear secuencia que representar? la cantidad
# de NA's con las que probar? los algoritmos
tabla_resultados_st <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_st) <- c('cantidad_nas', 'rmse_st', 'cor_st','mase_st','smape_st')
serie.temporal <- bagre_imput_st # CAMBIAR AQUI la serie con la que voy a trabajar

# Eliminar la serie de datos inventada
tabla_resultados_st <- tabla_resultados_st[-1,]

for (x in cuantos_na) {
  for (i in 1:cantidad) {
    # 1. Crear y - Con datos faltantes (NA) al azar
    y_na <- set_na(serie.temporal,na=sample(bagre_imput_st,x,replace = FALSE))

    # 2. Calcular los valores a imputar con los algoritmos de Modelos Estructurales y auto arima
    y_impute_st <- na_kalman(y_na,model = "StructTS",nit=-1)

    # 3. Calcular las metricas de desempe?o
    rmse_st <- rmse(serie.temporal,y_impute_st)
    cor_st <- cor(serie.temporal,y_impute_st)^2
    mase_st <- mase(serie.temporal,y_impute_st)
    smape_st <- smape(serie.temporal,y_impute_st)

    # 4. Almacenar en una tabla temporal las metricas
    temporal1 <- data.frame(x, rmse_st, cor_st,mase_st,smape_st)
    colnames(temporal1)[colnames(temporal1) == 'x'] <- 'cantidad_nas'

    # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
    tabla_resultados_st <- bind_rows(tabla_resultados_st, temporal1)
  }
}

#crear tabla en archivo txt
write.table(tabla_resultados_st,"tabla_resultados_st.txt")



#para Auto Arima----------------------------------------------------------------

cuantos_na <- seq(5,50,5) # Crear secuencia que representar? la cantidad
# de NA's con las que probar? los algoritmos
tabla_resultados_aa <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_aa) <- c('cantidad_nas', 'rmse_aa', 'cor_aa','mase_aa','smape_aa')
serie.temporal <- bagre_imput_aa # CAMBIAR AQUI la serie con la que voy a trabajar

# Eliminar la serie de datos inventada
tabla_resultados_aa <- tabla_resultados_aa[-1,]

for (x in cuantos_na) {
  for (i in 1:cantidad) {
    # 1. Crear y - Con datos faltantes (NA) al azar
    y_na <- set_na(serie.temporal,na=sample(bagre_imput_aa,x,replace = FALSE))

    # 2. Calcular los valores a imputar con el algoritmos auto.arima
    y_impute_aa <- na_kalman(y_na,model = "auto.arima",nit=-1)

    # 3. Calcular las metricas de desempe?o
    rmse_aa <- rmse(serie.temporal,y_impute_aa)
    cor_aa <- cor(serie.temporal,y_impute_aa)^2
    mase_aa <- mase(serie.temporal,y_impute_aa)
    smape_aa <- smape(serie.temporal,y_impute_aa)

    # 4. Almacenar en una tabla temporal las metricas
    temporal1 <- data.frame(x, rmse_aa, cor_aa,mase_aa,smape_aa)
    colnames(temporal1)[colnames(temporal1) == 'x'] <- 'cantidad_nas'

    # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
    tabla_resultados_aa <- bind_rows(tabla_resultados_aa, temporal1)
  }
}


#crear tabla en archivo txt
write.table(tabla_resultados_aa,"tabla_resultados_aa.txt")



#para Interpolacion Lineal----------------------------------------------------------------

cuantos_na <- seq(5,50,5) # Crear secuencia que representar? la cantidad
# de NA's con las que probar? los algoritmos
tabla_resultados_itl <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_itl) <- c('cantidad_nas', 'rmse_itl', 'cor_itl','mase_itl','smape_itl')
serie.temporal <- bagre_imput_itl # CAMBIAR AQUI la serie con la que voy a trabajar

# Eliminar la serie de datos inventada
tabla_resultados_itl <- tabla_resultados_itl[-1,]

for (x in cuantos_na) {
  for (i in 1:cantidad) {
    # 1. Crear y - Con datos faltantes (NA) al azar
    y_na <- set_na(serie.temporal,na=sample(bagre_imput_itl,x,replace = FALSE))

    # 2. Calcular los valores a imputar el algoritmo interpolacion linear
    y_impute_itl <- na_interpolation(y_na,option="linear")

    # 3. Calcular las metricas de desempe?o
    rmse_itl <- rmse(serie.temporal,y_impute_itl)
    cor_itl <- cor(serie.temporal,y_impute_itl)^2
    mase_itl <- mase(serie.temporal,y_impute_itl)
    smape_itl <- smape(serie.temporal,y_impute_itl)

    # 4. Almacenar en una tabla temporal las metricas
    temporal1 <- data.frame(x, rmse_itl, cor_itl,mase_itl,smape_itl)
    colnames(temporal1)[colnames(temporal1) == 'x'] <- 'cantidad_nas'

    # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
    tabla_resultados_itl <- bind_rows(tabla_resultados_itl, temporal1)
  }
}

#crear tabla en archivo txt
write.table(tabla_resultados_itl,"tabla_resultados_itl.txt")



#para Interpolacion Spline----------------------------------------------------------------

cuantos_na <- seq(5,50,5) # Crear secuencia que representar? la cantidad
# de NA's con las que probar? los algoritmos
tabla_resultados_its <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_its) <- c('cantidad_nas', 'rmse_its', 'cor_its','mase_its','smape_its')
serie.temporal <- bagre_imput_its # CAMBIAR AQUI la serie con la que voy a trabajar

# Eliminar la serie de datos inventada
tabla_resultados_its <- tabla_resultados_its[-1,]

for (x in cuantos_na) {
  for (i in 1:cantidad) {
    # 1. Crear y - Con datos faltantes (NA) al azar
    y_na <- set_na(serie.temporal,na=sample(bagre_imput_its,x,replace = FALSE))

    # 2. Calcular los valores a imputar con interpolacion spline
    y_impute_its <- na_interpolation(y_na,option="spline")

    # 3. Calcular las metricas de desempe?o
    rmse_its <- rmse(serie.temporal,y_impute_its)
    cor_its <- cor(serie.temporal,y_impute_its)^2
    mase_its <- mase(serie.temporal,y_impute_its)
    smape_its <- smape(serie.temporal,y_impute_its)

    # 4. Almacenar en una tabla temporal las metricas
    temporal1 <- data.frame(x, rmse_its, cor_its,mase_its,smape_its)
    colnames(temporal1)[colnames(temporal1) == 'x'] <- 'cantidad_nas'

    # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
    tabla_resultados_its <- bind_rows(tabla_resultados_its, temporal1)
  }
}

#crear tabla en archivo txt
write.table(tabla_resultados_its,"tabla_resultados_its.txt")



#para Interpolacion Stineman----------------------------------------------------------------

cuantos_na <- seq(5,50,5) # Crear secuencia que representar? la cantidad
# de NA's con las que probar? los algoritmos
tabla_resultados_itst <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_itst) <- c('cantidad_nas', 'rmse_itst', 'cor_itst','mase_itst','smape_itst')
serie.temporal <- bagre_imput_itst # CAMBIAR AQUI la serie con la que voy a trabajar

# Eliminar la serie de datos inventada
tabla_resultados_itst <- tabla_resultados_itst[-1,]

for (x in cuantos_na) {
  for (i in 1:cantidad) {
    # 1. Crear y - Con datos faltantes (NA) al azar
    y_na <- set_na(serie.temporal,na=sample(bagre_imput_itst,x,replace = FALSE))

    # 2. Calcular los valores a imputar con interpolacion stineman
    y_impute_itst <- na_interpolation(y_na,option="stine")

    # 3. Calcular las metricas de desempe?o
    rmse_itst <- rmse(serie.temporal,y_impute_itst)
    cor_itst <- cor(serie.temporal,y_impute_itst)^2
    mase_itst <- mase(serie.temporal,y_impute_itst)
    smape_itst <- smape(serie.temporal,y_impute_itst)

    # 4. Almacenar en una tabla temporal las metricas
    temporal1 <- data.frame(x, rmse_itst, cor_itst,mase_itst,smape_itst)
    colnames(temporal1)[colnames(temporal1) == 'x'] <- 'cantidad_nas'

    # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
    tabla_resultados_itst <- bind_rows(tabla_resultados_itst, temporal1)
  }
}

#crear tabla en archivo txt
write.table(tabla_resultados_itst,"tabla_resultados_itst.txt")





#para Moving Agerage Simple----------------------------------------------------------------

cuantos_na <- seq(5,50,5) # Crear secuencia que representar? la cantidad
# de NA's con las que probar? los algoritmos
tabla_resultados_mas <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_mas) <- c('cantidad_nas', 'rmse_mas', 'cor_mas','mase_mas','smape_mas')
serie.temporal <- bagre_imput_mas # CAMBIAR AQUI la serie con la que voy a trabajar

# Eliminar la serie de datos inventada
tabla_resultados_mas <- tabla_resultados_mas[-1,]

for (x in cuantos_na) {
  for (i in 1:cantidad) {
    # 1. Crear y - Con datos faltantes (NA) al azar
    y_na <- set_na(serie.temporal,na=sample(bagre_imput_mas,x,replace = FALSE))

    # 2. Calcular los valores a imputar con moving average simple
    y_impute_mas <- na_ma(y_na,k=4,weighting = "simple")

    # 3. Calcular las metricas de desempe?o
    rmse_mas <- rmse(serie.temporal,y_impute_mas)
    cor_mas <- cor(serie.temporal,y_impute_mas)^2
    mase_mas <- mase(serie.temporal,y_impute_mas)
    smape_mas <- smape(serie.temporal,y_impute_mas)

    # 4. Almacenar en una tabla temporal las metricas
    temporal1 <- data.frame(x, rmse_mas, cor_mas,mase_mas,smape_mas)
    colnames(temporal1)[colnames(temporal1) == 'x'] <- 'cantidad_nas'

    # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
    tabla_resultados_mas <- bind_rows(tabla_resultados_mas, temporal1)
  }
}

#crear tabla en archivo txt
write.table(tabla_resultados_mas,"tabla_resultados_mas.txt")



#para Moving Agerage Linear----------------------------------------------------------------

cuantos_na <- seq(5,50,5) # Crear secuencia que representar? la cantidad
# de NA's con las que probar? los algoritmos
tabla_resultados_mal <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_mal) <- c('cantidad_nas', 'rmse_mal', 'cor_mal','mase_mal','smape_mal')
serie.temporal <- bagre_imput_mal # CAMBIAR AQUI la serie con la que voy a trabajar

# Eliminar la serie de datos inventada
tabla_resultados_mal <- tabla_resultados_mal[-1,]

for (x in cuantos_na) {
  for (i in 1:cantidad) {
    # 1. Crear y - Con datos faltantes (NA) al azar
    y_na <- set_na(serie.temporal,na=sample(bagre_imput_mal,x,replace = FALSE))

    # 2. Calcular los valores a imputar con moving average simple
    y_impute_mal <- na_ma(y_na,k=4,weighting = "linear")

    # 3. Calcular las metricas de desempe?o
    rmse_mal <- rmse(serie.temporal,y_impute_mal)
    cor_mal <- cor(serie.temporal,y_impute_mal)^2
    mase_mal <- mase(serie.temporal,y_impute_mal)
    smape_mal <- smape(serie.temporal,y_impute_mal)

    # 4. Almacenar en una tabla temporal las metricas
    temporal1 <- data.frame(x, rmse_mal, cor_mal,mase_mal,smape_mal)
    colnames(temporal1)[colnames(temporal1) == 'x'] <- 'cantidad_nas'

    # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
    tabla_resultados_mal <- bind_rows(tabla_resultados_mal, temporal1)
  }
}

#crear tabla en archivo txt
write.table(tabla_resultados_mal,"tabla_resultados_mal.txt")




#para Moving Agerage Exponential----------------------------------------------------------------

cuantos_na <- seq(5,50,5) # Crear secuencia que representar? la cantidad
# de NA's con las que probar? los algoritmos
tabla_resultados_mae <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_mae) <- c('cantidad_nas', 'rmse_mae', 'cor_mae','mase_mae','smape_mae')
serie.temporal <- bagre_imput_mae # CAMBIAR AQUI la serie con la que voy a trabajar

# Eliminar la serie de datos inventada
tabla_resultados_mae <- tabla_resultados_mae[-1,]

for (x in cuantos_na) {
  for (i in 1:cantidad) {
    # 1. Crear y - Con datos faltantes (NA) al azar
    y_na <- set_na(serie.temporal,na=sample(bagre_imput_mae,x,replace = FALSE))

    # 2. Calcular los valores a imputar con moving average exponential
    y_impute_mae <- na_ma(y_na,k=4,weighting = "exponential")

    # 3. Calcular las metricas de desempe?o
    rmse_mae <- rmse(serie.temporal,y_impute_mae)
    cor_mae <- cor(serie.temporal,y_impute_mae)^2
    mase_mae <- mase(serie.temporal,y_impute_mae)
    smape_mae <- smape(serie.temporal,y_impute_mae)

    # 4. Almacenar en una tabla temporal las metricas
    temporal1 <- data.frame(x, rmse_mae, cor_mae,mase_mae,smape_mae)
    colnames(temporal1)[colnames(temporal1) == 'x'] <- 'cantidad_nas'

    # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
    tabla_resultados_mae <- bind_rows(tabla_resultados_mae, temporal1)
  }
}

#crear tabla en archivo txt
write.table(tabla_resultados_mae,"tabla_resultados_mae.txt")




#para SEADEC----------------------------------------------------------------

cuantos_na <- seq(5,50,5) # Crear secuencia que representar? la cantidad
# de NA's con las que probar? los algoritmos
tabla_resultados_seadec <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_seadec) <- c('cantidad_nas', 'rmse_seadec', 'cor_seadec','mase_seadec','smape_seadec')
serie.temporal <- bagre_imput_seadec # CAMBIAR AQUI la serie con la que voy a trabajar

# Eliminar la serie de datos inventada
tabla_resultados_seadec <- tabla_resultados_seadec[-1,]

for (x in cuantos_na) {
  for (i in 1:cantidad) {
    # 1. Crear y - Con datos faltantes (NA) al azar en la serie original con los NA sacados
    y_na <- set_na(serie.temporal,na=sample(bagre_imput_seadec,x,replace = FALSE))

    # 2. Calcular los valores a imputar con moving average exponential
    y_impute_seadec <- na_seadec(y_na,algorithm = "kalman")

    # 3. Calcular las metricas de desempe?o
    rmse_seadec <- rmse(serie.temporal,y_impute_seadec)
    cor_seadec <- cor(serie.temporal,y_impute_seadec)^2
    mase_seadec <- mase(serie.temporal,y_impute_seadec)
    smape_seadec <- smape(serie.temporal,y_impute_seadec)

    # 4. Almacenar en una tabla temporal las metricas
    temporal1 <- data.frame(x, rmse_seadec, cor_seadec,mase_seadec,smape_seadec)
    colnames(temporal1)[colnames(temporal1) == 'x'] <- 'cantidad_nas'

    # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
    tabla_resultados_seadec <- bind_rows(tabla_resultados_seadec, temporal1)
  }
}

#crear tabla en archivo txt
write.table(tabla_resultados_seadec,"tabla_resultados_seadec.txt")


#para LOCF----------------------------------------------------------------

cuantos_na <- seq(5,50,5) # Crear secuencia que representar? la cantidad
# de NA's con las que probar? los algoritmos
tabla_resultados_locf <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_locf) <- c('cantidad_nas', 'rmse_locf', 'cor_locf','mase_locf','smape_locf')
serie.temporal <- bagre_imput_locf # CAMBIAR AQUI la serie con la que voy a trabajar

# Eliminar la serie de datos inventada
tabla_resultados_locf <- tabla_resultados_locf[-1,]

for (x in cuantos_na) {
  for (i in 1:cantidad) {
    # 1. Crear y - Con datos faltantes (NA) al azar en la serie original con los NA sacados
    y_na <- set_na(serie.temporal,na=sample(bagre_imput_locf,x,replace = FALSE))

    # 2. Calcular los valores a imputar con moving average exponential
    y_impute_locf <- na_locf(y_na,option = "locf", na_remaining = "rev")

    # 3. Calcular las metricas de desempe?o
    rmse_locf <- rmse(serie.temporal,y_impute_locf)
    cor_locf <- cor(serie.temporal,y_impute_locf)^2
    mase_locf <- mase(serie.temporal,y_impute_locf)
    smape_locf <- smape(serie.temporal,y_impute_locf)

    # 4. Almacenar en una tabla temporal las metricas
    temporal1 <- data.frame(x, rmse_locf, cor_locf,mase_locf,smape_locf)
    colnames(temporal1)[colnames(temporal1) == 'x'] <- 'cantidad_nas'

    # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
    tabla_resultados_locf <- bind_rows(tabla_resultados_locf, temporal1)
  }
}

#crear tabla en archivo txt
write.table(tabla_resultados_locf,"tabla_resultados_locf.txt")






#para Interpolation con STL decomposition----------------------------------------------------------------

cuantos_na <- seq(5,50,5) # Crear secuencia que representar? la cantidad
# de NA's con las que probar? los algoritmos
tabla_resultados_stl <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_stl) <- c('cantidad_nas', 'rmse_stl', 'cor_stl','mase_stl','smape_stl')
serie.temporal <- bagre_imput_stl # CAMBIAR AQUI la serie con la que voy a trabajar

# Eliminar la serie de datos inventada
tabla_resultados_stl <- tabla_resultados_stl[-1,]

for (x in cuantos_na) {
  for (i in 1:cantidad) {
    # 1. Crear y - Con datos faltantes (NA) al azar en la serie original con los NA sacados
    y_na <- set_na(serie.temporal,na=sample(bagre_imput_stl,x,replace = FALSE))

    # 2. Calcular los valores a imputar con moving average exponential
    y_impute_stl <- na.interp(y_na)

    # 3. Calcular las metricas de desempe?o
    rmse_stl <- rmse(serie.temporal,y_impute_stl)
    cor_stl <- cor(serie.temporal,y_impute_stl)^2
    mase_stl <- mase(serie.temporal,y_impute_stl)
    smape_stl <- smape(serie.temporal,y_impute_stl)

    # 4. Almacenar en una tabla temporal las metricas
    temporal1 <- data.frame(x, rmse_stl, cor_stl,mase_stl,smape_stl)
    colnames(temporal1)[colnames(temporal1) == 'x'] <- 'cantidad_nas'

    # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
    tabla_resultados_stl <- bind_rows(tabla_resultados_stl, temporal1)
  }
}

#crear tabla en archivo txt
write.table(tabla_resultados_stl,"tabla_resultados_stl.txt")



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#Load data and create a whole table for each species

val_st=read.table("tabla_resultados_st.txt",h=T)
val_aa=read.table("tabla_resultados_aa.txt",h=T)
val_itl=read.table("tabla_resultados_itl.txt",h=T)
val_its=read.table("tabla_resultados_its.txt",h=T)
val_itst=read.table("tabla_resultados_itst.txt",h=T)
val_mas=read.table("tabla_resultados_mas.txt",h=T)
val_mal=read.table("tabla_resultados_mal.txt",h=T)
val_mae=read.table("tabla_resultados_mae.txt",h=T)
val_seadec=read.table("tabla_resultados_seadec.txt",h=T)
val_locf=read.table("tabla_resultados_locf.txt",h=T)
val_stl=read.table("tabla_resultados_stl.txt",h=T)


tabla=data.frame(val_aa$cantidad_nas,val_aa$rmse_aa,
                 val_aa$cor_aa,val_aa$mase_aa,val_aa$smape_aa,
                 val_itl$rmse_itl,val_itl$cor_itl,val_itl$mase_itl,
                 val_itl$smape_itl,val_its$rmse_its,val_its$cor_its,
                 val_its$mase_its,val_its$smape_its,val_itst$rmse_itst,
                 val_itst$cor_itst,val_itst$mase_itst,val_itst$smape_itst,
                 val_mae$rmse_mae,val_mae$cor_mae,val_mae$mase_mae,val_mae$smape_mae,
                 val_mal$rmse_mal,val_mal$cor_mal,val_mal$mase_mal,
                 val_mal$smape_mal,val_mas$rmse_mas,val_mas$cor_mas,
                 val_mas$mase_mas,val_mas$smape_mas,val_st$rmse_st,
                 val_st$cor_st,val_st$mase_st,val_st$smape_st,val_seadec$rmse_seadec,
                 val_seadec$cor_seadec,val_seadec$mase_seadec,val_seadec$smape_seadec,
                 val_locf$rmse_locf,val_locf$cor_locf,val_locf$mase_locf,val_locf$smape_locf,
                 val_stl$rmse_stl,val_stl$cor_stl,val_stl$mase_stl,val_stl$smape_stl)

colnames(tabla)=c("nas","aa_rmse","aa_cor","aa_mase","aa_smape","itl_rmse",
                  "itl_cor","itl_mase","itl_smape","its_rmse","its_cor","its_mase",
                  "its_smape","itst_rmse","itst_cor","itst_mase","itst_smape","mae_rmse","mae_cor",
                  "mae_mase","mae_smape","mal_rmse","mal_cor","mal_mase",
                  "mal_smape","mas_rmse","mas_cor","mas_mase","mas_smape",
                  "st_rmse","st_cor","st_mase","st_smape","seadec_rmse","seadec_cor","seadec_mase","seadec_smape",
                  "locf_rmse","locf_cor","locf_mase","locf_smape",
                  "stl_rmse","stl_cor","stl_mase","stl_smape")


write.table(tabla,"tabla_final.txt")

###########################################





################################################################################
#################### 2. VALIDATION WHEN SIMULATING THE NA'S STRUCTURE OF DATA IN CHUNKSSS ###############
#Con los mismos datos pero los imputados
#Hacerlo para cada m?todo por separado
#HAcerlo para cada especie por separado
#cambiar especie y m?todo en <-serie.temporal



#cantidad de veces a repetir el loop (same for all)
cantidad=1000
#espaciado en la serie donde van los gaps (same for all)
gap1=seq(1,12,1);gap2=seq(13,25,1);gap3=seq(26,38,1);gap4=seq(39,51,1);gap5=seq(52,64,1);gap6=seq(65,77,1);gap7=seq(78,90,1);gap8=seq(91,101,1)




#### PARA STRUCTS

tabla_resultados_st <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_st) <- c('nas_st','rmse_st', 'cor_st','mase_st','smape_st')
serie.temporal <- pelada_imput_st  # CAMBIAR AQUI la serie con la que voy a trabajar
# Eliminar la serie de datos inventada
tabla_resultados_st <- tabla_resultados_st[-1,]

for (i in 1:cantidad) {
  # 1. Crear y - Con datos faltantes (NA) al azar en la serie original con los NA sacados
  y_na <- set_na(serie.temporal, na=  c(serie.temporal[sample(gap1,1):sample(gap1,1)],serie.temporal[sample(gap2,1):sample(gap2,1)],serie.temporal[sample(gap3,1):sample(gap3,1)],serie.temporal[sample(gap4,1):sample(gap4,1)],serie.temporal[sample(gap5,1):sample(gap5,1)],serie.temporal[sample(gap6,1):sample(gap6,1)],serie.temporal[sample(gap7,1):sample(gap7,1)],serie.temporal[sample(gap8,1):sample(gap8,1)]))

  # 2. Calcular los valores a imputar con moving average exponential
  y_impute_st <- na_kalman(y_na,model="StructTS",nit = -1)

  # 3. Calcular las metricas de desempe?o
  nas_st <-sum(is.na(y_na))
  rmse_st <- rmse(serie.temporal,y_impute_st)
  cor_st <- cor(serie.temporal,y_impute_st)^2
  mase_st <- mase(serie.temporal,y_impute_st)
  smape_st <- smape(serie.temporal,y_impute_st)


  # 4. Almacenar en una tabla temporal las metricas
  temporal1 <- data.frame(nas_st,rmse_st, cor_st,mase_st,smape_st)
  colnames(temporal1)=c("nas_st","rmse_st","cor_st","mase_st","smape_st")

  # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
  tabla_resultados_st <- bind_rows(tabla_resultados_st, temporal1)
}

#crear tabla en archivo txt
write.table(tabla_resultados_st,"tabla_resultados_st2.txt")







#### PARA AUTO.ARIMA


tabla_resultados_aa <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_aa) <- c('nas_aa','rmse_aa', 'cor_aa','mase_aa','smape_aa')
serie.temporal <- pelada_imput_aa # CAMBIAR AQUI la serie con la que voy a trabajar
# Eliminar la serie de datos inventada
tabla_resultados_aa <- tabla_resultados_aa[-1,]

for (i in 1:cantidad) {
  # 1. Crear y - Con datos faltantes (NA) al azar en la serie original con los NA sacados
  y_na <- set_na(serie.temporal, na=  c(serie.temporal[sample(gap1,1):sample(gap1,1)],serie.temporal[sample(gap2,1):sample(gap2,1)],serie.temporal[sample(gap3,1):sample(gap3,1)],serie.temporal[sample(gap4,1):sample(gap4,1)],serie.temporal[sample(gap5,1):sample(gap5,1)],serie.temporal[sample(gap6,1):sample(gap6,1)],serie.temporal[sample(gap7,1):sample(gap7,1)],serie.temporal[sample(gap8,1):sample(gap8,1)]))

  # 2. Calcular los valores a imputar con moving average exponential
  y_impute_aa <- na_kalman(y_na,model="auto.arima",nit = -1)

  # 3. Calcular las metricas de desempe?o
  nas_aa <-sum(is.na(y_na))
  rmse_aa <- rmse(serie.temporal,y_impute_aa)
  cor_aa <- cor(serie.temporal,y_impute_aa)^2
  mase_aa <- mase(serie.temporal,y_impute_aa)
  smape_aa <- smape(serie.temporal,y_impute_aa)


  # 4. Almacenar en una tabla temporal las metricas
  temporal1 <- data.frame(nas_aa,rmse_aa, cor_aa,mase_aa,smape_aa)
  colnames(temporal1)=c("nas_aa","rmse_aa","cor_aa","mase_aa","smape_aa")

  # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
  tabla_resultados_aa <- bind_rows(tabla_resultados_aa, temporal1)
}

#crear tabla en archivo txt
write.table(tabla_resultados_aa,"tabla_resultados_aa2.txt")







#### PARA LINEAR INTERPOLATION


tabla_resultados_itl <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_itl) <- c('nas_itl','rmse_itl', 'cor_itl','mase_itl','smape_itl')
serie.temporal <- pelada_imput_itl # CAMBIAR AQUI la serie con la que voy a trabajar
# Eliminar la serie de datos inventada
tabla_resultados_itl <- tabla_resultados_itl[-1,]

for (i in 1:cantidad) {
  # 1. Crear y - Con datos faltantes (NA) al azar en la serie original con los NA sacados
  y_na <- set_na(serie.temporal, na=  c(serie.temporal[sample(gap1,1):sample(gap1,1)],serie.temporal[sample(gap2,1):sample(gap2,1)],serie.temporal[sample(gap3,1):sample(gap3,1)],serie.temporal[sample(gap4,1):sample(gap4,1)],serie.temporal[sample(gap5,1):sample(gap5,1)],serie.temporal[sample(gap6,1):sample(gap6,1)],serie.temporal[sample(gap7,1):sample(gap7,1)],serie.temporal[sample(gap8,1):sample(gap8,1)]))

  # 2. Calcular los valores a imputar con moving average exponential
  y_impute_itl <- na_interpolation(y_na,option="linear")

  # 3. Calcular las metricas de desempe?o
  nas_itl <-sum(is.na(y_na))
  rmse_itl <- rmse(serie.temporal,y_impute_itl)
  cor_itl <- cor(serie.temporal,y_impute_itl)^2
  mase_itl <- mase(serie.temporal,y_impute_itl)
  smape_itl <- smape(serie.temporal,y_impute_itl)


  # 4. Almacenar en una tabla temporal las metricas
  temporal1 <- data.frame(nas_itl,rmse_itl, cor_itl,mase_itl,smape_itl)
  colnames(temporal1)=c("nas_itl","rmse_itl","cor_itl","mase_itl","smape_itl")

  # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
  tabla_resultados_itl <- bind_rows(tabla_resultados_itl, temporal1)
}

#crear tabla en archivo txt
write.table(tabla_resultados_itl,"tabla_resultados_itl2.txt")




#### PARA INTERPOLATION SPLINE


tabla_resultados_its <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_its) <- c('nas_its','rmse_its', 'cor_its','mase_its','smape_its')
serie.temporal <- pelada_imput_its # CAMBIAR AQUI la serie con la que voy a trabajar
# Eliminar la serie de datos inventada
tabla_resultados_its <- tabla_resultados_its[-1,]

for (i in 1:cantidad) {
  # 1. Crear y - Con datos faltantes (NA) al azar en la serie original con los NA sacados
  y_na <- set_na(serie.temporal, na=  c(serie.temporal[sample(gap1,1):sample(gap1,1)],serie.temporal[sample(gap2,1):sample(gap2,1)],serie.temporal[sample(gap3,1):sample(gap3,1)],serie.temporal[sample(gap4,1):sample(gap4,1)],serie.temporal[sample(gap5,1):sample(gap5,1)],serie.temporal[sample(gap6,1):sample(gap6,1)],serie.temporal[sample(gap7,1):sample(gap7,1)],serie.temporal[sample(gap8,1):sample(gap8,1)]))

  # 2. Calcular los valores a imputar con moving average exponential
  y_impute_its <- na_interpolation(y_na,option="spline")

  # 3. Calcular las metricas de desempe?o
  nas_its <-sum(is.na(y_na))
  rmse_its <- rmse(serie.temporal,y_impute_its)
  cor_its <- cor(serie.temporal,y_impute_its)^2
  mase_its <- mase(serie.temporal,y_impute_its)
  smape_its <- smape(serie.temporal,y_impute_its)


  # 4. Almacenar en una tabla temporal las metricas
  temporal1 <- data.frame(nas_its,rmse_its, cor_its,mase_its,smape_its)
  colnames(temporal1)=c("nas_its","rmse_its","cor_its","mase_its","smape_its")

  # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
  tabla_resultados_its <- bind_rows(tabla_resultados_its, temporal1)
}

#crear tabla en archivo txt
write.table(tabla_resultados_its,"tabla_resultados_its2.txt")






#### PARA INTERPOLATION STINEMAN


tabla_resultados_itst <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_itst) <- c('nas_itst','rmse_itst', 'cor_itst','mase_itst','smape_itst')
serie.temporal <- pelada_imput_itst # CAMBIAR AQUI la serie con la que voy a trabajar
# Eliminar la serie de datos inventada
tabla_resultados_itst <- tabla_resultados_itst[-1,]

for (i in 1:cantidad) {
  # 1. Crear y - Con datos faltantes (NA) al azar en la serie original con los NA sacados
  y_na <- set_na(serie.temporal, na=  c(serie.temporal[sample(gap1,1):sample(gap1,1)],serie.temporal[sample(gap2,1):sample(gap2,1)],serie.temporal[sample(gap3,1):sample(gap3,1)],serie.temporal[sample(gap4,1):sample(gap4,1)],serie.temporal[sample(gap5,1):sample(gap5,1)],serie.temporal[sample(gap6,1):sample(gap6,1)],serie.temporal[sample(gap7,1):sample(gap7,1)],serie.temporal[sample(gap8,1):sample(gap8,1)]))

  # 2. Calcular los valores a imputar con moving average exponential
  y_impute_itst <- na_interpolation(y_na,option="stine")

  # 3. Calcular las metricas de desempe?o
  nas_itst <-sum(is.na(y_na))
  rmse_itst <- rmse(serie.temporal,y_impute_itst)
  cor_itst <- cor(serie.temporal,y_impute_itst)^2
  mase_itst <- mase(serie.temporal,y_impute_itst)
  smape_itst <- smape(serie.temporal,y_impute_itst)


  # 4. Almacenar en una tabla temporal las metricas
  temporal1 <- data.frame(nas_itst,rmse_itst, cor_itst,mase_itst,smape_itst)
  colnames(temporal1)=c("nas_itst","rmse_itst","cor_itst","mase_itst","smape_itst")

  # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
  tabla_resultados_itst <- bind_rows(tabla_resultados_itst, temporal1)
}

#crear tabla en archivo txt
write.table(tabla_resultados_itst,"tabla_resultados_itst2.txt")




#### PARA MOVING AVERAGE SIMPLE


tabla_resultados_mas <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_mas) <- c('nas_mas','rmse_mas', 'cor_mas','mase_mas','smape_mas')
serie.temporal <- pelada_imput_mas # CAMBIAR AQUI la serie con la que voy a trabajar
# Eliminar la serie de datos inventada
tabla_resultados_mas <- tabla_resultados_mas[-1,]

for (i in 1:cantidad) {
  # 1. Crear y - Con datos faltantes (NA) al azar en la serie original con los NA sacados
  y_na <- set_na(serie.temporal, na=  c(serie.temporal[sample(gap1,1):sample(gap1,1)],serie.temporal[sample(gap2,1):sample(gap2,1)],serie.temporal[sample(gap3,1):sample(gap3,1)],serie.temporal[sample(gap4,1):sample(gap4,1)],serie.temporal[sample(gap5,1):sample(gap5,1)],serie.temporal[sample(gap6,1):sample(gap6,1)],serie.temporal[sample(gap7,1):sample(gap7,1)],serie.temporal[sample(gap8,1):sample(gap8,1)]))

  # 2. Calcular los valores a imputar con moving average exponential
  y_impute_mas <- na_ma(y_na, weighting = "simple")

  # 3. Calcular las metricas de desempe?o
  nas_mas <-sum(is.na(y_na))
  rmse_mas <- rmse(serie.temporal,y_impute_mas)
  cor_mas <- cor(serie.temporal,y_impute_mas)^2
  mase_mas <- mase(serie.temporal,y_impute_mas)
  smape_mas <- smape(serie.temporal,y_impute_mas)


  # 4. Almacenar en una tabla temporal las metricas
  temporal1 <- data.frame(nas_mas,rmse_mas, cor_mas,mase_mas,smape_mas)
  colnames(temporal1)=c("nas_mas","rmse_mas","cor_mas","mase_mas","smape_mas")

  # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
  tabla_resultados_mas <- bind_rows(tabla_resultados_mas, temporal1)
}

#crear tabla en archivo txt
write.table(tabla_resultados_mas,"tabla_resultados_mas2.txt")




#### PARA MOVING AVERAGE LINEAR


tabla_resultados_mal <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_mal) <- c('nas_mal','rmse_mal', 'cor_mal','mase_mal','smape_mal')
serie.temporal <- pelada_imput_mal # CAMBIAR AQUI la serie con la que voy a trabajar
# Eliminar la serie de datos inventada
tabla_resultados_mal <- tabla_resultados_mal[-1,]

for (i in 1:cantidad) {
  # 1. Crear y - Con datos faltantes (NA) al azar en la serie original con los NA sacados
  y_na <- set_na(serie.temporal, na=  c(serie.temporal[sample(gap1,1):sample(gap1,1)],serie.temporal[sample(gap2,1):sample(gap2,1)],serie.temporal[sample(gap3,1):sample(gap3,1)],serie.temporal[sample(gap4,1):sample(gap4,1)],serie.temporal[sample(gap5,1):sample(gap5,1)],serie.temporal[sample(gap6,1):sample(gap6,1)],serie.temporal[sample(gap7,1):sample(gap7,1)],serie.temporal[sample(gap8,1):sample(gap8,1)]))

  # 2. Calcular los valores a imputar con moving average exponential
  y_impute_mal <- na_ma(y_na,weighting = "linear")

  # 3. Calcular las metricas de desempe?o
  nas_mal <-sum(is.na(y_na))
  rmse_mal <- rmse(serie.temporal,y_impute_mal)
  cor_mal <- cor(serie.temporal,y_impute_mal)^2
  mase_mal <- mase(serie.temporal,y_impute_mal)
  smape_mal <- smape(serie.temporal,y_impute_mal)


  # 4. Almacenar en una tabla temporal las metricas
  temporal1 <- data.frame(nas_mal,rmse_mal, cor_mal,mase_mal,smape_mal)
  colnames(temporal1)=c("nas_mal","rmse_mal","cor_mal","mase_mal","smape_mal")

  # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
  tabla_resultados_mal <- bind_rows(tabla_resultados_mal, temporal1)
}

#crear tabla en archivo txt
write.table(tabla_resultados_mal,"tabla_resultados_mal2.txt")








#### PARA MOVING AVERAGE EXPONENTIAL


tabla_resultados_mae <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_mae) <- c('nas_mae','rmse_mae', 'cor_mae','mase_mae','smape_mae')
serie.temporal <- pelada_imput_mae # CAMBIAR AQUI la serie con la que voy a trabajar
# Eliminar la serie de datos inventada
tabla_resultados_mae <- tabla_resultados_mae[-1,]

for (i in 1:cantidad) {
  # 1. Crear y - Con datos faltantes (NA) al azar en la serie original con los NA sacados
  y_na <- set_na(serie.temporal, na=  c(serie.temporal[sample(gap1,1):sample(gap1,1)],serie.temporal[sample(gap2,1):sample(gap2,1)],serie.temporal[sample(gap3,1):sample(gap3,1)],serie.temporal[sample(gap4,1):sample(gap4,1)],serie.temporal[sample(gap5,1):sample(gap5,1)],serie.temporal[sample(gap6,1):sample(gap6,1)],serie.temporal[sample(gap7,1):sample(gap7,1)],serie.temporal[sample(gap8,1):sample(gap8,1)]))

  # 2. Calcular los valores a imputar con moving average exponential
  y_impute_mae <- na_ma(y_na,weighting = "exponential")

  # 3. Calcular las metricas de desempe?o
  nas_mae <-sum(is.na(y_na))
  rmse_mae <- rmse(serie.temporal,y_impute_mae)
  cor_mae <- cor(serie.temporal,y_impute_mae)^2
  mase_mae <- mase(serie.temporal,y_impute_mae)
  smape_mae <- smape(serie.temporal,y_impute_mae)


  # 4. Almacenar en una tabla temporal las metricas
  temporal1 <- data.frame(nas_mae,rmse_mae, cor_mae,mase_mae,smape_mae)
  colnames(temporal1)=c("nas_mae","rmse_mae","cor_mae","mase_mae","smape_mae")

  # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
  tabla_resultados_mae <- bind_rows(tabla_resultados_mae, temporal1)
}

#crear tabla en archivo txt
write.table(tabla_resultados_mae,"tabla_resultados_mae2.txt")





#### PARA SEADEC


tabla_resultados_seadec <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_seadec) <- c('nas_seadec','rmse_seadec', 'cor_seadec','mase_seadec','smape_seadec')
serie.temporal <- pelada_imput_seadec # CAMBIAR AQUI la serie con la que voy a trabajar
# Eliminar la serie de datos inventada
tabla_resultados_seadec <- tabla_resultados_seadec[-1,]

for (i in 1:cantidad) {
  # 1. Crear y - Con datos faltantes (NA) al azar en la serie original con los NA sacados
  y_na <- set_na(serie.temporal, na=  c(serie.temporal[sample(gap1,1):sample(gap1,1)],serie.temporal[sample(gap2,1):sample(gap2,1)],serie.temporal[sample(gap3,1):sample(gap3,1)],serie.temporal[sample(gap4,1):sample(gap4,1)],serie.temporal[sample(gap5,1):sample(gap5,1)],serie.temporal[sample(gap6,1):sample(gap6,1)],serie.temporal[sample(gap7,1):sample(gap7,1)],serie.temporal[sample(gap8,1):sample(gap8,1)]))

  # 2. Calcular los valores a imputar con moving average exponential
  y_impute_seadec <- na_seadec(y_na,algorithm = "kalman")

  # 3. Calcular las metricas de desempe?o
  nas_seadec <-sum(is.na(y_na))
  rmse_seadec <- rmse(serie.temporal,y_impute_seadec)
  cor_seadec <- cor(serie.temporal,y_impute_seadec)^2
  mase_seadec <- mase(serie.temporal,y_impute_seadec)
  smape_seadec <- smape(serie.temporal,y_impute_seadec)


  # 4. Almacenar en una tabla temporal las metricas
  temporal1 <- data.frame(nas_seadec,rmse_seadec, cor_seadec,mase_seadec,smape_seadec)
  colnames(temporal1)=c("nas_seadec","rmse_seadec","cor_seadec","mase_seadec","smape_seadec")

  # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
  tabla_resultados_seadec <- bind_rows(tabla_resultados_seadec, temporal1)
}

#crear tabla en archivo txt
write.table(tabla_resultados_seadec,"tabla_resultados_seadec2.txt")




#### PARA LOCF


tabla_resultados_locf <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_locf) <- c('nas_locf','rmse_locf', 'cor_locf','mase_locf','smape_locf')
serie.temporal <- pelada_imput_locf # CAMBIAR AQUI la serie con la que voy a trabajar
# Eliminar la serie de datos inventada
tabla_resultados_locf <- tabla_resultados_locf[-1,]

for (i in 1:cantidad) {
  # 1. Crear y - Con datos faltantes (NA) al azar en la serie original con los NA sacados
  y_na <- set_na(serie.temporal, na=  c(serie.temporal[sample(gap1,1):sample(gap1,1)],serie.temporal[sample(gap2,1):sample(gap2,1)],serie.temporal[sample(gap3,1):sample(gap3,1)],serie.temporal[sample(gap4,1):sample(gap4,1)],serie.temporal[sample(gap5,1):sample(gap5,1)],serie.temporal[sample(gap6,1):sample(gap6,1)],serie.temporal[sample(gap7,1):sample(gap7,1)],serie.temporal[sample(gap8,1):sample(gap8,1)]))

  # 2. Calcular los valores a imputar con moving average exponential
  y_impute_locf <- na_locf(y_na,option = "locf", na_remaining = "rev")

  # 3. Calcular las metricas de desempe?o
  nas_locf <-sum(is.na(y_na))
  rmse_locf <- rmse(serie.temporal,y_impute_locf)
  cor_locf <- cor(serie.temporal,y_impute_locf)^2
  mase_locf <- mase(serie.temporal,y_impute_locf)
  smape_locf <- smape(serie.temporal,y_impute_locf)


  # 4. Almacenar en una tabla temporal las metricas
  temporal1 <- data.frame(nas_locf,rmse_locf, cor_locf,mase_locf,smape_locf)
  colnames(temporal1)=c("nas_locf","rmse_locf","cor_locf","mase_locf","smape_locf")

  # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
  tabla_resultados_locf <- bind_rows(tabla_resultados_locf, temporal1)
}

#crear tabla en archivo txt
write.table(tabla_resultados_locf,"tabla_resultados_locf2.txt")





#### PARA NA.INTERP con STL decomposition


tabla_resultados_stl <- data.frame(1,2,3,4,5) #Creo una tabla con una serie de datos
#Colocarle nombres a las columnas de la tabla
colnames(tabla_resultados_stl) <- c('nas_stl','rmse_stl', 'cor_stl','mase_stl','smape_stl')
serie.temporal <- pelada_imput_stl # CAMBIAR AQUI la serie con la que voy a trabajar
# Eliminar la serie de datos inventada
tabla_resultados_stl <- tabla_resultados_stl[-1,]

for (i in 1:cantidad) {
  # 1. Crear y - Con datos faltantes (NA) al azar en la serie original con los NA sacados
  y_na <- set_na(serie.temporal, na=  c(serie.temporal[sample(gap1,1):sample(gap1,1)],serie.temporal[sample(gap2,1):sample(gap2,1)],serie.temporal[sample(gap3,1):sample(gap3,1)],serie.temporal[sample(gap4,1):sample(gap4,1)],serie.temporal[sample(gap5,1):sample(gap5,1)],serie.temporal[sample(gap6,1):sample(gap6,1)],serie.temporal[sample(gap7,1):sample(gap7,1)],serie.temporal[sample(gap8,1):sample(gap8,1)]))

  # 2. Calcular los valores a imputar con moving average exponential
  y_impute_stl <- na.interp(y_na)

  # 3. Calcular las metricas de desempe?o
  nas_stl <-sum(is.na(y_na))
  rmse_stl <- rmse(serie.temporal,y_impute_stl)
  cor_stl <- cor(serie.temporal,y_impute_stl)^2
  mase_stl <- mase(serie.temporal,y_impute_stl)
  smape_stl <- smape(serie.temporal,y_impute_stl)


  # 4. Almacenar en una tabla temporal las metricas
  temporal1 <- data.frame(nas_stl,rmse_stl, cor_stl,mase_stl,smape_stl)
  colnames(temporal1)=c("nas_stl","rmse_stl","cor_stl","mase_stl","smape_stl")

  # 5. Agregar las metricas de la tabla temporal, a la tabla de resultados final
  tabla_resultados_stl <- bind_rows(tabla_resultados_stl, temporal1)
}

#crear tabla en archivo txt
write.table(tabla_resultados_stl,"tabla_resultados_stl2.txt")












#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#Load data and create a whole table for each species

val_st=read.table("tabla_resultados_st2.txt",h=T);val_st_o=order(val_st$nas_st);val_st=val_st[val_st_o,]
val_aa=read.table("tabla_resultados_aa2.txt",h=T);val_aa_o=order(val_aa$nas_aa);val_aa=val_aa[val_aa_o,]
val_itl=read.table("tabla_resultados_itl2.txt",h=T);val_itl_o=order(val_itl$nas_itl);val_itl=val_itl[val_itl_o,]
val_its=read.table("tabla_resultados_its2.txt",h=T);val_its_o=order(val_its$nas_its);val_its=val_its[val_its_o,]
val_itst=read.table("tabla_resultados_itst2.txt",h=T);val_itst_o=order(val_itst$nas_itst);val_itst=val_itst[val_itst_o,]
val_mas=read.table("tabla_resultados_mas2.txt",h=T);val_mas_o=order(val_mas$nas_mas);val_mas=val_mas[val_mas_o,]
val_mal=read.table("tabla_resultados_mal2.txt",h=T);val_mal_o=order(val_mal$nas_mal);val_mal=val_mal[val_mal_o,]
val_mae=read.table("tabla_resultados_mae2.txt",h=T);val_mae_o=order(val_mae$nas_mae);val_mae=val_mae[val_mae_o,]
val_seadec=read.table("tabla_resultados_seadec2.txt",h=T);val_seadec_o=order(val_seadec$nas_seadec);val_seadec=val_seadec[val_seadec_o,]
val_locf=read.table("tabla_resultados_locf2.txt",h=T);val_locf_o=order(val_locf$nas_locf);val_locf=val_locf[val_locf_o,]
val_stl=read.table("tabla_resultados_stl2.txt",h=T);val_stl_o=order(val_stl$nas_stl);val_stl=val_stl[val_stl_o,]

tabla=data.frame(val_aa$nas_aa,val_aa$rmse_aa,
                 val_aa$cor_aa,val_aa$mase_aa,val_aa$smape_aa,val_itl$nas_itl,
                 val_itl$rmse_itl,val_itl$cor_itl,val_itl$mase_itl,
                 val_itl$smape_itl,val_its$nas_its,val_its$rmse_its,val_its$cor_its,
                 val_its$mase_its,val_its$smape_its,val_itst$nas_itst,val_itst$rmse_itst,
                 val_itst$cor_itst,val_itst$mase_itst,val_itst$smape_itst,
                 val_mae$nas_mae,val_mae$rmse_mae,val_mae$cor_mae,val_mae$mase_mae,val_mae$smape_mae,
                 val_mal$nas_mal,val_mal$rmse_mal,val_mal$cor_mal,val_mal$mase_mal,
                 val_mal$smape_mal,val_mas$nas_mas,val_mas$rmse_mas,val_mas$cor_mas,
                 val_mas$mase_mas,val_mas$smape_mas,val_st$nas_st,val_st$rmse_st,
                 val_st$cor_st,val_st$mase_st,val_st$smape_st,val_seadec$nas_seadec,val_seadec$rmse_seadec,
                 val_seadec$cor_seadec,val_seadec$mase_seadec,val_seadec$smape_seadec,
                 val_locf$nas_locf,val_locf$rmse_locf,val_locf$cor_locf,val_locf$mase_locf,val_locf$smape_locf,
                 val_stl$nas_stl,val_stl$rmse_stl,val_stl$cor_stl,val_stl$mase_stl,val_stl$smape_stl)

colnames(tabla)=c("nas_aa","aa_rmse","aa_cor","aa_mase","aa_smape","nas_itl","itl_rmse",
                  "itl_cor","itl_mase","itl_smape","nas_its","its_rmse","its_cor","its_mase",
                  "its_smape","nas_itst","itst_rmse","itst_cor","itst_mase","itst_smape","nas_mae",
                  "mae_rmse","mae_cor","mae_mase","mae_smape","nas_mal","mal_rmse","mal_cor","mal_mase",
                  "mal_smape","nas_mas","mas_rmse","mas_cor","mas_mase","mas_smape","nas_st","st_rmse",
                  "st_cor","st_mase","st_smape","nas_seadec","seadec_rmse","seadec_cor","seadec_mase",
                  "seadec_smape","nas_locf","locf_rmse","locf_cor","locf_mase","locf_smape",
                  "nas_stl","stl_rmse","stl_cor","stl_mase","stl_smape")


write.table(tabla,"tabla_final2.txt")






