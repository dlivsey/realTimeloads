#' Computes sediment load from optical and acoustic backscatter measurements
#'
#' Computes sediment load per guideline from optical and acoustic backscatter measurements combined to the "Sediment Composition Index" (SCI) see Livsey et al (2023)
#' @returns total load with uncertainty computed from estimates of concentration from SCI
#' @seealso
#' \code{\link{realTimeloads}} Package help file
#' @examples
#' \dontrun{
#' Output <- ExampleCodeSCI()
#' }
#' @references GUIDLINE
#' Livsey et al (2023): https://doi.org/10.1029/2022WR033982
#' @author Daniel Livsey, September (2023), livsey.daniel@@gmail.com, ORCID: 0000-0002-2028-6128
#' @export
#'
ExampleCodeSCI <- function() {
# Load data
Site <- realTimeloads::ExampleData$Site
ADCP <- realTimeloads::ExampleData$ADCP
Echo_Intensity_Beam_1 <- realTimeloads::ExampleData$Echo_Intensity
# # example code assumes backscatter is equal across beams
Echo_Intensity_Beam_2 <- realTimeloads::ExampleData$Echo_Intensity
Sonde <- realTimeloads::ExampleData$Sonde
Height <- realTimeloads::ExampleData$Height
Discharge <- realTimeloads::ExampleData$Discharge

### Process acoustic backscatter data  --------------------------------------
ADCPOutput <- acoustic_backscatter_processing(Site,ADCP,Height,Sonde,Echo_Intensity_Beam_1,Echo_Intensity_Beam_2)

# Compute SCI and Kt
# timeseries
to <- realTimeloads::ExampleData$Sonde$time
Co <- realTimeloads::ExampleData$Sediment_Samples$SSCpt_mg_per_liter
OBSo <- realTimeloads::ExampleData$Sonde$Turbidity_FNU
SNRo <- ADCPOutput$mean_sediment_corrected_backscatter_beam_1_and_beam_2_dB
Kto <- SNRo - 10*log10(Co)
SCIo <- SNRo - 10*log10(OBSo)

# calibration, randomly sample synthetic data to simulate sampling of hydrograph
n <- 12 # arbitrary number
ind <-sample(1:nrow(ADCP),n,replace=FALSE)
t<- realTimeloads::ExampleData$Sonde$time[ind]
OBS <- realTimeloads::ExampleData$Sonde$Turbidity_FNU[ind]
C <- realTimeloads::ExampleData$Sediment_Samples$SSCpt_mg_per_liter[ind]
Cxs <- realTimeloads::ExampleData$Sediment_Samples$SSCxs_mg_per_liter[ind]
SNR <- ADCPOutput$mean_sediment_corrected_backscatter_beam_1_and_beam_2_dB[ind]
Kt <- SNR - 10*log10(C)
SCI <- SNR - 10*log10(OBS)

# Fit C(dB)
bootC <- bootstrap_regression(data.frame(SNR,C),'log10(C)~SNR')
# Fit Kt(SCI)
bootKt <- bootstrap_regression(data.frame(SCI,Kt),'Kt~SCI')

# predict C(dB) and C(dB,SCI) w/ uncertainty
# C(dB)
iters <- nrow(bootC$regression_parameters)
regC <- bootC$regression_parameters
resC <- matrix(sample(bootC$fit$residuals,iters*length(to),replace=TRUE),nrow = iters,ncol=length(to)) # randomly sample residuals

# C(dB,SCI)
regKt <- bootKt$regression_parameters
resKt <- matrix(sample(bootC$fit$residuals,iters*length(to),replace=TRUE),nrow = iters,ncol=length(to)) # randomly sample residuals

# used for predictions in for loop below
fitC <- bootC$fit
dfC <- data.frame(SNRo)
names(dfC)<- 'SNR'
nfitC <-names(fitC$coefficient)
fitKt <- bootKt$fit
dfKt <- data.frame(SCIo)
names(dfKt)<- 'SCI'
nfitKt <-names(fitKt$coefficient)


Cpoi =  matrix(NA, nrow = iters, ncol = length(SCIo))
Ktpi =  matrix(NA, nrow = iters, ncol = length(SCIo))
Cpi =  matrix(NA, nrow = iters, ncol = length(SCIo))
for (i in 1:iters) {
  # C(dB)
  fitC$coefficient[nfitC] <- as.numeric(bootC$regression_parameters[i,]) # update coefficients for each i
  Cpoi[i,] <- 10^(predict(fitC,dfC)+resC[i,])*bootC$BCF[i]

  # C(dB,SCI)
  fitKt$coefficient[nfitKt] <- as.numeric(bootKt$regression_parameters[i,]) # update coefficients for each i
  Ktpi[i,] <- predict(fitKt,dfKt) + resKt[i,]
  Cpi[i,] <- 10^(0.1*SNRo-0.1*Ktpi[i,])
}
quants <- c(0.0527, 0.1587, 0.5, 0.8414, 0.9473) # +/- 1 and 2 sigma and median (i.e., reported) estimate

# quantiles on predicted C
Cpo <- t(apply(Cpoi , 2 , quantile , probs = quants , na.rm = TRUE ))
Cp <- t(apply(Cpi , 2 , quantile , probs = quants , na.rm = TRUE ))

# check outputs
#plot(Co,col = "black",ylim = c(100,1600))
#lines(Cpo[,1],col = "red")
#lines(Cpo[,5],col = "red")
#lines(Cp[,1],col = "blue")
#lines(Cp[,5],col = "blue")

# Note C(dB) missing peak in C while C(dB,SCI) does note
plot(to,Co,col='black',ylim = c(100,1600),xlab = "time (AEST)",ylab="Analyte concentration (mg/l)",main = "Estimated versus actual concentration")
lines(to,Cp[,3],col='red')
lines(to,Cpo[,3],col=c(4))
legend("topright",legend = c("Measured concentration (C)", "C predicted from dB and SCI", "C predicted from dB"), lty = c(0,1,1),col = c(1, 2,4),pch = c(1,-1,-1))

# compute load Qs (kt)
# compute dt (second) for load integration
dt =  c()
dt[2:length(to)] <- as.numeric(difftime(to[2:length(to)],to[1:length(to)-1],units = "secs"))
dt[1] = median(dt,na.rm=TRUE) # assume time step 1 using median dt

Qsoi =  matrix(NA, nrow = iters, ncol = length(to))
Qsi =  matrix(NA, nrow = iters, ncol = length(to))
Q <- Discharge$Discharge_m_cubed_per_s
for (i in 1:iters) {
  Qsi[i,] <- Cpi[i,]*Q*dt*1e-9 # assumes concentration is mg/l and discharge is cubic meters per sec, dt is in seconds
  Qsoi[i,] <- Cpoi[i,]*Q*dt*1e-9 # assumes concentration is mg/l and discharge is cubic meters per sec, dt is in seconds
  }

### Compute uncertainty on concentration and load timeseries ---------------
cQs = rowSums(Qsi,na.rm=TRUE) # total load for each iteration
Total_load_kt_from_SCI <- data.frame(t(quantile(cQs,quants,na.rm=TRUE)))
colnames(Total_load_kt_from_SCI) = c('minus_two_sigma_confidence','minus_one_sigma_confidence','median_confidence','plus_one_sigma_confidence','plus_two_sigma_confidence')

cQso = rowSums(Qsoi,na.rm=TRUE) # total load for each iteration
Total_load_kt_from_dB <- data.frame(t(quantile(cQso,quants,na.rm=TRUE)))
colnames(Total_load_kt_from_dB) = c('minus_two_sigma_confidence','minus_one_sigma_confidence','median_confidence','plus_one_sigma_confidence','plus_two_sigma_confidence')

# Percent difference in uncertainty quantiles relative to actual load
(Total_load_kt_from_SCI-Total_load_kt_from_dB)/sum(Co*Q*dt*1e-9,na.rm=TRUE)*100

return(Total_load_kt_from_SCI)
}
