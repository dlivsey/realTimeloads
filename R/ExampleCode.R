#' Computes sediment load per guideline from ExampleData
#'
#' Computes sediment load per guideline from ExampleData
#' @returns list with data frames of estimated concentration and flux along with data used in regression and surrogate timeseries
#' @seealso
#' \code{\link{realTimeloads}} Package help file
#' @examples
#' \dontrun{
#' Output <- ExampleCode()
#' }
#' @references GUIDLINE
#' @author Daniel Livsey, September (2023), livsey.daniel@@gmail.com, ORCID: 0000-0002-2028-6128
#' @export
#'
ExampleCode <- function() {
### Load data
InputData <- realTimeloads::ExampleData
Site <- InputData$Site
ADCP <- InputData$ADCP
Echo_Intensity_Beam_1 <- InputData$Echo_Intensity
Echo_Intensity_Beam_2 <- InputData$Echo_Intensity # example code assumes backscatter is equal across beams
Sonde <- InputData$Sonde
Height <- InputData$Height
Discharge <- InputData$Discharge
Sediment_Samples <- InputData$Sediment_Samples

### Process acoustic backscatter data  --------------------------------------
ADCPOutput <- acoustic_backscatter_processing(Site,ADCP,Height,Sonde,Echo_Intensity_Beam_1,Echo_Intensity_Beam_2)


### Compute estimate of analyte timeseries using https://doi.org/10.1029/2007WR006088 --------

# threshold (minutes) for interpolation of Surrogate to Analyte timeseries
threshold <- 30

# Surrogate "X" (e.g, acoustic backscatter and/or turbidity)
tx <- ADCPOutput$time
x <-ADCPOutput$mean_sediment_corrected_backscatter_beam_1_and_beam_2_dB

# Analyte "y" (e.g., TSS, NOx, etc)
ty<-Sediment_Samples$time
y<-Sediment_Samples$SSCpt_mg_per_liter

# interpolate surrogate onto analyte timeseries
# Randomly sample analyte timeseries n times to simulate pump samples
n <- 100
ind_calibration <-sample(1:length(ty),n,replace=FALSE)
calibration <- surrogate_to_analyte_interpolation(tx,x,ty[ind_calibration],y[ind_calibration],30)

# Specify calibration variable names for later storage in regression data
names(calibration)[grepl('surrogate',names(calibration))] <- 'SNR_dB'
names(calibration)[grepl('analyte',names(calibration))] <- 'SSCpt_mg_per_liter'

# compute bootstrap regression parameters per https://doi.org/10.1029/2007WR006088
boot <- bootstrap_regression(calibration,"log10(SSCpt_mg_per_liter) ~ SNR_dB+I(SNR_dB^2)")

res <- as.numeric(boot$fit$residuals)
iters <- nrow(boot$regression_parameters)
resC <- matrix(sample(res,iters*length(tx),replace=TRUE),nrow = iters,ncol=length(tx)) # randomly sample residuals
yts =  matrix(NA, nrow = iters, ncol = length(tx)) # each row provides an estimated timeseries of the analyte

fit <- boot$fit
df <- data.frame(x)
names(df) <- 'SNR_dB'
nfit <-names(fit$coefficient)

for (i in 1:iters) {
  #fit$coefficient[nfit] <- as.numeric(boot$regression_parameters[i,]) # update coefficients for each i
  yts[i,] <- 10^(predict(fit,df)+resC[i,])*boot$BCF[i]
}
rm(fit,nfit)

### calculate load (kt) ------------------------------------------------------------
# Interpolate discharge onto surrogate time series (e.g., Backscatter t.s)
output <- linear_interpolation_with_time_limit(Discharge$time,Discharge$Discharge_m_cubed_per_s,tx,threshold)
Q <- output$x_interpolated

# compute dt (second) for load integration
dt =  c()
dt[2:length(tx)] <- as.numeric(difftime(tx[2:length(tx)],tx[1:length(tx)-1],units = "secs"))
dt[1] = median(dt,na.rm=TRUE) # assume time step 1 using median dt

# compute load Qs (kt)
Qsi =  matrix(NA, nrow = iters, ncol = length(x))
for (i in 1:iters) {
  Qsi[i,] <- yts[i,]*Q*dt*1e-9 # assumes concentration is mg/l and discharge is cubic meters per sec, dt is in seconds
}

### Compute uncertainty on concentration and load timeseries ---------------
cQs = rowSums(Qsi,na.rm=TRUE) # total load for each iteration
quants <- c(0.0527, 0.1587, 0.5, 0.8414, 0.9473) # +/- 1 and 2 sigma and median (i.e., reported) estimate

Total_load_kt <- data.frame(t(quantile(cQs,quants,na.rm=TRUE)))
colnames(Total_load_kt) = c('minus_two_sigma_confidence','minus_one_sigma_confidence','median_confidence','plus_one_sigma_confidence','plus_two_sigma_confidence')

Reported_real_time_load_estimate <- Total_load_kt$median_confidence

Analyte_concentration_timeseries_mg_per_liter <- data.frame(t(apply(yts , 2 , quantile , probs = quants , na.rm = TRUE )))
colnames(Analyte_concentration_timeseries_mg_per_liter) = c('minus_two_sigma_confidence','minus_one_sigma_confidence','median_confidence','plus_one_sigma_confidence','plus_two_sigma_confidence')

Analyte_flux_timeseries_kt <- data.frame(t(apply(Qsi, 2 , quantile , probs = quants , na.rm = TRUE )))
colnames(Analyte_flux_timeseries_kt) = c('minus_two_sigma_confidence','minus_one_sigma_confidence','median_confidence','plus_one_sigma_confidence','plus_two_sigma_confidence')


### Plots and calculation checks ----------------------------------------------------

# actual loads (kt)
Cxs <- Sediment_Samples$SSCxs_mg_per_liter
Cpt <- Sediment_Samples$SSCpt_mg_per_liter
Qspt <- Cpt*Q*dt*1e-9
Qsxs <- Cxs*Q*dt*1e-9

indp <- seq(from = 1, to = nrow(ADCP), by = 1) # vector for plotting subset of timeseries

# concentration plot
plot(tx[indp],Analyte_concentration_timeseries_mg_per_liter$median_confidence[indp],col='red',type = "l",pch = c(21),xlab = "time (AEST)",ylab="Analyte concentration (mg/l)",main = "Estimated versus actual concentration")
#lines(tx[indp],Analyte_concentration_timeseries_mg_per_liter$minus_two_sigma_confidence[indp],col='green')
#lines(tx[indp],Analyte_concentration_timeseries_mg_per_liter$plus_two_sigma_confidence[indp],col='green')
points(ty[ind_calibration],y[ind_calibration],col = 'black',pch = c(21))
legend("topright",legend = c("Estimated concentration", "Actual concentration"), lty = c(1, 0),col = c(2, 1),pch = c(-1,21))

# Flux plot
plot(tx[indp],Analyte_flux_timeseries_kt$median_confidence[indp]/dt[indp]*1e3,col='red',type = "l",pch = c(21),xlab = "time (AEST)",ylab="Analyte flux (ton per second)",main = "Estimated versus actual flux")
#lines(tx[indp],Analyte_flux_timeseries_kt$minus_two_sigma_confidence[indp]/dt[indp]*1e3,col='green')
#lines(tx[indp],Analyte_flux_timeseries_kt$plus_two_sigma_confidence[indp]/dt[indp]*1e3,col='green')
points(tx[ind_calibration],Qspt[ind_calibration]/dt[ind_calibration]*1e3,col = 'black',pch = c(21))
legend("topright",legend = c("Estimated analyte flux", "Actual analyte flux"), lty = c(1, 0),col = c(2, 1),pch = c(-1,21))

# check estimated load / modeled load
#mean(cQs)/sum(Qspt,na.rm=TRUE) # relative to Cpt (biased) load
#mean(cQs)/sum(Qsxs,na.rm=TRUE) # relative to Cxs (actual) load

Output <- list("time" = tx,"surrogate_timeseries_used_for_prediction_of_analyte"= df,"regression_data" = calibration,"regression_parameters_estimated_from_bootstrap_resampling" = boot,"Analyte_concentration_timeseries_mg_per_liter"= Analyte_concentration_timeseries_mg_per_liter,"Dicharge" = Discharge,"Analyte_flux_timeseries_kt" =Analyte_flux_timeseries_kt)

return(Output)
}
