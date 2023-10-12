#' Compute load with uncertainty on concentration estimates
#'
#' Compute load with uncertainty on concentration estimates from bootstrap regression after Rustomji and Wilkinson (2008)
#' @param Surrogate data frame with time (PosixCt) and surrogate(s) (x,...)
#' @param Discharge data frame with time (PosixCt) and discharge in cubic meters per second
#' @param Regression data frame from bootstrap_regression() that determines analyte(surrogate)
#' @param period two element vector time (PosixCt) indicating period over which load is computed
#' @returns list with data frames of estimated concentration and flux used to compute load (i.e., the sum of flux)
#' @section Note:
#' Surrogate and Discharge time series can be on different time steps
#'
#' If period is NULL, computes load over time in Surrogate
#' @section Warning:
#' Discharge should be in cubic meters per second
#'
#' Analyte concentration estimated from surrogate should be in milligrams per second
#' @examples
#' \dontrun{
#' Turbidity_FNU <- realTimeloads::ExampleData$Sonde$Turbidity
#' TSS_mg_per_l <- realTimeloads::ExampleData$Sediment_Samples$SSCpt_mg_per_liter
#' Discharge <- realTimeloads::ExampleData$Discharge
#' Calibration <- data.frame(Turbidity_FNU,TSS_mg_per_l)
#' time <- realTimeloads::ExampleData$Sonde$time
#' Surrogate <- data.frame(time,Turbidity_FNU)
#' Regression = bootstrap_regression(Calibration,'TSS_mg_per_l~Turbidity_FNU')
#' period <- c(as.POSIXct("2000-02-16 AEST"),as.POSIXct("2000-03-16 AEST"))
#' Output <- compute_load(Surrogate,Discharge,Regression,period)
#' }
#' @references
#' Rustomji, P., & Wilkinson, S. N. (2008). Applying bootstrap resampling to quantify uncertainty in fluvial suspended sediment loads estimated using rating curves. Water resources research, 44(9).https://doi.org/10.1029/2007WR006088
#'
#' Helsel, D.R., Hirsch, R.M., Ryberg, K.R., Archfield, S.A., and Gilroy, E.J., 2020, #' Statistical methods in water resources: U.S. Geological Survey Techniques and Methods, book 4, chap. A3, 458 p. https://doi.org/10.3133/tm4a3
#' @author Daniel Livsey (2023) ORCID: 0000-0002-2028-6128
#' @export
#'
compute_load <- function(Surrogate,Discharge,Regression,period = NULL) {

if (is.null(period)) {
period <- as.POSIXct(c())
period[1] <- min(Surrogate$time)
period[2] <- max(Surrogate$time)
}

fit <- Regression$fit
tx <- Surrogate[,grepl('time',names(Surrogate))]
n <- names(Surrogate)
df <- data.frame(Surrogate[,n[!grepl('time',names(Surrogate))]])
colnames(df)<- n[!grepl('time',n)]
nfit <-names(fit$coefficient)
iters <- nrow(Regression$regression_parameters)
resC <- matrix(sample(Regression$fit$residuals,iters*length(tx),replace=TRUE),nrow = iters,ncol=length(tx)) # randomly sample residuals
yts =  matrix(NA, nrow = iters, ncol = length(tx)) # each row provides an estimated timeseries of the analyte

# y(x,...)
if (sum(grepl("BCF",names(Regression)))==0) {
  for (i in 1:iters) {
    fit$coefficient[nfit] <- as.numeric(Regression$regression_parameters[i,]) # update coefficients for each i
    yts[i,] <- predict(fit,df)+resC[i,]
  }
}

#log(y)(x,...)
if (sum(grepl("BCF",names(Regression)))==1) {
for (i in 1:iters) {
  fit$coefficient[nfit] <- as.numeric(Regression$regression_parameters[i,]) # update coefficients for each i
  yts[i,] <- 10^(predict(fit,df)+resC[i,])*Regression$BCF[i]
}
}

# Interpolate discharge onto surrogate timeseries
threshold<- 60
output <- linear_interpolation_with_time_limit(Discharge$time,Discharge$Discharge_m_cubed_per_s,Surrogate$time,threshold)
Q <- output$x_interpolated

# compute dt (second) for load integration
dt =  c()
dt[2:length(tx)] <- as.numeric(difftime(tx[2:length(tx)],tx[1:length(tx)-1],units = "secs"))
dt[1] = median(dt,na.rm=TRUE) # assume time step 1 using median

# compute load (i.e, load = C*Q*dt) from flux Qs (flux = C*Q)
Qsi =  matrix(NA, nrow = iters, ncol = length(tx))
for (i in 1:iters) {
  Qsi[i,] <- yts[i,]*Q*dt*1e-9 # assumes concentration is mg/l and discharge is cubic meters per sec, dt is in seconds
}

### Compute uncertainty on concentration and load timeseries ---------------
cQs = rowSums(Qsi,na.rm=TRUE) # total load for each iteration
quants <- c(0.0527, 0.1587, 0.5, 0.8414, 0.9473) # +/- 1 and 2 sigma and median (i.e., reported) estimate

Analyte_concentration_timeseries_mg_per_liter <- data.frame(t(apply(yts , 2 , quantile , probs = quants , na.rm = TRUE )))
colnames(Analyte_concentration_timeseries_mg_per_liter) = c('minus_two_sigma_confidence','minus_one_sigma_confidence','median_confidence','plus_one_sigma_confidence','plus_two_sigma_confidence')

Analyte_load_timeseries_kt <- data.frame(t(apply(Qsi, 2 , quantile , probs = quants , na.rm = TRUE )))
colnames(Analyte_load_timeseries_kt) = c('minus_two_sigma_confidence','minus_one_sigma_confidence','median_confidence','plus_one_sigma_confidence','plus_two_sigma_confidence')

# Compute load over period indicated
ind <- Surrogate$time>=min(period) & Surrogate$time<=max(period)
cQso = rowSums(Qsi[,ind],na.rm=TRUE) # total load for each iteration
Total_load_kt <- data.frame(t(quantile(cQso,quants,na.rm=TRUE)))
colnames(Total_load_kt) = c('minus_two_sigma_confidence','minus_one_sigma_confidence','median_confidence','plus_one_sigma_confidence','plus_two_sigma_confidence')


Output <- list("time" = tx,"Analyte_concentration_timeseries_mg_per_liter"=Analyte_concentration_timeseries_mg_per_liter,"Analyte_load_timeseries_kt"=Analyte_load_timeseries_kt,"surrogate"= Surrogate,"regression_data" = Regression,"Discharge"=Discharge,'period_over_which_load_is_computed'=period)

return(Output)
}
