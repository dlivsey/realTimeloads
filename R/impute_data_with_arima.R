#' Returns x with gaps imputed using forecast::auto.arima with optional uncertainty estimation using Monte Carlo resampling
#'
#' Returns x with gaps imputed using forecast::auto.arima with optional uncertainty estimation using Monte Carlo resampling. Uncertainty on imputed data is estimated using using Monte Carlo (MC) resampling adapting methods of Rustomji and Wilkinson (2008). forecast::auto.arima searches for best ARIMA model and allows for step-changes in the ARIMA model of x(Xreg)
#'
#' @param time time for x (time, POSIXct)
#' @param x any quantity (double)
#' @param Xreg additional predictors for ARIMA
#' @param ti time vector for interpolation (time, POSIXct)
#' @param MC number of Monte Carlo simulations for uncertainty estimation
#' @param ptrain proportion of data used for training and testing model
#' @returns list with x imputed at time or ti, if given. Uncertainty estimated from Monte Carlo simulations
#' @note
#' For infilling missing suspended-sediment concentration data (e.g., TSS or SSC) I suggest using the Froude number and discharge for Xreg. For tidally affected sites, consider using tidally filtered discharge.
#' @examples
#' # Impute non-tidal data
#' time <- realTimeloads::ExampleData$Sediment_Samples$time
#' xo <- realTimeloads::ExampleData$Sediment_Samples$SSCxs_mg_per_liter
#' Q <- realTimeloads::ExampleData$Discharge$Discharge_m_cubed_per_s
#' idata <- sample(1:length(xo),round(length(xo)*0.5),replace=FALSE)
#' x <- rep(NA,length(xo))
#' x[idata] <- xo[idata] # simulated samples
#' Output <- impute_data_with_arima(time,x,Xreg = Q)
#'
#' @author Daniel Livsey (2023) ORCID: 0000-0002-2028-6128
#' @references
#' Rustomji, P., & Wilkinson, S. N. (2008). Applying bootstrap resampling to quantify uncertainty in fluvial suspended sediment loads estimated using rating curves. Water resources research, 44(9).
#'
#'Hyndman R, Athanasopoulos G, Bergmeir C, Caceres G, Chhay L, O'Hara-Wild M, Petropoulos F, Razbash S, Wang E, Yasmeen F (2023). forecast: Forecasting functions for time series and linear models. R package version 8.21.1, https://pkg.robjhyndman.com/forecast/.
#'
#'Hyndman RJ, Khandakar Y (2008). “Automatic time series forecasting: the forecast package for R.” Journal of Statistical Software, 26(3), 1–22. doi:10.18637/jss.v027.i03.
#'
#' @export
#'
impute_data_with_arima <- function(time,x,Xreg = NULL,ti = NULL,MC = 1,ptrain = 0.8) {

  if (is.null(ti)) ti <- time

  if (length(unique(diff(ti)))>1) {
    msg <-'time (or ti if given) must be regualarly spaced for ARIMA, no imputation undertaken'
    stop(msg)
  }

  # ensure inputs do not have missing values, duplicate times, or unsorted time
  igood <- !duplicated(time)
  time <- time[igood]
  x <- x[igood]
  df <- data.frame(time,x)
  df <- df[order(df$time),]
  time<- df$time
  x<- df$x
  rm(df)
  ti <- sort(ti)

  x<- as.vector(x)
  # Xreg must have column names, calling cbind() gives default colnames

  Xreg <-cbind(Xreg)

  # linearly interpolate vector if needed and permit 1 hr to nearest data point, ti is equally spaced time-vector
  if (length(ti)==length(time)) {
    if (as.double(sum(ti-time))==0) xi <- x # if ti=time, no need to interpolate
    if (as.double(sum(ti-time))!=0) {
      xi <-as.double(linear_interpolation_with_time_limit(time,x,ti,1)$x_interpolated)
    }
  }
  if (length(ti)!=length(time)) {
    xi <-as.double(linear_interpolation_with_time_limit(time,x,ti,1)$x_interpolated)
  }

  missing_data <- !is.finite(xi)

      ### ARIMA prediction with bootstrap uncertainty estimation  ----

      nit <- -1
      predictions <- matrix(nrow=length(xi),ncol = MC)
      validation_data <- list()
      # run MC to get validation residuals
      for (i in 1:MC) {
        #print(i)
        ival <- is.element(xi,sample(xi[is.finite(xi)],round(sum(is.finite(xi))*(1-ptrain))))
        xtrain <- xi
        xtrain[ival] <- NA
        # train validation model
        mod <- forecast::auto.arima(xtrain,xreg = Xreg,stepwise=FALSE,approximation=FALSE,num.cores=1,parallel=TRUE)$model
        # get validation data
        kal <- stats::KalmanSmooth(xtrain, mod, nit = -1)
        erg <- kal$smooth
        karima <- erg[,,drop = FALSE] %*% as.matrix(mod$Z)

        validation_data[[i]] <- data.frame(validation_predictions = karima[ival],validation_data = xi[ival])

        # get predictions for all data
        kal <- stats::KalmanSmooth(xi, mod, nit = -1)
        erg <- kal$smooth
        karima <- erg[,,drop = FALSE] %*% as.matrix(mod$Z)

        predictions[,i] <- karima

        # # rm() temporary variables
        rm(kal)
        rm(erg)
        rm(karima)
      }

      validation_data <- do.call(rbind,validation_data)
      #plot(validation_data$validation_predictions,validation_data$validation_data)

      res <- validation_data$validation_data-validation_data$validation_predictions
      iters <- 2000
      yp <- rowMeans(predictions)

      quants <- c(0.0527, 0.1587, 0.5, 0.8414, 0.9473) # +/- 1 and 2 sigma and median (i.e., reported) estimate
      tse <- data.frame(matrix(nrow = length(xi),ncol=length(quants)))
      colnames(tse) = c('x_at_minus_two_sigma_confidence','x_at_minus_one_sigma_confidence','x_at_median_confidence','x_at_plus_one_sigma_confidence','x_at_plus_two_sigma_confidence')

      for (i in 1:length(yp)) {
        tse[i,] <- quantile(yp[i]+sample(res,iters,replace = TRUE),probs = quants , na.rm = TRUE)
      }


      ### format outputs to match realTimeloads::impute_data ----

      # find data gaps exceeding one day
      readings_per_hour <- round(60/as.double(difftime(ti[2],ti[1],units = 'mins')))
      zy <- imputeTS::na_interpolation(xi,maxgap = readings_per_hour*24)
      tse[is.na(zy),] <- NA

      x = tse$x_at_median_confidence
      x[is.finite(xi)] <- xi[is.finite(xi)]

      Imputed_data = data.frame(time = ti,x,imputed = !is.finite(xi),data_gap_exceeds_one_day = is.na(zy),tse)

      colnames(validation_data) <- c("validation_data","validation_data_predictions")
      fit_summary <- summary(lm('validation_data~validation_data_predictions',validation_data))
      validation_r_squared <- fit_summary$adj.r.squared
      # root- squared error (units of validation_data)
      validation_rmse <- sqrt(mean(fit_summary$residuals^2))
      # Model standard percentage error (MSPE) (Rasmussen et al., 2009) in units of validation_data
      validation_mspe <- validation_rmse/mean(validation_data$validation_data)*100

      Validation_data_and_statistics <- list(validation_data=validation_data,validation_r_squared=validation_r_squared,validation_rmse=validation_rmse,validation_mspe=validation_mspe,proportion_of_data_held_out_for_validation = 1-ptrain,number_of_Monte_Carlo_simulations_for_uncertainty_estimation = MC)

      predictors <- colnames(Xreg)

      imputation_code <- "realTimeloads::impute_data_with_arima, Livsey (2023)"


      Output <- list(Imputed_data=Imputed_data,Validation_data_and_statistics=Validation_data_and_statistics,predictors=predictors,imputation_code=imputation_code)


      return(Output)

    }

