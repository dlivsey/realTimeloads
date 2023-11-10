#' Returns x with gaps imputed using ARIMA and Decision Trees, optional uncertainty estimation using Monte Carlo resampling
#'
#' Returns x with gaps imputed using ARIMA and Decision Trees with option to use  harmonic model as predictors for x in decision tree algorithm. Uncertainty on imputed data is estimated using using Monte Carlo (MC) resampling adapting methods of Rustomji and Wilkinson (2008)
#'
#' @param time time for x (time, POSIXct)
#' @param x any quantity (double)
#' @param Xreg additional predictors for decision tree, required if harmonic is FALSE (rows = time, or if given, ti)
#' @param ti time vector for interpolation (time, POSIXct)
#' @param hfit model object from TideHarmonics::ftide
#' @param harmonic logical if x exhibits tidal or diurnal variability
#' @param only_use_Xreg logical for using Xreg only in decision tree
#' @param MC number of Monte Carlo simulations for uncertainty estimation
#' @param ptrain proportion of data used for training and testing model
#' @note If MC == 1, uncertainty is not evaluated. If ptrain == 1, uncertainty and validation accuracy are not computed
#' @returns list with x imputed at time or ti, if given. Uncertainty estimated from Monte Carlo simulations
#' @examples
#' # Impute non-tidal data
#' time <- realTimeloads::ExampleData$Sediment_Samples$time
#' xo <- realTimeloads::ExampleData$Sediment_Samples$SSCxs_mg_per_liter
#' Q <- realTimeloads::ExampleData$Discharge$Discharge_m_cubed_per_s
#' idata <- sample(1:length(xo),round(length(xo)*0.5),replace=FALSE)
#' x <- rep(NA,length(xo))
#' x[idata] <- xo[idata] # simulated samples
#' flow_concentrtion_ratio <- imputeTS::na_interpolation(Q/x)
#' Xreg <- cbind(Q,flow_concentrtion_ratio)
#' Output <- impute_data(time,x,Xreg,MC = 10,ptrain = 0.8)
#'
#' # Impute tidal data
#' time <-TideHarmonics::Portland$DateTime[1:(24*90)]
#' xo <-TideHarmonics::Portland$SeaLevel[1:(24*90)]
#' idata <- sample(1:length(xo),round(length(xo)*0.5),replace=FALSE)
#' x <- rep(NA,length(xo))
#' x[idata] <- xo[idata] # simulated samples
#' Output <- impute_data(time,x,harmonic = TRUE,MC = 10,ptrain = 0.8)
#' @author Daniel Livsey (2023) ORCID: 0000-0002-2028-6128
#' @references
#' Rustomji, P., & Wilkinson, S. N. (2008). Applying bootstrap resampling to quantify uncertainty in fluvial suspended sediment loads estimated using rating curves. Water resources research, 44(9).
#'
#' van Buuren S, Groothuis-Oudshoorn K (2011). “mice: Multivariate Imputation by Chained Equations in R.” Journal of Statistical Software, 45(3), 1-67. doi:10.18637/jss.v045.i03.
#'
#'Stephenson AG (2016). Harmonic Analysis of Tides Using TideHarmonics. https://CRAN.R-project.org/package=TideHarmonics.
#'
#'Moritz S, Bartz-Beielstein T (2017). “imputeTS: Time Series Missing Value Imputation in R.” The R Journal, 9(1), 207–218. doi:10.32614/RJ-2017-009.
#'
#' @export
#'
impute_data <- function(time,x,Xreg = NULL,ti = NULL,hfit = NULL,harmonic=FALSE,only_use_Xreg=FALSE,MC = 1,ptrain = 1) {

  #st <- Sys.time() # used in testing for inspecting run time of code blocks

  if (!harmonic & is.null(Xreg)) {
    msg <-'Xreg must be provided if harmonic is FALSE, no imputation undertaken'
    stop(msg)
  }

  if (is.null(ti)) ti <- time

  if (length(unique(diff(ti)))>1) {
    msg <-'time (or ti if given) must be regualarly spaced for ARIMA, no imputation undertaken'
    stop(msg)
  }

  # if user desires to use all data in training model no uncertainty or validation accuracy estimated
  if (ptrain==1 & MC !=1) {
    message('All data used in training model (ptrain=1), setting MC = 1, no uncertainty or validation accuracy estimated')
    MC = 1
  }

  Output <- NULL
  if (length(unique(diff(ti)))==1) {
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

    # x must be a vector for mice::mice.impute.cart
    x<- as.vector(x)
    # Xreg must have column names, calling cbind() gives default colnames
    Xreg <-cbind(Xreg)

    readings_per_hour <- round(60/as.double(difftime(ti[2],ti[1],units = 'mins')))
    readings_per_day <- readings_per_hour*24
    dt <- as.double(difftime(ti[2],ti[1],units = 'hours'))

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

    ##message(paste('number of finite data after linear interpolation: '),sum(is.finite(xi)))
    ##message(paste('number of missing data after linear interpolation: '),sum(!is.finite(xi)))

    ### Fit harmonic model for harmonic data ----
    if (harmonic & !only_use_Xreg) {

      igood <- is.finite(x) # ftide does not accept NA
      if (is.null(hfit)) hfit <- TideHarmonics::ftide(x[igood],time[igood],TideHarmonics::hc7)
      # sum of harmonic amplitudes and msl
      xp <- predict(hfit,from = min(ti),to = max(ti),by = dt)

      # individual harmonic amplitudes can be used if desired
      xph <- data.frame(t(predict(hfit,from = min(ti),to = max(ti),by = dt,split=TRUE)))
      # range for each harmonic
      xpr<-data.frame(sapply(xph, range)[2,]-sapply(xph, range)[1,])
      colnames(xpr) <- "Harmonic_Amplitude"
      iorder <-order(xpr$Harmonic_Amplitude, decreasing = TRUE)
      rn <- rownames(xpr)[iorder]
      xpr<-data.frame(xpr[iorder,]) # range for each harmonic sorted in descending amplitude
      colnames(xpr) <- "Harmonic_Amplitude_m"
      rownames(xpr) <- rn
      xph <- xph[,iorder] # order harmonic table timeseries by amplitude as well
      if (ncol(xph)>7) {
        ih <- xpr$Harmonic_Amplitude/max(xpr$Harmonic_Amplitude)>0.2
      }
      if (ncol(xph)<=7) {
        ih <- !vector("logical",ncol(xph))
      }
      # get estimate of non-tidal signal
      # use harmonic model predictions to infill missing data and reduce data lost to filter ringing in xnt
      #Xarima <- cbind(xp)
      #xii <- imputeTS::na_kalman(xi,model ="auto.arima",xreg = as.matrix(Xarima),maxgap = readings_per_hour*24)
      # call to: imputeTS::na_kalman(xi,model ="auto.arima",xreg = as.matrix(Xarima),maxgap = readings_per_hour*24), adds longer run-time, using linear interpolation to speed up code

      xii <- imputeTS::na_interpolation(xi,maxgap = readings_per_hour*24) # infill data to prevent filter ringing

      xnt <- butterworth_tidal_filter(ti,xii)  # Non-tidal stage using USGS butterworth
      xnt <- imputeTS::na_interpolation(xnt,maxgap = readings_per_hour*24*7)
      xnt[!is.finite(xnt)] <- median(xnt,na.rm=TRUE)

    }

    ### Impute up to 3 hr gaps using ARIMA ----
    # Decided to skip y(ARIMA(X)) and go to y(ARIMA(y))
    # if (is.null(X)) Xarima <- data.frame(xp) # one may use more than one predictor, e.g., an upstream gauging
    # if (!is.null(X)) Xarima <- data.frame(xp,X) # one may use more than one predictor, e.g., an upstream gauging
    #
    # # automatically searches for best arima using auto.arima in forecast
    # xf <- na_kalman(xi,model ="auto.arima",xreg = as.matrix(Xarima),maxgap = readings_per_hour*3)
    #
    # # found that na_kalman can predict spurious values at start and of record
    # dXdt <- abs(c(median(diff(xf)),diff(xf)))
    # anomalies <-is.element(dXdt,boxplot(dXdt,plot=FALSE)[["out"]])
    # # print(paste('Number of suspect data points =',sum(anomalies)))
    # # print(paste('Number of data points imputed =',sum(is.finite(xf))-sum(is.finite(xi))))
    # #remove suspect anomalies
    # xf[anomalies] <- NA

    ### Impute up to 3 hrs gaps using y(ARIMA(y)) ---
    xf <- imputeTS::na_kalman(xi,maxgap = readings_per_hour*3)


    ### Impute largest gaps using regression trees ----
    if (only_use_Xreg & !is.null(Xreg)) Xtree <- scale(Xreg)
    if (!only_use_Xreg & !is.null(Xreg) & !harmonic) Xtree <- scale(Xreg)
    if (!only_use_Xreg & is.null(Xreg) & harmonic) Xtree <- scale(cbind(xp,xnt,as.matrix(xph[,ih])))
    if (!only_use_Xreg & !is.null(Xreg) & harmonic) Xtree <- scale(cbind(Xreg,xp,xnt))

    # use "maxgap" argument to find gaps that exceed 1 day duration
    data_gap_exceeds_one_day <- !is.finite(imputeTS::na_kalman(xf,maxgap = readings_per_day*1))

    # remaining gaps
    missing_data_2 <- !is.finite(xf)

    # defunct code block
    # fit cart once and smooth imputed data if desired
    #xf[missing_data_2] <- mice::mice.impute.cart(y=xf,ry = !missing_data_2, x=Xtree, wy = missing_data_2, minbucket = 5)
    ## smooth imputed data from trees using moving average
    ## window is centered w/ k steps before and after reading
    #smoothing_window <- readings_per_hour*3
    #w <- rep(1/smoothing_window,smoothing_window) # simple moving average
    # w <- 1/(2^smoothing_window) # exponential moving average
    #xs <- as.matrix(stats::filter(xf,w, sides=2))
    ## do not smooth-over observed data
    #xs[!missing_data_2] <- xf[!missing_data_2]

    # testing run time of blocks
    ##message('time_to_tree:')
    ##print(Sys.time() - st)
    ##st <- Sys.time() # reset time

    # uncertainty is not estimated if MC or ptrain = 1
    uncertainty <-imputation_uncertainty(x=xf,Xtree=Xtree,MC=MC,ptrain)

    # testing run time of blocks
    ##message('tree_run_time_with_uncertainty:')
    ##print(Sys.time() - st)

    if (ptrain==1) {
      x_imputed = uncertainty$prediction$x_imputed
      x_imputed[!missing_data] <- xf[!missing_data] # ensure observed data is not overwritten by imputed estimates
      df <- data.frame(time = ti,x=x_imputed,imputed=missing_data,data_gap_exceeds_one_day)
    }

    if (ptrain!=1) {
      x_imputed = uncertainty$prediction$x_imputed
      x_imputed[!missing_data] <- xf[!missing_data] # ensure observed data is not overwritten by imputed estimates
      df <- cbind(data.frame(time = ti,x = x_imputed,imputed = missing_data,data_gap_exceeds_one_day),uncertainty$prediction[,names(uncertainty$prediction)!='x_imputed'])
    }

    Output <- list("Imputed_data"=df,"Validation_data_and_statistics"= uncertainty[names(uncertainty)!='prediction'],'predictors'= colnames(Xtree),'imputation_code'= 'realTimeLoads::impute_data, Livsey (2023)')


    return(Output)
  }

  return(Output)
}

### functions for uncertainty estimation ----
imputation_uncertainty <- function(x,Xtree,MC,ptrain) {

  # ptrain 0.7 gives: 70 percent train/test, 30 percent hold-out validation for each simulation
  n_train_test <- round(sum(is.finite(x))*ptrain)

  # limit max_attempts for large MC to ensure code does not run too long
  if (MC<=100) {
    max_attempts <- 100
  }
  if (MC>100) {
    max_attempts <- 10
  }
  if (ptrain==1) {
    max_attempts <- 1 # if using all data set max_attempts to 1
  }

  bootstrap_data <- list()
  estimates <- list()

  for (i in 1:MC) {
    #print(paste('MC simulation =',i))
    out <- NULL
    attempt <- 1
    # for each i resample data if mice::mice.impute.cart() fails
    while(is.null(out) && attempt <= max_attempts) {
      if (length(x)>10000 & MC >100) {
        message(paste(paste('Monte Carlo Simulation =',i),paste('mice::mice.impute.cart() attempt =',attempt)))
      }
      attempt <- attempt + 1
      try(out <- mc_model(x,n_train_test,Xtree,MC,ptrain))
    }

    # store timeseries
    estimates[i] <- list(out$xi)

    # store validation data from random sample
    df <- data.frame(validation_data=out$validation_data,validation_data_predictions=out$xi[out$ival])
    bootstrap_data[i] <- list(df)
  }

  if (ptrain==1) {
    # Validation R-squared, RMSE, and Model standard percentage error (MSPE) not applicable when ptrain == 1, used all data to train model
    validation_r_squared <- NULL
    # root- squared error (units of validation_data)
    validation_rmse <- NULL
    # Model standard percentage error (MSPE) (Rasmussen et al., 2009) in units of validation_data
    validation_mspe <- NULL
    # Time series with no uncertainty, length(estimates)==1
    tse <- data.frame(estimates)
    colnames(tse) <- 'x_imputed'
  }

  if (ptrain<1) {
    # Validation R-squared, RMSE, and Model standard percentage error (MSPE)
    df <- do.call(rbind,bootstrap_data)
    fit_summary <- summary(lm('validation_data~validation_data_predictions',df))
    validation_r_squared <- fit_summary$adj.r.squared
    # root- squared error (units of validation_data)
    validation_rmse <- sqrt(mean(fit_summary$residuals^2))
    # Model standard percentage error (MSPE) (Rasmussen et al., 2009) in units of validation_data
    validation_mspe <- validation_rmse/mean(df$validation_data)*100

    if (MC==1) {
      # Time series with no uncertainty, length(estimates)==1
      tse <- data.frame(estimates)
      colnames(tse) <- 'x_imputed'
    }

    # Time series with uncertainty
    if (MC>1) {
      quants <- c(0.0527, 0.1587, 0.5, 0.8414, 0.9473) # +/- 1 and 2 sigma and median (i.e., reported) estimate
      tse<-data.frame(t(apply(as.matrix(do.call(rbind,estimates)), 2 , quantile , probs = quants , na.rm = TRUE )))
      colnames(tse) = c('x_at_minus_two_sigma_confidence','x_at_minus_one_sigma_confidence','x_at_median_confidence','x_at_plus_one_sigma_confidence','x_at_plus_two_sigma_confidence')

      tse$x_imputed <- tse$x_at_median_confidence

    }
  }

  # Store outputs
  Output <- list('prediction' = tse,'validation_data'=df,'validation_r_squared'=validation_r_squared,'validation_rmse'= validation_rmse,'validation_mspe'=validation_mspe,'proportion_of_data_held_out_for_validation' = 1-ptrain,'number_of_Monte_Carlo_simulations_for_uncertainty_estimation' = MC)

  return(Output)
}
# sometimes optimization routines in mice::mice.impute.cart() can fail depending on data set, using try() above to resample if mice::mice.impute.cart() fails
mc_model <- function(x,n_train_test,Xtree,MC,ptrain) {
  # randomly sample
  train_test_data <- sample(x[is.finite(x)],n_train_test)
  itrain <- is.element(x,train_test_data)
  ival <- !itrain & is.finite(x)
  validation_data <- x[ival]
  # train and test model on itrain and predict at ival

  if (ptrain<1) {
  ##st <- Sys.time() # reset time
  xi <- mice::mice.impute.cart(y=x,ry = itrain, x=Xtree, wy = !logical(length(x)), minbucket = 5)
  ##message('time_to_tree_predicting_all:')
  ##print(Sys.time() - st)
  }

  if (ptrain==1) {
    ##st <- Sys.time() # reset time
    validation_data <- NULL
    xii <- mice::mice.impute.cart(y=x,ry = itrain, x=Xtree, wy = !is.finite(x), minbucket = 5)
    xi <- x
    xi[!is.finite(x)] <- xii
    ##message('time_to_tree_predicting_subset:')
    ##print(Sys.time() - st)
  }

  out <- list("xi"=xi,"ival"=ival,"validation_data"=validation_data)
  return(out)
}

