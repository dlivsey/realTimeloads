#' Compute timeseries with uncertainty from bootstrap regression
#'
#' Compute uncertainty on timeseries from bootstrap regression after Rustomji and Wilkinson (2008)
#' @param Surrogate data frame with time (PosixCt) and surrogate(s) (x,...)
#' @param Regression data frame from bootstrap_regression() that determines analyte(surrogate)
#' @returns list with inputs and uncertainty on timeseries estimated from Regression
#' @examples
#' \dontrun{
#' Turbidity_FNU <- realTimeloads::ExampleData$Sonde$Turbidity
#' TSS_mg_per_l <- realTimeloads::ExampleData$Sediment_Samples$SSCpt_mg_per_liter
#' Calibration <- data.frame(Turbidity_FNU,TSS_mg_per_l)
#' time <- realTimeloads::ExampleData$Sonde$time
#' Surrogate <- data.frame(time,Turbidity_FNU)
#' Regression = bootstrap_regression(Calibration,'TSS_mg_per_l~Turbidity_FNU')
#' Output <- estimate_timeseries(Surrogate,Regression)
#' }
#' @references
#' Rustomji, P., & Wilkinson, S. N. (2008). Applying bootstrap resampling to quantify uncertainty in fluvial suspended sediment loads estimated using rating curves. Water resources research, 44(9).https://doi.org/10.1029/2007WR006088
#'
#' Helsel, D.R., Hirsch, R.M., Ryberg, K.R., Archfield, S.A., and Gilroy, E.J., 2020, #' Statistical methods in water resources: U.S. Geological Survey Techniques and Methods, book 4, chap. A3, 458 p. https://doi.org/10.3133/tm4a3
#' @author Daniel Livsey (2023) ORCID: 0000-0002-2028-6128
#' @export
#'
estimate_timeseries <- function(Surrogate,Regression) {

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

  ### Compute uncertainty on timeseries ---------------
  quants <- c(0.0527, 0.1587, 0.5, 0.8414, 0.9473) # +/- 1 and 2 sigma and median (i.e., reported) estimate

  estimated_timeseries_quantiles <- data.frame(t(apply(yts , 2 , quantile , probs = quants , na.rm = TRUE )))
  colnames(estimated_timeseries_quantiles) = c('minus_two_sigma_confidence','minus_one_sigma_confidence','median_confidence','plus_one_sigma_confidence','plus_two_sigma_confidence')

  Output <- list("time" = tx,"estimated_timeseries_quantiles"=estimated_timeseries_quantiles,"estimated_timeseries" = yts,"surrogate"= Surrogate,"regression_data" = Regression)

  return(Output)
}
