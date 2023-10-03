#' Interpolate timeseries x(t) onto y(t)
#'
#' Interpolate timeseries x(t) onto y(t) with temporal threshold on interpolation
#' @param tx time for x "surrogate" (time, POSIXct)
#' @param x quantity used to estimate y, for example, accoustic backscatter
#' @param ty time for y "analyte" (time, POSIXct)
#' @param y measured quantity, for example, an analyte such as suspended-sediment concentration
#' @param threshold maximum duration where interpolation is allowed (minutes)
#' @returns a data frame with surrogate (x) interpolated onto timestep of analyte (y), interpolated values exceeding threshold are excluded from the output
#' @examples
#' tx <- as.POSIXct(seq(0,24*60^2,60*1), origin = "2000-01-01",tz = "Australia/Brisbane")
#' x <- sin(1:length(t))
#' ty <- as.POSIXct(seq(0,24*60^2,60*15), origin = "2000-01-01",tz = "Australia/Brisbane")
#' y <- seq(0,24*60^2,60*15)
#' threshold <- 10
#' calibration <- surrogate_to_analyte_interpolation(tx,x,ty,y,threshold)
#' @author Daniel Livsey, September (2023), livsey.daniel@@gmail.com, ORCID: 0000-0002-2028-6128
#' @export
#'
surrogate_to_analyte_interpolation <- function(tx,x,ty,y,threshold) {
  # Interpolate timeseries x(tx) onto y(ty)
  # Assumes all time is in POSIXct
  # threshold is limit of interpolation in minutes
  # Daniel Livsey, 20/07/2023

  # Example:
  # surrogate (e.g., turbidity, acoustic backsactter)
  #tx <- Processed_Backscatter$time
  #x <- Processed_Backscatter$Mean_sediment_corrected_backscatter_dB
  # analyte (e.g., concentration of suspended-sediment (SSC, TSS))
  #ty<-Sediment_Samples$time
  #y<-Sediment_Samples$SSCpt_mg_per_liter

  surrogate <- data.frame(tx,x)
  analyte <- data.frame(ty,y)

  # Interpolate data
  # For the code example, tx and ty are on the same time step. For some samples tx and ty may not be aligned, when tx and ty are not aligned one should ensure interpolation is not made beyond some time duration
  x_interpolated <- data.frame(approx(surrogate$tx,surrogate$x, xout = analyte$ty,
                                      rule = 1, method = "linear", ties = mean))

  # QAQC compute time to antecedent, nearest, and subsequent interpolated value
  dt =  matrix(NA, nrow = nrow(analyte), ncol = 3)
  for (i in 1:length(ty)) {
    st <- ty[i] # analyte time
    ind <- which(abs(tx - st) == min(abs(tx - st))) # index of nearest ADCP reading

    if (ind==1) {
      # time to nearest time stamp
      dt[i,2] <-as.numeric(difftime(ty[i], tx[ind], units = "mins"))
      # time to subsequent time stamp
      dt[i,3] <-as.numeric(difftime(tx[ind+1],ty[i], units = "mins"))
    }

    if (ind>1) {
      # time to antecedent time stamp
      dt[i,1] <-as.numeric(difftime(tx[ind-1],ty[i], units = "mins"))
      # time to nearest time stamp
      dt[i,2] <-as.numeric(difftime(ty[i], tx[ind], units = "mins"))
      # time to subsequent time stamp
      dt[i,3] <-as.numeric(difftime(tx[ind+1],ty[i], units = "mins"))
    }

    if (ind==length(tx)) {
      # time to antecedent time stamp
      dt[i,1] <-as.numeric(difftime(tx[ind-1],ty[i], units = "mins"))
      # time to nearest time stamp
      dt[i,2] <-as.numeric(difftime(ty[i], tx[ind], units = "mins"))
    }

  }

  # makes some rule for interpolation
  dt <- abs(dt)
  # reject pair if closest matching read is beyond X mins away and if antecedent or subsequent surrogate read is beyond X mins away
  ibad <- dt[,2]>threshold & pmax(dt[,1],dt[,3])>threshold
  # Number of rejected pairs
  Number_of_rejected_pairs <- sum(as.numeric(ibad))

  # Store table w/ acceptable pairs
  time <- ty[!ibad]
  surrogate <- x_interpolated$y[!ibad]
  analyte <- y[!ibad]

  # Table for regression
  calibration <- data.frame(time,surrogate,analyte)

  # Select only finite pairs
  calibration<-calibration[complete.cases(calibration),]

  return(calibration)
}
