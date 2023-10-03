#' Linearly interpolate timeseries t(x) onto new timesetep ti
#'
#' Linear interpolation limited by time since previous or following reading
#' @param t time for x (time, POSIXct)
#' @param x any quantity, for example discharge (double)
#' @param ti time where t(x) will be interpolated to (time, POSIXct)
#' @param threshold maximum duration where interpolation is allowed (hours)
#' @returns a data frame with time (ti), x interpolated from t(x) onto ti, and logical (ibad) if interpolation exceeded threshold
#' @examples
#' InputData <- realTimeloads::ExampleData
#' ADCP <- InputData$ADCP
#' Height <- InputData$Height
#' # Interpolate river height to ADCP time
#' t <- realTimeloads::ExampleData$Height$time
#' x <- realTimeloads::ExampleData$Height$Height_m
#' ti <-realTimeloads::ExampleData$ADCP$time
#' threshold <- 1
#' Output<- linear_interpolation_with_time_limit(t,x,ti,threshold)
#' @author Daniel Livsey, September (2023), livsey.daniel@@gmail.com, ORCID: 0000-0002-2028-6128
#' @export
#'
linear_interpolation_with_time_limit <- function(t,x,ti,threshold) {
  # linearly interpolate timeseries x(t) onto ti
  # Assumes all time is in POSIXct
  # threshold is limit of interpolation to nearest read in minutes
  # Daniel Livsey, 20/07/2023

    # when duplicate "t" occurs code does not work as expected
    ind<-!duplicated(t)
    t<-t[ind]
    x<-x[ind]

    x_interpolated <- approx(t,x,ti)$y

    d1 <- data.table::data.table(t)
    d2 <- data.table::data.table(ti)
    d1$RollDate<-t
    d2$RollDate<-ti
    data.table::setkey(d1,"RollDate")
    data.table::setkey(d2,"RollDate")
    rb<-d1[d2, roll = -Inf]
    rf<-d1[d2, roll = Inf]
    dtb <- as.double(abs(difftime(rb$ti,rb$t,units="hours")))
    dtf <- as.double(abs(difftime(rf$ti,rf$t,units="hours")))
    # prevent interpolation when nearest finite value is not within X hours
    ibad <- dtb>threshold | dtf>threshold

    x_interpolated[ibad] <- NA

    # Output table
    time <- ti
    Output <- data.frame(time,x_interpolated,ibad)
    return(Output)
}
