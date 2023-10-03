#' Load data from csv files and prepare ExampleData.rda
#'
#' Code used to import csv files to R, can be used as template to read-in csv files with same format.
#' @returns list for inputs needed in ExampleCode.R
#' @section Warning:
#' Package arguments require variable names and units to match the names and variable units provided in csv files
#'
#' Save csv files in .txt format to ensure time format is not changed when editing/saving csv in Excel
#' @examples
#' ExampleData <- ExampleCode()
#' @references GUIDLINE
#' @author Daniel Livsey, September (2023), livsey.daniel@@gmail.com, ORCID: 0000-0002-2028-6128
#' @export
#'
importCSVfiles <- function() {

# defunct argument used when txt files are at 6 min timestep in txt files used in guideline. to reduce package size, 6 min data was down-sampled to hourly data
#if  (missing(step)) {
#  step <- 9 # default to hourly data for example
#}


## use hourly data for example code in package, 6 minute example data used in GUIDLINE available on Github

# LOAD DATA -----------------------------------------------------
# Site
fpath <- system.file("extdata", "ADCP_elevations_and_site_datums.txt", package="realTimeloads")
Site <- read.csv(fpath)
Site$Start_date_and_time <- as.POSIXct(Site$Start_date_and_time,format = "%d-%b-%Y %H:%M:%S",tz = "Australia/Brisbane")
Site$End_date_and_time <- as.POSIXct(Site$End_date_and_time,format = "%d-%b-%Y %H:%M:%S",tz = "Australia/Brisbane")

# Stage referenced to local datum in Site data frame
fpath <- system.file("extdata", "Height.txt", package="realTimeloads")
Height <- read.csv(fpath)
Height$time <- as.POSIXct(Height$time,format = "%d-%b-%Y %H:%M:%S",tz = "Australia/Brisbane")

# Sediment concentration
fpath <- system.file("extdata", "Sediment_Samples.txt", package="realTimeloads")
Sediment_Samples <- read.csv(fpath)
Sediment_Samples$time <- as.POSIXct(Sediment_Samples$time,format = "%d-%b-%Y %H:%M:%S",tz = "Australia/Brisbane")

# Discharge data
fpath <- system.file("extdata", "Discharge.txt", package="realTimeloads")
Discharge <- read.csv(fpath)
Discharge$time <- as.POSIXct(Discharge$time,format = "%d-%b-%Y %H:%M:%S",tz = "Australia/Brisbane")

# Sonde data
fpath <- system.file("extdata", "Sonde.txt", package="realTimeloads")
Sonde <- read.csv(fpath)
Sonde$time <- as.POSIXct(Sonde$time,format = "%d-%b-%Y %H:%M:%S",tz = "Australia/Brisbane")

# ADCP data except backscatter
fpath <- system.file("extdata", "ADCP_Data.txt", package="realTimeloads")
ADCP <- read.csv(fpath)
ADCP$time <- as.POSIXct(ADCP$time,format = "%d-%b-%Y %H:%M:%S",tz = "Australia/Brisbane")

# Acoustic backscatter
fpath <- system.file("extdata", "ADCP_Echo_Intensity.txt", package="realTimeloads")
Echo_Intensity_Beam_1 <- read.csv(fpath)
Echo_Intensity_Beam_1$time <- as.POSIXct(Echo_Intensity_Beam_1$time,format = "%d-%b-%Y %H:%M:%S",tz = "Australia/Brisbane")

# For example, beam 1 and beam 2 have equal backscatter, however for actual deployments counts in beam 1 and 2 will differ
Echo_Intensity_Beam_2 <- Echo_Intensity_Beam_1

# Export data -------------------------------------------------------------

# modify if one desires data lower frequency
step <- 1
indices_for_Example_data <- seq(1,nrow(ADCP),step)
ind <- indices_for_Example_data

ExampleData <- list("Site"=Site,"ADCP"=ADCP[ind,],"Echo_Intensity"=Echo_Intensity_Beam_1[ind,],"Sonde"=Sonde[ind,],"Height" = Height[ind,],"Discharge"=Discharge[ind,],"Sediment_Samples"=Sediment_Samples[ind,])

# Update ExampleData.rda data !!! for package development only !!!
#usethis::use_data(ExampleData, overwrite = TRUE) # !!! for package development only !!!

return(ExampleData)
}
