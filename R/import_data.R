#' Load data from comma-delimited .txt files to list to be used in function hADCPLoads()
#'
#' Imports csv files to R, file names, variable names (and units) in csv text files must match variable names used in ExampleData.rda
#' @param data_folder file path to folder containing .txt csv files with format that matches files in extdata package folder
#' @returns list with data frames used in package code, see ?ExampleData for list format
#' @section Warning:
#' Synthetic data used in ExampleData only has backscatter for one beam ("ADCP_Echo_Intensity.txt"), for user data, one should have backscatter for two beams with following names: "ADCP_Echo_Intensity_Beam_1.txt" and "ADCP_Echo_Intensity_Beam_2.txt"
#'
#' Package arguments require variable names and units to match the names and variable units provided (see ?ExampleData, or .txt files in extdata folder)
#'
#' Suggest saving all csv files in .txt format to ensure time format is not changed when editing/saving csv in Excel
#' @examples
#' \dontrun{
#' InputData <- import_data() # loads text files provided in package folder "extdata"
#' }
#' @seealso
#' \code{\link{hADCPLoads}} Process acoustic backscatter from hADCP and compute load using InputData from import_Data()
#' @references GUIDLINE
#' @author Daniel Livsey, September (2023), livsey.daniel@@gmail.com, ORCID: 0000-0002-2028-6128
#' @export
#'
import_data <- function(data_folder) {
  # LOAD DATA -----------------------------------------------------
  # file paths
  # Load text files provided in package if no data_folder provided
  if (missing(data_folder)) {
    print("LOADING PACKAGE EXAMPLE DATA")

    fpathSite <- system.file("extdata", "ADCP_elevations_and_site_datums.txt", package="realTimeloads")

    fpathHeight <- system.file("extdata", "Height.txt", package="realTimeloads")

    fpathSed <- system.file("extdata", "Sediment_Samples.txt", package="realTimeloads")

    fpathDischarge <- system.file("extdata", "Discharge.txt", package="realTimeloads")

    fpathSonde <- system.file("extdata", "Sonde.txt", package="realTimeloads")

    fpathADCP <- system.file("extdata", "ADCP_Data.txt", package="realTimeloads")

    fpathADCPBackscatterBeam1 <- system.file("extdata", "ADCP_Echo_Intensity.txt", package="realTimeloads")

    fpathADCPBackscatterBeam2 <- system.file("extdata", "ADCP_Echo_Intensity.txt", package="realTimeloads")
  }
  # Load text files provided in package if data_folder provided
  if (!missing(data_folder)) {
    fpathSite <- paste(data_folder,"/ADCP_elevations_and_site_datums.txt",sep = "")

    fpathHeight <- paste(data_folder,"/Height.txt",sep = "")

    fpathSed <- paste(data_folder,"/Sediment_Samples.txt",sep = "")

    fpathDischarge <- paste(data_folder,"/Discharge.txt",sep = "")

    fpathSonde <- paste(data_folder,"/Sonde.txt",sep = "")

    fpathADCP <- paste(data_folder,"/ADCP_Data.txt",sep = "")

    fpathADCPBackscatterBeam1 <- paste(data_folder, "/ADCP_Echo_Intensity_Beam_1.txt",sep = "")

    fpathADCPBackscatterBeam2 <- paste(data_folder,"/ADCP_Echo_Intensity_Beam_2.txt",sep = "")
  }
  # Site
  Site <- read.csv(fpathSite)
  Site$Start_date_and_time <- as.POSIXct(Site$Start_date_and_time,format = "%d-%b-%Y %H:%M:%S",tz = "Australia/Brisbane")
  Site$End_date_and_time <- as.POSIXct(Site$End_date_and_time,format = "%d-%b-%Y %H:%M:%S",tz = "Australia/Brisbane")

  # Stage referenced to local datum in Site data frame
  Height <- read.csv(fpathHeight)
  Height$time <- as.POSIXct(Height$time,format = "%d-%b-%Y %H:%M:%S",tz = "Australia/Brisbane")

  # Sediment concentration
  Sediment_Samples <- read.csv(fpathSed)
  Sediment_Samples$time <- as.POSIXct(Sediment_Samples$time,format = "%d-%b-%Y %H:%M:%S",tz = "Australia/Brisbane")

  # Discharge data
  Discharge <- read.csv(fpathDischarge)
  Discharge$time <- as.POSIXct(Discharge$time,format = "%d-%b-%Y %H:%M:%S",tz = "Australia/Brisbane")

  # Sonde data
  Sonde <- read.csv(fpathSonde)
  Sonde$time <- as.POSIXct(Sonde$time,format = "%d-%b-%Y %H:%M:%S",tz = "Australia/Brisbane")

  # ADCP data except backscatter
  ADCP <- read.csv(fpathADCP)
  ADCP$time <- as.POSIXct(ADCP$time,format = "%d-%b-%Y %H:%M:%S",tz = "Australia/Brisbane")

  # Acoustic backscatter from ADCP
  # Beam 1
  Echo_Intensity_Beam_1 <- read.csv(fpathADCPBackscatterBeam1)
  Echo_Intensity_Beam_1$time <- as.POSIXct(Echo_Intensity_Beam_1$time,format = "%d-%b-%Y %H:%M:%S",tz = "Australia/Brisbane")
  # Beam 2
  Echo_Intensity_Beam_2 <- read.csv(fpathADCPBackscatterBeam2)
  Echo_Intensity_Beam_2$time <- as.POSIXct(Echo_Intensity_Beam_2$time,format = "%d-%b-%Y %H:%M:%S",tz = "Australia/Brisbane")

  # Export data -------------------------------------------------------------
  # modify step if one desires to use lower-frequency data
  step <- 1
  indices_for_Input_data <- seq(1,nrow(ADCP),step)
  ind <- indices_for_Input_data

  InputData <- list("Site"=Site,"ADCP"=ADCP[ind,],"Echo_Intensity_Beam_1"=Echo_Intensity_Beam_1[ind,],"Echo_Intensity_Beam_2"=Echo_Intensity_Beam_2[ind,],"Sonde"=Sonde[ind,],"Height" = Height[ind,],"Discharge"=Discharge[ind,],"Sediment_Samples"=Sediment_Samples[ind,])

  return(InputData)
}
