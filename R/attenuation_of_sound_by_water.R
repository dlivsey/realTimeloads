#' Compute attenuation of sound in water given frequency, temperature, and salinity
#'
#' Computes attenuation of sound in water per Ainslie and McColm (1998)
#' @param freq frequency of sound (Hz)
#' @param temp Water temperature (degrees C)
#' @param sal Salinity (PSU)
#' @returns attenuation of sound in water (dB/m), divide by 20*log10(exp(1)) to convert to Nepers/m
#' @examples
#' InputData <- realTimeloads::ExampleData
#' freq <- InputData$ADCP$Accoustic_Frequency_kHz*1000
#' cond <-InputData$Sonde$Conductivity_uS_per_cm
#' temp <- InputData$Sonde$Water_Temperature_degC
#' dbar <- InputData$Sonde$Pressure_dbar
#' sal <- ctd2sal(cond,temp,dbar)
#' aw <- attenuation_of_sound_by_water(freq,temp,sal) # dB/m
#' awNp <- attenuation_of_sound_by_water(freq,temp,sal)/(20*log10(exp(1))) # Np/m
#' @references
#' Ainslie, M. A., & McColm, J. G. (1998). A simplified formula for viscous and chemical absorption in sea water. The Journal of the Acoustical Society of America, 103(3), 1671-1672.
#'
#' Author modified Matlab code from David Schoellhamer
#' @author Daniel Livsey (2023) ORCID: 0000-0002-2028-6128
#' @export
#'
attenuation_of_sound_by_water  <- function(freq,temp,sal) {
  # -------------------------------------------------------------------------
  # Description: Function compute sound absorption in sea water based on
  # Ainslie and McColm (1998), as outlined at:
  # http://resource.npl.co.uk/acoustics/techguides/seaabsorption/physics.html#
  # CHECKED CALC AGAINST CALCULATOR ON
  # http://resource.npl.co.uk/acoustics/techguides/seaabsorption (DNL,Dec 8, 2020)
  # Inputs:
  #    freq: frequency in Hz
  #    temp: temperature in degress C
  #    sal: salinity in ppt or psu
  #
  # Outputs:
  #   aw: absorption in dB/m (For conversion to Nepers/m use 1 Np =
  #   (20*log10(exp(1))) dB
  # See: Appendix 1: Misra, D. K. (2012). Radio-frequency and microwave communication circuits: analysis and design. John Wiley & Sons.
  #
  # Assumptions:
  #    Depth ~= 0 (i.e. near the surface)
  #    pH=8
  #
  #
  #-------------------------------------------------------------------------
  #
  freq = freq*1e-3 # Convert Hz to kHz
  #
  visc=0.00049*freq^2*exp(-1*temp/27) #viscous absorption in dB/km
  #
  f1=0.78*sqrt(sal/35)*exp(temp/26) #boric acid
  #
  f2=42*exp(temp/17) #Magnesium sulphate
  #
  boric=0.106*f1*freq^2/(f1^2+freq^2)
  #
  magsul=0.52*(1+temp/43)*(sal/35)*f2*freq^2/(f2^2+freq^2)
  #
  aw=(boric+magsul+visc)/1000 #db/m

  return(aw)
}
