#' Compute attenuation of sound in water given frequency, temperature, and salinity
#'
#' Computes attenuation of sound in water per Ainslie and McColm (1998)
#' @param f frequency of sound (Hz)
#' @param T Water temperature (degrees C)
#' @param S Salinity (PSU)
#' @returns attenuation of sound in water (dB/m), divide by 20*log10(exp(1)) to convert to Nepers/m
#' @examples
#' InputData <- realTimeloads::ExampleData
#' f <- InputData$ADCP$Accoustic_Frequency_kHz*1000
#' cond <-InputData$Sonde$Conductivity_uS_per_cm
#' temp <- InputData$Sonde$Water_Temperature_degC
#' dbar <- InputData$Sonde$Pressure_dbar
#' S <- ctd2sal(cond,temp,dbar)
#' aw <- attenuation_of_sound_by_water(f,temp,S) # dB/m
#' awNp <- attenuation_of_sound_by_water(f,temp,S)/(20*log10(exp(1))) # Np/m
#' @references Ainslie and McColm (1998): https://doi.org/10.1121/1.421258
#'
#' Author modified matlab code from David Schoellhamer
#' @author Daniel Livsey, September (2023), livsey.daniel@@gmail.com, ORCID: 0000-0002-2028-6128
#' @export
#'
attenuation_of_sound_by_water  <- function(f,T,S) {
  # -------------------------------------------------------------------------
  # Description: Function compute sound absorption in sea water based on
  # Ainslie and McColm (1998), as outlined at:
  # http://resource.npl.co.uk/acoustics/techguides/seaabsorption/physics.html#
  # CHECKED CALC AGAINST CALCULATOR ON
  # http://resource.npl.co.uk/acoustics/techguides/seaabsorption (DNL,Dec 8, 2020)
  # Inputs:
  #    f: frequency in Hz
  #    T: temperature in degress C
  #    S: salinity in ppt or psu
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
  f = f*1e-3 # Convert Hz to kHz
  #
  visc=0.00049*f^2*exp(-1*T/27) #viscous absorption in dB/km
  #
  f1=0.78*sqrt(S/35)*exp(T/26) #boric acid
  #
  f2=42*exp(T/17) #Magnesium sulphate
  #
  boric=0.106*f1*f^2/(f1^2+f^2)
  #
  magsul=0.52*(1+T/43)*(S/35)*f2*f^2/(f2^2+f^2)
  #
  aw=(boric+magsul+visc)/1000 #db/m

  return(aw)
}
