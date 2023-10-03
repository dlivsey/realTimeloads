#' Compute speed of sound in water given salinity, temperature, and depth
#'
#' Computes speed of sound in water per Del grosso (1974)
#' @param S Salinity (PSU)
#' @param T Water temperature (degrees C)
#' @param D Water depth (m) or pressure (dBar)
#' @returns Speed of sound in water (m/s)
#' @examples
#' InputData <- realTimeloads::ExampleData
#' Sonde<- InputData$Sonde
#' S <- ctd2sal(Sonde$Conductivity_uS_per_cm,Sonde$Water_Temperature_degC,Sonde$Pressure_dbar)
#' c <- speed_of_sound(S,Sonde$Water_Temperature_degC,Sonde$Pressure_dbar)
#' @references Del grosso (1974): https://doi.org/10.1121/1.1903388
#'
#' Author modified matlab code from David Schoellhamer
#' @author Daniel Livsey, September (2023), livsey.daniel@@gmail.com, ORCID: 0000-0002-2028-6128
#' @export
#'
speed_of_sound <- function(S,T,D) {
  # S salinity (ppt or psu)
  # T water temperature (deg C)
  # D water depth (m)
  # Del grosso uses pressure in kg/cm^2. To get to this from dbars one must
  # divide by "g". From the UNESCO algorithms (referring to ANON (1970)
  # BULLETIN GEODESIQUE) we have this formula for g as a function of latitude
  # and pressure. We set latitude to 45 degrees for convenience!
  #   Del Grosso, "A New Equation for the speed of sound in Natural
  #            Waters", J. Acoust. Soc. Am. 56#4 (1974).

  #D <- 1
  #S <- 30
  #T <- 15

  XX <- sin(45*pi/180)
  GR  <-  9.780318*(1+(5.2788E-3+2.36E-5*XX)*XX) + 1.092E-6*D
  P <- D/GR;

  C000  <-  1402.392;
  DCT  <-  (0.501109398873e1-(0.550946843172e-1 - 0.221535969240e-3*T)*T)*T
  DCS  <-  (0.132952290781e1 + 0.128955756844e-3*S)*S
  DCP  <-  (0.156059257041e0 + (0.244998688441e-4 - 0.883392332513e-8*P)*P)*P
  DCSTP  <-  -0.127562783426e-1*T*S + 0.635191613389e-2*T*P +0.265484716608e-7*T*T*P*P- 0.159349479045e-5*T*P*P+0.522116437235e-9*T*P*P*P - 0.438031096213e-6*T*T*T*P
  DCSTP <- DCSTP - 0.161674495909e-8*S*S*P*P + 0.968403156410e-4*T*T*S+0.485639620015e-5*T*S*S*P - 0.340597039004e-3*T*S*P
  ssp <-  C000 + DCT + DCS + DCP + DCSTP

  return(ssp)
}
