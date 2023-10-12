#' Compute salinity (PSU) from conductivity, water temperature, and depth
#'
#' Computes salinity from conductivity, water temperature, and depth.
#' @param cond Conductance (microS/cm)
#' @param temp Water temperature (degrees C)
#' @param dbar Pressure (dBar) or water depth (m)
#' @returns Salinity in PSU
#' @section Warning:
#' If specific conductivity is returned from the sonde, the temperature at which specific conductivity is computed should be utilized
#' @examples
#' Sonde <- realTimeloads::ExampleData$Sonde
#' S <- ctd2sal(Sonde$Conductivity_uS_per_cm,Sonde$Water_Temperature_degC,Sonde$Pressure_dbar)
#' @references
#' Fofonoff, N. P., & Millard Jr, R. C. (1983). Algorithms for the computation of fundamental properties of seawater.
#'
#' Chen, C. T. A., & Millero, F. J. (1986). Thermodynamic properties for natural waters covering only the limnological range 1. Limnology and Oceanography, 31(3), 657-662.
#'
#' Hill, K., Dauphinee, T., & Woods, D. (1986). The extension of the Practical Salinity Scale 1978 to low salinities. IEEE Journal of Oceanic Engineering, 11(1), 109-112.
#'
#' Author modified Matlab code from David Schoellhamer
#' @author Daniel Livsey (2023) ORCID: 0000-0002-2028-6128
#' @export
#'
ctd2sal  <-  function(cond,temp,dbar) {
# Returns salinity from cond (microS/cm), temp (deg C),
# and bar (dbar),
# 1 dbar ~= 1 m depth,note affect of depth < 100 m is minimal
# if cond is specific conductivity T must be set to 25 deg
# From UNESCO technical paper 44, 1983, with revisions to fresh
# water from Chen and Millero, 1986.
#
# Low salinity correction from Hill and others (1986)
#
# modified by David Schoellhamer, 6/22/2001
# Translated to R by Daniel Livsey, 13/07/2023

# Without low salinity correction and commented out R:
  # Test Values          Cond       Temp     Bar        Sal
#                 41.4691359462    15       0       35.0000
#                 49.7629631354    20      200      37.2456
#                 26.9549383650     5      150      27.9953
# ******************************************************************

  ##### define constants for salinity calculation
sala0 <- 0.0080
sala1 <- -0.1692
sala2 <- 25.3851
sala3 <- 14.0941
sala4 <- -7.0261
sala5 <- 2.7081
salb0 <- 0.0005
salb1 <- -0.0056
salb2 <- -0.0066
salb3 <- -0.0375
salb4 <- 0.0636
salb5 <- -0.0144
salk <- 0.0162
salc0 <- 0.6766097
salc1 <- 2.00564e-2
salc2 <- 1.104259e-4
salc3 <- -6.9698e-7
salc4 <- 1.0031e-9
sale1 <- 2.070e-5
sale2 <- -6.370e-10
sale3 <- 3.989e-15
sald1 <- 3.426e-2
sald2 <- 4.464e-4
sald3 <- 4.215e-1
sald4 <- -3.107e-3

####convert micro s to milli s
cond  <-  cond/1000

####default to 25 deg and 0 dbar
if (length(temp)==0) {
temp  <-  25*rep(1, length(c))
}

if (length(dbar)==0) {
dbar <- rep(0, length(c))
}

##### temperature correction to conductivity
rt  <-  salc0+salc1*temp+salc2*temp^2+salc3*temp^3+salc4*temp^4

##### salinity correction to conductivity (41.4691359462 is cond at (SAL,TEMP,PRESS) 35,20,0)
#R  <-  cond/41.4691359462
#note: Jon Burau uses 42.922, Hill and others use 42.914
R  <-  cond/42.922

##### pressure correction to conductivity
Rp  <-  1+(dbar*(sale1+sale2*dbar+sale3*dbar^2))/ (1+sald1*temp+sald2*temp^2+R*(sald3+sald4*temp))
Rt  <-  R/(Rp*rt)

##### calculate the salinity
#use sqrt to avoid complex numbers
root <- sqrt(Rt)
#delS  <-  ((temp-15)/(1+salk*(temp-15)))*(salb0+salb1*Rt^(1/2)+salb2*Rt+salb3*Rt^(3/2)+salb4*Rt^2+salb5*Rt^(5/2))
#sal  <-  sala0+sala1*Rt^(1/2)+sala2*Rt+sala3*Rt^(3/2)+sala4*Rt^2+sala5*Rt^(5/2)+delS

delS  <-  ((temp-15)/(1+salk*(temp-15)))*(salb0+salb1*root+salb2*Rt+salb3*root^3+salb4*Rt^2+salb5*root^5)
sal  <-  sala0+sala1*root+sala2*Rt+sala3*root^3+sala4*Rt^2+sala5*root^5+delS

#low salinity correction
x <- 400*Rt
y <- 100*Rt
rooty <- sqrt(y)
fT <- (temp-15)/(1+salk*(temp-15))
sal <- sal-sala0/(1+1.5*x+x^2)-salb0*fT/(1+rooty+y+rooty^3)
sal <- round(sal,2)

}
