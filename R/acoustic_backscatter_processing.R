#' Process acoustic backscatter from hADCP
#'
#' Processes acoustic backscatter from horizontally profiling ADCP (hADCP). Returns attenuation of sound due to water and suspended-sediment. Applies all corrections to acoustic backscatter detailed in the guideline.
#'
#' @param Site Data frame with site, local vertical datum, and ADCP elevation information
#' \describe{
#'   \item{Site_name}{Site name (string)}
#'   \item{Site_number}{Unique site code (string)}
#'   \item{ADCP_elevation_above_bed_m}{Elevation of the ADCP above the bed (m)}
#'   \item{ADCP_elevation_above_gauge_datum_m}{Elevation of the ADCP above local gauge datum (m)}
#'   \item{Distance_of_gauge_datum_below_thalweg_m}{Distance from local gauge datum to lower point in cross-section (m)}
#'   \item{Start_date_and_time}{Installation date of ADCP (time, POSIXct)}
#'   \item{End_date_and_time}{Date if/when ADCP is moved vertically (time, POSIXct)}
#'   \item{Comment}{User comment (string)}
#'   }
#' @param ADCP Data frame with various readings from ADCP
#' \describe{
#'   \item{Site_number}{Unique site code (string)}
#'   \item{time}{Date and time (time, POSIXct)}
#'   \item{Ensemble}{Measurment ensemble number (integer)}
#'   \item{Accoustic_Frequency_kHz}{Acoustic frequency of ADCP (kHz)}
#'   \item{Transducer_radius_m}{Radius of ADCP transducer (m)}
#'   \item{Beam_angle_degrees}{Angle of beam relative to normal (degrees)}
#'   \item{Beam_aspect_ratio}{Ratio of beam radius to beam length (-)}
#'   \item{Range_to_bed_of_acoustic_beams_m}{Normal range to bed, optional (m)}
#'   \item{Range_to_water_surface_of_acoustic_beams_m}{Normal range to water surface, optional (m)}
#'   \item{Number_of_Cells}{Number of measurement cells along beam (integer)}
#'   \item{Bin_Size_m}{Cell width measured normal to ADCP (m)}
#'   \item{Blanking_distance_m}{Blanking distance measured normal to ADCP (m)}
#'   \item{Instrument_serial_number}{Serial number of ADCP instrument (string)}
#'   \item{CPU_serial_number}{Serial number of ADCP CPU (string)}
#'   \item{Ambient_Noise_Level_Beam_1_Counts}{Ambient noise level for beam 1, optional (counts)}
#'   \item{Ambient_Noise_Level_Beam_2_Counts}{Ambient noise level for beam 2, optional (counts)}
#'   \item{Distance_to_Bin_1_mid_point_m}{Reported distance normal to ADCP to midpoint of bin/cell (m)}
#'   \item{Distance_to_surface_m}{Reported depth of ADCP from vertical beam, optional (m)}
#'   \item{Speed_of_sound_m_per_s}{Speed of sound used by ADCP in the field (m/s)}
#'   \item{Temperature_degC}{Temperature recorded by ADCP (degrees C)}
#'   \item{Pressure_dbar}{Pressure recorded by ADCP (dBar)}
#'   \item{Salinity_PSU}{Salinity in PSU recorded or assumed in ADCP data file, optional (PSU)}
#'   \item{Distance_to_surface_m}{Distance to water surface reported by vertical beam of ADCP (m)}
#'   \item{Power_supply_voltage}{Power to ADCP (V)}
#'   }
#'@param Height Data frame with timeseries of river height
#' \describe{
#'   \item{time}{Date and time (time, POSIXct)}
#'   \item{Height_m}{Water surface elevation above gauge datum (m)}
#'   \item{Site_number}{Unique site code (string)}
#'   }
#' @param Sonde Data frame with timeseries of conductivity, temperature, and depth from sonde
#' \describe{
#'   \item{time}{Date and time (time, POSIXct)}
#'   \item{Water_Temperature_degC}{Temperature (degrees C)}
#'   \item{Conductivity_uS_per_cm}{Conductivity (microS/cm)}
#'   \item{Pressure_dbar}{Pressure (dbar)}
#'   \item{Site_number}{Unique site code (string)}
#'   }
#' @param Echo_Intensity_Beam_1 Data frame of acoustic backscatter measurements from beam 2
#' \describe{
#'   \item{Site_number}{Unique site code (string)}
#'   \item{time}{Date and time (time, POSIXct)}
#'   \item{Echo_Intensity_Counts_cell_n}{Acoustic backscatter in nth cell (counts)}
#'   }
#' @param Echo_Intensity_Beam_2 Data frame of acoustic backscatter measurements from beam 2
#' \describe{
#'   \item{Site_number}{Unique site code (string)}
#'   \item{time}{Date and time (time, POSIXct)}
#'   \item{Echo_Intensity_Counts_cell_n}{Acoustic backscatter in nth cell (counts)}
#'   }
#' @param Instrument_Noise_Level Estimate of noise level, recommended if ambient noise level is not recorded (counts)
#' @param Include_Rayleigh Logical to include data within Rayleigh Distance for processing of acoustic backsactter
#' @param Include_near_field_correction Logical to include near-field correction of Downing et al (1995)
#' @returns List with processed data, all variable names and units are written-out in list items, see GUIDLINE for details of each variable
#' @examples
#' \dontrun{
#' InputData <- realTimeloads::ExampleData
#' Site <- InputData$Site
#' ADCP <- InputData$ADCP
#' Height <- InputData$Height
#' Sonde <- InputData$Sonde
#' EIa <- InputData$Echo_Intensity
#' # example code assumes backscatter is equal across beams
#' EIb <- InputData$Echo_Intensity
#' Output <- acoustic_backscatter_processing(Site,ADCP,Height,Sonde,EIa,EIb)
#' }
#' @references
#' Livsey, D.N. (in review). National Industry Guidelines for hydrometric monitoring–Part 12: Application of acoustic Doppler velocity meters to measure suspended-sediment load. Bureau of Meteorology. Melbourne, Australia.
#'
#' @author Daniel Livsey (2023) ORCID: 0000-0002-2028-6128
#' @export
#'
acoustic_backscatter_processing  <-  function(Site,ADCP,Height,Sonde,Echo_Intensity_Beam_1,Echo_Intensity_Beam_2,Instrument_Noise_Level = NULL,Include_Rayleigh=FALSE,Include_near_field_correction=TRUE) {

# Site datum and ADCP elevation above bed
mab <- Site$ADCP_elevation_above_bed_m # ADCP height above bed (m,up positive)
Vertical_Distance_From_Gauge_to_Thalweg = Site$Distance_of_gauge_datum_below_thalweg_m  # distance between stage datum and thalweg (m,up positive)
ADCP_m_above_datum <- Site$ADCP_elevation_above_gauge_datum_m # distance from ADCP elevation to stage datum (m, up positive)

# ADCP instrument characteristics and settings
f <- ADCP$Accoustic_Frequency_kHz*1000 # kHz to Hz, from manufacturer, note a 600 KhZ ADCP has frequency close to 600 kHz but the actual value should be read from the instrument
at <- ADCP$Transducer_radius_m # m, from manufacturer
binSize <- ADCP$Bin_Size_m # m, set by user
AR <- ADCP$Beam_aspect_ratio # ratio of the distance of ADCP beam from the ADCP to the distance to the nearest boundary, accounts for conical shape of beam as it spreads-out get value from ADCP manual

t <- ADCP$time
nt <- length(t)
ncells <- ADCP$Number_of_Cells[1]

# Only process data if cell settings are held constant
test_constant_cell_settings <- length(unique(ADCP$Blanking_distance_m))+length(unique(ADCP$Number_of_Cells))+length(unique(ADCP$Bin_Size_m))

if (test_constant_cell_settings != 3) {
  stop('ERROR! Code assumes constant cell settings for processing')
  Output <- 'ERROR! Code assumes constant cell settings for processing'
  }

if (test_constant_cell_settings == 3) {

  # Depth of ADCP and total depth
  H <- Height$Height_m - Vertical_Distance_From_Gauge_to_Thalweg # Total depth of channel at deepest point (m)
  H <- linear_interpolation_with_time_limit(Height$time,H,t,Inf)$x_interpolated
  D <- H - ADCP_m_above_datum # Depth of the ADCP, one can use ADCP pressure (dBar~= 1 m water depth) or vertical distance water surface from vertical beam (if available)

  # prefer vertical beam depth from ADCP if available
  if (sum(colnames(ADCP)=="Distance_to_surface_m")) {
    Dv <- ADCP$Distance_to_surface_m
    Dv[Dv>1000] <- NA # on RDI ADCP Depth is ~ 6000 m when reading is invalid
    if (sum(is.finite(Dv))>2) {
    D <- imputeTS::na_ma(Dv)
    }
  }

# Speed of sound in water (c, m/s)
T <- ADCP$Temperature_degC # as record by ADCP

# Conductivity to Salinity (PSU), if Conductivity is specific conductivity T in ctd2sal must be set to 25, T and D MUST be the T and D where the conductivity was measured, this may or may not be at ADCP position (e.g, in flow cell)
# if sonde is in flow cell, set pressure to 0.
S <- ctd2sal(Sonde$Conductivity_uS_per_cm,Sonde$Temperature_degC,Sonde$Pressure_dBar)
# allow  linear interpolation over all time gaps
S <- linear_interpolation_with_time_limit(Sonde$time,S,t,Inf)
S <- S$x_interpolated

# Infill missing salinity data
if (sum(is.finite(S)) >= 2) {
  S <- imputeTS::na_ma(S)
}
if (sum(is.finite(S)) < 2) {
  if (is.element('Salinity_PSU',colnames(ADCP))) {
    S <- ADCP$Salinity_PSU
    warning('No salinity data from sonde, setting salinity to user-programmed value')
  }
  if (!is.element('Salinity_PSU',colnames(ADCP))) {
  S[] <- 0
  warning('No salinity data from sonde or ADCP, assuming salinity of 0 PSU')
  }
}

# S,T, and D is that of the ADCP not necessarily that of the Sonde (e.g., if the sonde is in a flow cell)
# ADCP may use specified speed of sound to estimate pulse length (tau), user should check if ADCP uses constant speed of sound or is estimating c w/ time
co <- ADCP$Speed_of_sound_m_per_s
Dc <- D
Dc[D<0] <-0
c <- speed_of_sound(S,T,Dc) # actual speed of sound in water (m/s)

# Rayleigh distance (m)
lamba <- c/(f) # wavelength of sound (m)
k = (2*pi)/lamba # acoustic wave number (1/m)
Rayleigh_Distance = (pi*at^2)/lamba # Rayleigh distance (m)


# compute echo intensity scalar (kc) per RDI (deciBels per count)
# T is temperature of the ADCP (deg C) NOT necessarily that of the Sonde
kc <- 127.3/(T+273)

# attenuation of sound in water (aw, dB/m)
# T is temperature of the ADCP (deg C) NOT necessarily that of the Sonde
aw <- attenuation_of_sound_by_water(f,T,S)

# Compute maximum profiling distance considering beam spreading
# assumes the ADCP is horizontal
Rmax <- D*AR
Rbed <- mab*AR # range to bed
Rmax[Rmax>Rbed] <- Rbed[Rmax>Rbed]

# If Range_to_bed_of_acoustic_beams_m and Range_to_water_surface_of_acoustic_beams_m provided use these variables to define Rmax
# Range_to_bed_of_acoustic_beams_m is needed for sites with irregular channel bottom where beam may hit bed at bed elevation prior to thalweg elevation, pointed out by Stephen Wallace
# Range_to_water_surface_of_acoustic_beams_m can be computed from aspect ratio and depth reported by vertical beam available on some ADCPs
if (is.element('Range_to_bed_of_acoustic_beams_m',colnames(ADCP)) & is.element('Range_to_water_surface_of_acoustic_beams_m',colnames(ADCP))) {
  Rmax <- apply(cbind(ADCP$Range_to_bed_of_acoustic_beams_m,ADCP$Range_to_water_surface_of_acoustic_beams_m),1, FUN = min, na.rm = TRUE)
}

# Correct counts for changes in power supply voltage
V <- ADCP$Power_supply_voltage
Pc_counts = 20*log10(V/12)/kc # see  eq. 10 in guideline

# Check for Ambient Noise Level data, set to NA if no data provided
if (!is.element('Ambient_Noise_Level_Beam_1_Counts',colnames(ADCP))) {
  ADCP$Ambient_Noise_Level_Beam_1_Counts <- NA
}
if (!is.element('Ambient_Noise_Level_Beam_2_Counts',colnames(ADCP))) {
  ADCP$Ambient_Noise_Level_Beam_2_Counts <- NA
}

# determine if ambient noise level was measured
m_amb_1 <- is.finite(ADCP$Ambient_Noise_Level_Beam_1_Counts)
m_amb_2 <- is.finite(ADCP$Ambient_Noise_Level_Beam_2_Counts)
# where ambient noise is measured, one needs correction for changes in power supply
# where NL is not estimated impute with mean, note if one or more ambient noise levels are provided noise level is imputed to mean of those observations NOT Instrument_Noise_Level
if (sum(m_amb_1)>0) {
  NLa <- rep(NA,nrow(ADCP))
  NLa[m_amb_1] <- ADCP$Ambient_Noise_Level_Beam_1_Counts[m_amb_1] + Pc_counts[m_amb_1]
  NLa[!m_amb_1] <- mean(NLa[m_amb_1]) + Pc_counts[!m_amb_1] # infill missing data with mean ambient noise level and apply correction for changes in power supply, adding power supply correction b/c mean(NLa[m_amb_1]) is best estimate of noise level corrected to nominal voltage
}
if (sum(m_amb_2)>0) {
  NLb <- rep(NA,nrow(ADCP))
  NLb[m_amb_2] <- ADCP$Ambient_Noise_Level_Beam_2_Counts[m_amb_2] + Pc_counts[m_amb_2]
  NLb[!m_amb_2] <- mean(NLb[m_amb_2]) + Pc_counts[!m_amb_2] # infill missing data with mean ambient noise level and apply correction for changes in power supply, adding power supply correction b/c mean(NLa[m_amb_1]) is best estimate of noise level corrected to nominal voltage
}

# Set NLa to estimate of instrument noise if provided
if (sum(m_amb_1)==0&!is.null(Instrument_Noise_Level)) {
  NLa <- Instrument_Noise_Level*rep(1,nrow(ADCP)) + Pc_counts # assuming instrument noise level measured at nominal power supply in office/or lab a correction for changes in power supply in the field would be needed
}
# Set NLb to estimate of instrument noise if provided
if (sum(m_amb_2)==0&!is.null(Instrument_Noise_Level)) {
  NLb <- Instrument_Noise_Level*rep(1,nrow(ADCP)) + Pc_counts # assuming instrument noise level measured at nominal power supply in office/or lab a correction for changes in power supply in the field would be needed
}

# Case where no Ambient_Noise_Level or Instrument Noise Level is provided, warn user
if (sum(m_amb_1)==0 & is.null(Instrument_Noise_Level)) {
  # Noise level (counts) for each beam will be 0 if no Ambient or Instrument Noise Level data are provided
  # adding correction for power supply changes so that changes in power supply do not introduce artificial change in EIsnr2a and EIsnr2b below
  NLa <- rep(0,nrow(ADCP)) + Pc_counts
  # Ensure NLa and NLb do not go below zero when no ambient or instrument noise level data are given and Pc_counts is negative
  if (min(Pc_counts,na.rm=TRUE)<0) {
    NLa <-  NLa - min(Pc_counts,na.rm=TRUE)
  }
  warning('No Ambient_Noise_Level_Beam_1_Counts recorded and Instrument_Noise_Level is NULL, estimate of noise level should be provided in Instrument_Noise_Level, see Haught et al (2017): https://doi.org/10.1002/2016WR019695')
}
if (sum(m_amb_2)==0 & is.null(Instrument_Noise_Level)) {
  # Noise level (counts) for each beam will be 0 if no Ambient or Instrument Noise Level data are provided
  # adding correction for power supply changes so that changes in power supply do not introduce artificial change in EIsnr2a and EIsnr2b below
  NLb <- rep(0,nrow(ADCP)) + Pc_counts
  # Ensure NLa and NLb do not go below zero when no ambient or instrument noise level data are given and Pc_counts is negative
  if (min(Pc_counts,na.rm=TRUE)<0) {
    NLb <-  NLb - min(Pc_counts,na.rm=TRUE)
  }
  warning('No Ambient_Noise_Level_Beam_2_Counts recorded and Instrument_Noise_Level is NULL, estimate of noise level should be provided in Instrument_Noise_Level, see Haught et al (2017): https://doi.org/10.1002/2016WR019695')
}

# Echo intensity when SNR (linear units) equals 2
# One should not use velocity or backscatter for sediment when SNR ~=2
EIsnr2a = 10*log10(10^(0.1*kc*NLa)+1)/kc + NLa
EIsnr2b = 10*log10(10^(0.1*kc*NLb)+1)/kc + NLb

# Get counts in matrices
EIa<-as.matrix(Echo_Intensity_Beam_1[grepl("Counts",colnames(Echo_Intensity_Beam_1))])+Pc_counts # a is beam 1

EIb<-as.matrix(Echo_Intensity_Beam_2[grepl("Counts",colnames(Echo_Intensity_Beam_2))])+Pc_counts # b is beam 2


# Compute MBsnr, MBsnr are counts converted to decibels and corrected to noise floor, MBsnr is the signal-to-noise ratio in log units but one needs to QA/QC data using SNR in linear units
MBsnr_a <- 10*log10(10^(kc*(EIa-NLa)/10)-1)
MBsnr_b <- 10*log10(10^(kc*(EIb-NLb)/10)-1)
colnames(MBsnr_a) = gsub("Echo_Intensity", "MBsnr", colnames(MBsnr_a))
colnames(MBsnr_b) = gsub("Echo_Intensity", "MBsnr", colnames(MBsnr_b))

# NA or -Inf can be returned if user Instrument Noise Level is higher than recorded EI
MBsnr_a[abs(MBsnr_a)==Inf] <- NA
MBsnr_a[is.nan(MBsnr_a)] <- NA
MBsnr_b[abs(MBsnr_a)==Inf] <- NA
MBsnr_b[is.nan(MBsnr_a)] <- NA

# Compute range measured along beam for each ping to cell center
# And find indices of cells w/ echo intensity > SNR 2
# Bin size and blanking distance are measured perpendicular to the face of the ADCP r is measured along the beam, correct using Beam_angle_degrees
cosr <- cos(ADCP$Beam_angle_degrees*pi/180)
# actual bin 1 start (r1, m) and cell size (dr, m) (as measured along the beam), "o" indicates value used by instrument in the field
r1o = ADCP$Distance_to_Bin_1_mid_point_m/cosr
dro <- ADCP$Bin_Size_m/cosr
r1 = r1o*(c/co) # from mod of eq. 7 of RDI ADCP Principles of Operation
dr <- dro*(c/co)
rend = r1+dr*ncells
# Pre-allocate matrices
r =  matrix(NA, nrow = nt, ncol = ncells)
igooda =  matrix(FALSE, nrow = nt, ncol = ncells)
igoodb =  matrix(FALSE, nrow = nt, ncol = ncells)

for(i in 1:length(t)) {
  # print(i)
  r[i,] = seq(r1[i],r1[i]+dr[i]*(ncells-1),length.out = ncells)
  # logical to determine cell should be used to compute mean backscatter from WCB or SCB.
  # (r[i,]+dr[i]) dr[i] is added below to compute range to end of cells

  # added option to include or exclude backsactter data within Rayleigh distance

  # backsactter data w /in Rayleigh_Distance included in processing
  if (Include_Rayleigh) {
  # (r[i,]+dr[i]) dr[i] is added to compute range to end of cells
  igooda[i,] = EIa[i,]>EIsnr2a[i] & (r[i,]+dr[i]) < Rmax[i]
  igoodb[i,] = EIb[i,]>EIsnr2b[i] & (r[i,]+dr[i]) < Rmax[i]
  }
  # backsactter data w/in Rayleigh_Distance excluded in processing
  # (r[i,]+dr[i]) dr[i] is added to compute range to end of cells
  if (!Include_Rayleigh) {
    igooda[i,] = EIa[i,]>EIsnr2a[i] & (r[i,]+dr[i]) > Rayleigh_Distance[i] & (r[i,]+dr[i]) < Rmax[i]
    igoodb[i,] = EIb[i,]>EIsnr2b[i] & (r[i,]+dr[i]) > Rayleigh_Distance[i] & (r[i,]+dr[i]) < Rmax[i]
  }
}

# ensure no NA in logical matrices
igooda[is.na(igooda)] <- FALSE
igoodb[is.na(igoodb)] <- FALSE
# number of valid cells per ping
igna <- rowSums(igooda)
# number of valid cells per ping
ignb <- rowSums(igoodb)

# Near-field correction (dimensionless)
# added option to include or exclude Near-field correction
if (!Include_near_field_correction & Include_Rayleigh) {
  warning('Near-field correction likely needed when using data within Rayleigh Distace') }

psi <- near_field_correction(f,c,r,at)

if (!Include_near_field_correction) {
  psi[] <- 1 #
}

# Correction  for beam spreading and attenuation due to water (dB)
corr <- 2*r*aw + 20*log10(r*psi)
# If beyond the Rayleigh distance one can omit the near-field correction (psi)
# corr <- 2*r*aw* + 20*log10(r)

# Compute WCB (dB)
WCBa = MBsnr_a+corr
WCBb = MBsnr_b+corr

# Compute sediment attenuation coefficient (dB/m)
# Preallocate matrices
as_a <- rep(NA,length(t))
as_b <- rep(NA,length(t))
as_a_pvalue <- rep(NA,length(t))
as_b_pvalue <- rep(NA,length(t))
SCBa =  matrix(NA, nrow = nt, ncol = ncells)
SCBb =  matrix(NA, nrow = nt, ncol = ncells)
# mean SCB
muSCBa =  rep(NA,length(t))
muSCBb =  rep(NA,length(t))
#(i in 1:length(t))
for (i in 1:length(t)) {
#print(i)
# compute slope of WCB(r) for beam 1
if (igna[i]>3) {
x <- as.matrix(r[i,igooda[i,]])
y <- as.matrix(WCBa[i,igooda[i,]])
df <- data.frame(x,y)
# fit linear model
fit <- lm(y ~ x, data=df)
# get p-value of f statistic
f <- summary(fit)$fstatistic
as_a[i] = fit$coefficients[2]*-0.5
p <- pf(f[1],f[2],f[3],lower.tail=F)
attributes(p) <- NULL
as_a_pvalue[i] <- p
# Compute SCB only if as is positive and statistically significant
if (as_a[i]>0 & as_a_pvalue[i]<=0.05) {
SCBa[i,] = WCBa[i,]+2*r[i,]*as_a[i]
muSCBa[i] = mean(SCBa[i,igooda[i,]],na.rm = TRUE)
# view summary of linear model
#summary(fit)
# check code
#yp <- predict(fit,df)
#plot(r[i,],MBsnr_a[i,],ylim = c(min(MBsnr_a[i,])-1,max(MBsnr_a[i,])+10))
#points(r[i,],WCBa[i,],col='red')
#points(x, y, col='blue')
#lines(x,yp,col='green')
#points(r[i,],SCBa[i,],col='green')
}

if (as_a[i]<0 | as_a_pvalue[i]>0.05) {
  muSCBa[i] = mean(WCBa[i,igooda[i,]],na.rm = TRUE)
}

}

# compute slope of WCB(r) for beam 2
if (ignb[i]>3) {
  x <- as.matrix(r[i,igoodb[i,]])
  y <- as.matrix(WCBb[i,igoodb[i,]])
  df <- data.frame(x,y)
  # fit linear model
  fit <- lm(y ~ x, data=df)
  # get p-value of f statistic
  f <- summary(fit)$fstatistic
  as_b[i] = fit$coefficients[2]*-0.5
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  as_b_pvalue[i] <- p
  # Compute SCB only if as is positive and statistically significant
  if (as_b[i]>0 & as_b_pvalue[i]<=0.05) {
    SCBb[i,] = WCBb[i,]+2*r[i,]*as_b[i]
    muSCBb[i] = mean(SCBb[i,igoodb[i,]],na.rm = TRUE)

    # view summary of linear model
    #summary(fit)
    # check code
    #yp <- predict(fit,df)
    #plot(r[i,],MBsnr_b[i,],ylim = c(min(MBsnr_b[i,])-1,max(MBsnr_b[i,])+10))
    #points(r[i,],WCBb[i,],col='red')
    #points(x, y, col='blue')
    #lines(x,yp,col='green')
    #points(r[i,],SCBb[i,],col='green')
  }

  if (as_b[i]<0 | as_b_pvalue[i]>0.05) {
    muSCBb[i] = mean(WCBb[i,igoodb[i,]],na.rm = TRUE)
  }

}
}

# Compute mean dB from both beams
dB <- rowMeans(cbind(muSCBb,muSCBa),na.rm = TRUE)
# Double Checked code results by comparing dB to Matlab dB

# Save processed data
time <- t
Mean_sediment_corrected_backscatter_dB <- dB
Processed_Backscatter <- data.frame(time,Mean_sediment_corrected_backscatter_dB)

Output <- list("time" = t,"mean_sediment_corrected_backscatter_beam_1_and_beam_2_dB" = dB,"mean_sediment_corrected_backscatter_beam_1_dB" = muSCBa,"mean_sediment_corrected_backscatter_beam_2_dB" = muSCBb,"range_along_beam_m" = r,"echo_intensity_beam_1_counts"= EIa,"echo_intensity_beam_2_counts"= EIb,"noise_level_counts_beam_1" = NLa,"noise_level_counts_beam_2" = NLb,"measured_backscatter_in_SNR_beam_1_dB" = MBsnr_a,"measured_backscatter_in_SNR_beam_2_dB" = MBsnr_b,"water_corrected_backscatter_beam_1_dB" = WCBa,"water_corrected_backscatter_beam_2_dB" = WCBb,"sediment_corrected_backscatter_beam_1_dB" = SCBa,"sediment_corrected_backscatter_beam_2_dB" = SCBb,"attenuation_of_sound_due_to_water_dB_per_m" = aw,"attenuation_of_sound_due_to_sediment_beam_1_dB_per_m" = as_a,"attenuation_of_sound_due_to_sediment_beam_2_dB_per_m" = as_b,"kc_dB_per_counts" = kc,"speed_of_sound_m_per_sec" = c)

#write.csv(Processed_Backscatter, "Data/Processed_Backscatter.csv", row.names=FALSE)

}

return(Output)

}
