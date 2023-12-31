% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impute_data.R
\name{impute_data}
\alias{impute_data}
\title{Returns x with gaps imputed using ARIMA and Decision Trees, optional uncertainty estimation using Monte Carlo resampling}
\usage{
impute_data(
  time,
  x,
  Xreg = NULL,
  ti = NULL,
  hfit = NULL,
  harmonic = FALSE,
  only_use_Xreg = FALSE,
  MC = 1,
  ptrain = 1
)
}
\arguments{
\item{time}{time for x (time, POSIXct)}

\item{x}{any quantity (double)}

\item{Xreg}{additional predictors for decision tree, required if harmonic is FALSE (rows = time, or if given, ti)}

\item{ti}{time vector for interpolation (time, POSIXct)}

\item{hfit}{model object from TideHarmonics::ftide}

\item{harmonic}{logical if x exhibits tidal or diurnal variability}

\item{only_use_Xreg}{logical for using Xreg only in decision tree}

\item{MC}{number of Monte Carlo simulations for uncertainty estimation}

\item{ptrain}{proportion of data used for training and testing model}
}
\value{
list with x imputed at time or ti, if given. Uncertainty estimated from Monte Carlo simulations
}
\description{
Returns x with gaps imputed using ARIMA and Decision Trees with option to use  harmonic model as predictors for x in decision tree algorithm. Uncertainty on imputed data is estimated using using Monte Carlo (MC) resampling adapting methods of Rustomji and Wilkinson (2008)
}
\note{
If MC == 1, uncertainty is not evaluated. If ptrain == 1, uncertainty and validation accuracy are not computed
}
\examples{
# Impute non-tidal data
time <- realTimeloads::ExampleData$Sediment_Samples$time
xo <- realTimeloads::ExampleData$Sediment_Samples$SSCxs_mg_per_liter
Q <- realTimeloads::ExampleData$Discharge$Discharge_m_cubed_per_s
idata <- sample(1:length(xo),round(length(xo)*0.5),replace=FALSE)
x <- rep(NA,length(xo))
x[idata] <- xo[idata] # simulated samples
flow_concentrtion_ratio <- imputeTS::na_interpolation(Q/x)
Xreg <- cbind(Q,flow_concentrtion_ratio)
Output <- impute_data(time,x,Xreg,MC = 10,ptrain = 0.8)

# Impute tidal data
time <-TideHarmonics::Portland$DateTime[1:(24*90)]
xo <-TideHarmonics::Portland$SeaLevel[1:(24*90)]
idata <- sample(1:length(xo),round(length(xo)*0.5),replace=FALSE)
x <- rep(NA,length(xo))
x[idata] <- xo[idata] # simulated samples
Output <- impute_data(time,x,harmonic = TRUE,MC = 10,ptrain = 0.8)
}
\references{
Rustomji, P., & Wilkinson, S. N. (2008). Applying bootstrap resampling to quantify uncertainty in fluvial suspended sediment loads estimated using rating curves. Water resources research, 44(9).

van Buuren S, Groothuis-Oudshoorn K (2011). “mice: Multivariate Imputation by Chained Equations in R.” Journal of Statistical Software, 45(3), 1-67. doi:10.18637/jss.v045.i03.

Stephenson AG (2016). Harmonic Analysis of Tides Using TideHarmonics. https://CRAN.R-project.org/package=TideHarmonics.

Moritz S, Bartz-Beielstein T (2017). “imputeTS: Time Series Missing Value Imputation in R.” The R Journal, 9(1), 207–218. doi:10.32614/RJ-2017-009.
}
\author{
Daniel Livsey (2023) ORCID: 0000-0002-2028-6128
}
