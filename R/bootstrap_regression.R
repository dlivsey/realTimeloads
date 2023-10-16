#' Regression parameters estimated using bootstrap resampling
#'
#' Computes uncertainty in regression parameters of y(x) after Rustomji and Wilkinson (2008)
#' @param Calibration data frame with surrogate(s) followed by analyte in last column
#' @param fit_eq equation used to fit y(x), string (e.g, "y ~ x + x2", "y ~ x", "log10(y)~x')
#' @param fit_glm logical to use Generalized Linear Models for models with factor (i.e., categorical) predictors
#' @returns list with bootstrap regression parameters and list output from stats::lm()
#' @section Warning:
#' User should inspect regression residuals and relevant statistics to ensure model form is reasonable, suggested reading: regression diagnostics in Statistical Methods in Water Resources (https://doi.org/10.3133/tm4a3).
#'
#' One can call plot(fit) to view various regression diagnostic plots
#' @section Note:
#' Bias Correction Factor (BCF) is only relevant when analyte is transformed to log units, see https://doi.org/10.3133/tm4a3
#' to convert a model that used log(analyte) back to linear units use: analyte = 10^(f(surrogates)) x BCF
#' @examples
#' \donttest{
#' # linear model
#' x <- 1:10
#' y <- 0.5*x + 10
#' boot <- bootstrap_regression(data.frame(x,y),"y~x")
#' # polynomial model, call to I() needed for squaring x in equation string
#' x <- 1:10
#' y <- x + x^2
#' boot <- bootstrap_regression(data.frame(x,y),"y ~ x+I(x^2)")
#' # power law model
#' # BCF returned since y is transformed to log units
#' x <- 1:10
#' y <- x^0.3
#' boot <- bootstrap_regression(data.frame(x,y),"log10(y)~log10(x)")
#' # multivariate model
#' a <- 1:10
#' b <- a*2
#' c <- a^2*b^3
#' boot <- bootstrap_regression(data.frame(a,b,c),"log10(c)~log10(a)+log10(b)")
#' }
#' @references
#' Rustomji, P., & Wilkinson, S. N. (2008). Applying bootstrap resampling to quantify uncertainty in fluvial suspended sediment loads estimated using rating curves. Water resources research, 44(9).https://doi.org/10.1029/2007WR006088
#'
#' Helsel, D.R., Hirsch, R.M., Ryberg, K.R., Archfield, S.A., and Gilroy, E.J., 2020, #' Statistical methods in water resources: U.S. Geological Survey Techniques and Methods, book 4, chap. A3, 458 p. https://doi.org/10.3133/tm4a3
#' @author Daniel Livsey (2023) ORCID: 0000-0002-2028-6128
#' @export
#'
bootstrap_regression <- function(Calibration,fit_eq,fit_glm=FALSE) {

  # determine if y is transformed into log units
  Log_model <- grepl('log',sub("~.*", "",fit_eq))

  df <- Calibration
  # fit model
  fit <- lm(fit_eq, data=df)

  # For models using non-continous predictors (e.g., categorical or "factor" variables)
  if (fit_glm) {
  fit <- glm(fit_eq, data=df)
  }

  yp <- predict(fit,df)
  #summary(fit)
  res <- as.numeric(fit$residuals)

  n <- length(yp)

  iters <- 2000 # iters = 2000 in Rustomji and Wilkinson (2008)

  dfi <- df # y is updated in loop below
  BCF =  matrix(NA, nrow = iters, ncol = 1) # Bias correction factor for log(y) models, see https://doi.org/10.3133/tm4a3
  # Preallocate matrix with regression parameters
  # Columns: Intercept,m of x,m of x2, etc.
  regP =  data.frame(matrix(NA, nrow = iters, ncol = length(fit$coefficients)))
  names(regP)<-names(fit$coefficients)
  ny <- names(Calibration)[ncol(Calibration)]
  #ny <- strsplit(fit_eq,'~')[[1]][1]

  # Bootstrap regression parameters per Rustomji and Wilkinson (2008)
  if (!Log_model){
      for (i in 1:iters) {
        # get residuals
        resi =  sample(res,n,replace=TRUE)
        # update y for each iteration
        dfi[,ny] = yp+resi
        # fit model
        if (!fit_glm) {
          ifit <- lm(fit_eq,dfi)
        }
        if (fit_glm) {
          ifit <- glm(fit_eq,data=dfi)
        }


        # Store regression parameters
        regP[i,] <- ifit$coefficients
      }
    Regression <- list('regression_data'= Calibration,'fit_eq' = fit_eq,'fit' = fit,"regression_parameters" = regP)
  }


  if (Log_model){
      for (i in 1:iters) {
        # get residuals
        resi =  sample(res,n,replace=TRUE)
        # update y for each iteration
        dfi[,ny] = 10^(yp+resi)
        # fit model
        if (!fit_glm) {
          ifit <- lm(fit_eq,dfi)
        }
        if (fit_glm) {
          ifit <- glm(fit_eq,data=dfi)
        }

        # Store regression parameters
        regP[i,] <- ifit$coefficients
        # compute Bias correction factor for log10(y) models only
        BCF[i] = mean(10^ifit$residuals)
      }
      Regression <- list('regression_data'= Calibration, 'fit_eq'= fit_eq,'fit' = fit,"regression_parameters" = regP,'BCF'= BCF)
  }

  return(Regression)
}
