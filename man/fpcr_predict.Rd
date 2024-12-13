\name{fpcr_predict}
\alias{fpcr_predict}
\title{Prediction for functional Cox regression models with frailty}
\description{Generates predictions for a functional Cox regression model with frailty, using new functional and scalar predictors, and calculates the concordance index.}
\usage{
fpcr_predict(object, time, event, group = NULL, X, Z, nb, gp)
}
\arguments{
  \item{object}{An \code{fpcr} object, representing a fitted functional Cox regression model.}
  \item{time}{A numeric vector representing the time-to-event or censoring time for new observations.}
  \item{event}{A binary numeric vector indicating whether the event occurred (\code{1}) or was censored (\code{0}) for new observations.}
  \item{group}{(Optional) A grouping variable for frailty effects, such as cluster IDs for new data. If \code{NULL}, no frailty term is included.}
  \item{X}{A matrix of dimensions \( n \times j \), where \( n \) is the number of new observations and \( j \) is the number of grid points for the functional predictor.}
  \item{Z}{A matrix of dimensions \( n \times p \), where \( p \) is the number of scalar predictors for the new observations.}
  \item{nb}{An integer specifying the number of B-spline basis expansion functions used to approximate the functional principal components.}
  \item{gp}{A numeric vector representing the grid points corresponding to the functional predictor.}
}
\value{
A list with the following components:
  \item{predictions}{A numeric vector of predicted values (linear predictors) for the new observations.}
  \item{concordance}{The concordance index for the predictions, measuring the model's predictive accuracy.}
}
\author{
Deniz Inan, Ufuk Beyaztas, Carmen D. Tekwe, Xiwei Chen, and Roger S. Zoh
}
\examples{
## Example Usage
set.seed(123)
# Simulate training data
train_data <- data_generation(n = 100, j = 101, p = 6, gamma0 = rep(0.5, 6),
                            tau = 0.2, mev = 0.5, theta = 2)
# Fit fpcr model
fit <- fpcr(time = train_data$time, event = train_data$event, X = train_data$Xt,
            Z = train_data$Z, nb = 5, gp = train_data$gp)

# Simulate new data for prediction
test_data <- data_generation(n = 100, j = 101, p = 6, gamma0 = rep(0.5, 6),
                            tau = 0.2, mev = 0.5, theta = 2)

# Apply fpcr_predict
predictions <- fpcr_predict(object = fit, time = test_data$time, event = test_data$event,
                             X = test_data$Xt, Z = test_data$Z, nb = 5, gp = test_data$gp)

}
