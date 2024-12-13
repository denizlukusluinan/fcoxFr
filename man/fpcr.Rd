\name{fpcr}
\alias{fpcr}
\title{Functional principal component Cox regression with frailty}
\description{Estimates the coefficients for a functional Cox regression model with optional frailty, incorporating scalar covariates and functional principal components for the functional predictor.}
\usage{
fpcr(time, event, group = NULL, X, Z, weights = NULL, nb, gp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{A numeric vector representing the time-to-event or censoring time for each observation.}
  \item{event}{A binary numeric vector indicating whether the event occurred (\code{1}) or was censored (\code{0}).}
  \item{group}{(Optional) A grouping variable for frailty effects, such as cluster IDs. If \code{NULL}, no frailty term is included.}
  \item{X}{A matrix of dimensions \( n \times j \), where \( n \) is the number of observations, and \( j \) is the number of grid points for the functional predictor.}
  \item{Z}{A matrix of dimensions \( n \times p \), where \( p \) is the number of scalar predictors.}
  \item{weights}{(Optional) A vector of weights for each observation. If \code{NULL}, equal weights are used.}
  \item{nb}{An integer specifying the number of B-spline basis expansion functions used to approximate the functional principal components.}
  \item{gp}{A numeric vector representing the grid points corresponding to the functional predictor.}
}
\value{
A list with the following components:
  \item{bhat}{A numeric vector representing the estimated functional coefficient \(\hat{\beta}(s)\) evaluated at the grid points \code{gp}.}
  \item{gammahat}{A numeric vector of estimated scalar regression coefficients for the covariates \code{Z}.}
  \item{concordance}{The concordance index for the fitted Cox model, measuring the model's predictive accuracy.}
  \item{model}{The fitted \code{coxme} model object.}
  \item{fpca}{The result of the FPCA performed on the functional predictor, including principal component scores, eigenvalues, and eigenfunctions.}
}
\author{
Deniz Inan, Ufuk Beyaztas, Carmen D. Tekwe, Xiwei Chen, and Roger S. Zoh
}
\examples{
## Example Usage
set.seed(123)
# Simulate data
sim_data <- data_generation(n = 100, j = 101, p = 6, gamma0 = rep(0.5, 6),
                            tau = 0.2, mev = 0.5, theta = 2)
# Apply fpcr
fit <- fpcr(time = sim_data$time, event = sim_data$event, X = sim_data$Xt,
            Z = sim_data$Z, nb = 10, gp = sim_data$gp)
}
