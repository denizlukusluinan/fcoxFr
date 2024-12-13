\name{data_generation}
\alias{data_generation}
\title{Generate data for functional cox regression model with frailty}
\description{This function is used to simulate data for the functional cox
regression model with frailty \deqn{
h_i(t \vert \bm{Z}_i, X_i, \bm{w}) = h_0(t) \exp \left(\bm{Z}_i^\top \bm{\gamma}
+ \int_{\mathcal{I}} X_i(s) \beta(s) ds + \bm{U}_i^\top \bm{w} \right),} where \eqn{h_0(t)} is the baseline hazard function, \eqn{h_i(t \vert \bm{Z}_i, X_i, \bm{w})} is the hazard function incorporating frailty for subject \eqn{i}, \eqn{\bm{\gamma} = (\gamma_1, \ldots, \gamma_p )^\top}is the vector of regression coefficients, \eqn{\beta(s)} is the functional regression coefficient, and the frailty for subject \eqn{i}, \eqn{\bm{w}_i}, captures the effect of unobserved factors influencing the risk of the event for subject \eqn{i}. Here, \eqn{\bm{U}_i = (U_{1i}, \ldots, U_{ki})^\top} is a \eqn{k}-dimensional indicator vector that indicates the family of \eqn{i}-th individual.
}
\usage{
data_generation(n, j, p, gamma0, tau, mev, theta)
}
\arguments{
  \item{n}{An integer, specifying the number of observations to be generated.}
  \item{j}{An integer, denoting the number of grid points, i.e., a fine grid on the interval [0, 1].}
  \item{p}{An integer, denoting the number of scalar predictors to be generated.}
  \item{gamma0}{A vector, containing the parameters for the scalar predictors.}
  \item{tau}{A numeric value, denoting the rate parameter of exponential censoring distribution.}
  \item{mev}{A numeric value, denoting the measurement error variance.}
  \item{theta}{A numeric value, denoting the scale parameter of the gamma distribution.}
}
\value{
  \item{time}{A numeric vector of length \eqn{n} representing the observed time-to-event or censoring times for each individual. This is the minimum of the failure time and the censoring time.}
  \item{event}{A binary numeric vector of length nn indicating whether an event occurred (11) or if the observation was censored (00).}
  \item{Xt}{A matrix representing the functional predictor observations for each individual, including measurement error. Each row corresponds to the observed functional predictor values at \eqn{j} grid points.}
  \item{Xt_true}{A matrix representing the true (unobserved) functional predictor values for each individual, without measurement error.}
  \item{Z}{A matrix representing the scalar predictor values for each individual.}
  \item{beta}{A numeric vector of length \eqn{j} representing the true functional coefficient values \eqn{\beta(s)} evaluated at the \eqn{j} grid points.}
  \item{gamma}{A numeric vector containing the regression coefficients for the scalar predictors.}
  \item{gp}{A numeric vector of length \eqn{j} containing the grid points where the functional predictor and functional coefficient are evaluated. These values range from 0 to 1, equally spaced.}
}
\references{
Kong, D., Ibrahim, J. G., Lee, E. and Zhu, H. (2018), b
