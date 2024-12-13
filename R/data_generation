data_generation <- function(n, j, p, gamma0, tau, mev, theta)
{

  grid_points <- seq(0, 1, length.out = j)
  base_mat <- matrix(, nrow = 22, ncol = j)
  for(i in 1:10)
    base_mat[i,] <- sqrt(2) * sin(pi * (2*i-1) * grid_points)
  for(i in 11:20)
    base_mat[i,] <- sqrt(2) * cos(pi * (2*i-21) * grid_points)
  base_mat[21,] <- 1
  base_mat[22,] <- grid_points

  Sigma <- matrix(0, 22, 22)
  Sigma[1,1] <- 1
  Sigma[11,11] <- 1
  for(i in 2:10)
    Sigma[i,i] <- (i-1)^(-2)
  for(i in 12:20)
    Sigma[i,i] <- (i-11)^(-2)
  Sigma[21,21] <- 1
  Sigma[22,22] <- 1

  beta_fun <- function(t)
    {
    sigma <- 0.3
    out <- 0.3 * (sin(pi*t) - cos(pi*t) + sin(3*pi*t) - cos(3*pi*t) + sin(5*pi*t)/9 -
                    cos(5*pi*t)/9 + sin(7*pi*t)/16 - cos(7*pi*t)/16 + sin(9*pi*t)/25 -
                    cos(9*pi*t)/25 + (1/sqrt(2*pi)/sigma) * exp(-(t-0.5)^2/2/sigma^2))
    return(out)
  }

  beta_t <- beta_fun(grid_points)

  corr_mat <- matrix(0, nrow(Sigma), p)
  corr_mat[1, 1:p] <- 0.1
  corr_mat[1:p, 1] <- 0.1

  SigmaZ <- 0.5^t(sapply(1:p, function(i, k) abs(i-k), 1:p))
  Bigsigma <-rbind(cbind(Sigma, corr_mat), cbind(t(corr_mat), SigmaZ))
  mu <- rep(0, nrow(Bigsigma))
  data <- mvrnorm(n, mu, Bigsigma)
  us <- data[, 1:nrow(Sigma)]
  Xt_smooth <- us %*% base_mat
  options(warn = -1)
  msm_err <- matrix(rnorm(n*j, 0, sqrt(mev)), n, j)
  Xt_me <- Xt_smooth + msm_err
  Z <- data[, (nrow(Sigma) + 1):(nrow(Sigma) + p)]
  params <- exp(Xt_me %*% beta_t * (grid_points[2]-grid_points[1]) + Z %*% gamma0)
  frailty <- rgamma(n, shape = 1/theta, scale = theta)
  failure_time <- rexp(n, rate = frailty* params)
  cens_time <- rexp(n, rate=tau)
  event <- rep(0, n)
  event <- as.numeric(failure_time < cens_time)
  time <- failure_time * event + cens_time * (rep(1, n) - event)
  censrate<-1-sum(event)/n
  return(list(time=time, event=event, Xt = Xt_me, Xt_true=Xt_smooth, Z=Z,
              beta=beta_t, gamma=gamma0, gp=grid_points))
}
