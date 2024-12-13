fpcr <- function(time, event, group = NULL, X, Z, weights = NULL, nb, gp)
{

  n <- dim(X)[1]
  j <- dim(X)[2]

  fpca <- getPCA(data = X, nbasis = nb, gp = gp)
  fscore <- fpca$PCAscore
  evalbase <- fpca$evalbase
  PCAcoef <- fpca$PCAcoef

  model_mat <- cbind(fscore, Z)
  if(is.null(group)){
    id <- 1:n
  }else{
    id <- group
  }
  model_matrix <- data.frame(cbind(time, event, id=id, model_mat))
  for(i in 4:dim(model_matrix)[2])
    colnames(model_matrix)[i] = paste("V", (i-3), sep = "")

  colnames(model_matrix)[1] <- "time"

  var_name <- paste(c(colnames(model_matrix)[-(1:3)], "(1|id)"), collapse = "+")

  cox_formula <- as.formula(paste("Surv(time, event)~", var_name, sep = ""))
  if(is.null(weights)){
    cox_model <- coxme(cox_formula, data = model_matrix)
  }else{
    cox_model <- coxme(cox_formula, data = model_matrix, weights = weights,
                       control = coxme.control(
                         eps = 1e-05, toler.chol = .Machine$double.eps^0.75,
                         iter.max = 5000, sparse.calc = NULL,
                         optpar = list(method = "CG",
                                       control = list(reltol = 1e-3, maxit = 5000))))
  }

  cox_result <- summary(cox_model)
  cox_coefs <- cox_result$coefficients
  coef_f <- cox_coefs[,1]

  conc <- 1 - concordance(Surv(time, event) ~
                            cox_model$linear.predictor, model_matrix)$concordance

  bhat_t <- evalbase %*% (PCAcoef$coefs %*% as.matrix(coef_f[1:dim(fscore)[2]]))
  gamma_hat <- coef_f[-(1:dim(fscore)[2])]

  return(list(bhat = bhat_t, gammahat = gamma_hat, concordance = conc,
              model = cox_model, fpca = fpca))

  }
