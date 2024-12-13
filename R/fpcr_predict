fpcr_predict <- function(object, time, event, group = NULL, X, Z, nb, gp)
{

  n <- dim(X)[1]
  j <- dim(X)[2]

  fpca <- object$fpca
  fscore <- getPCA_test(fpca, X)
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

  cox_model <- object$model
  preds <- predict_coxme(cox_model, newdata=model_matrix)
  conc <- 1 - concordance(Surv(time, event) ~ preds, model_matrix)$concordance

  return(list(predictions = preds, concordance = conc))
}
