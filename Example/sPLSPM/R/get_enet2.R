#' @title Calculate elastic net
#'
#' @details
#' Calculate elastic net with glmnet package
#' 
#' @import glmnet
#' @template internals
#' @keywords internal

# X.m = X[,blockinds==j] 
# z = Z[,j]

get_enet2 = function(X.m, z, lambda, nonzero, alpha = 0.5){
  
  s_time <- system.time(
                fit <- glmnet(X.m, z, family = "gaussian", alpha = alpha) 
              )

  tLL <- fit$nulldev - deviance(fit)
  k <- fit$df
  n <- fit$nobs

  BIC<-log(n)*k - tLL

  b_hat <- fit$beta[,which(BIC==min(BIC))]

  b_hat
}