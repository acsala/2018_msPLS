#' @title Outer Weights
#'
#' @details
#' Internal function. \code{get_weights} is called by \code{plspm}
#'
#' @note
#' Calculate outer weights (under Lohmoller's algorithm)
#'
#' @param X scaled data
#' @param path_matrix matrix with path connections
#' @param blocks list with variables in each block
#' @param specs list with algorithm specifications
#' @return list of outer weights, ODM, iter
#' @export
#' @template internals
#' @keywords internal

get_weights_p <- function(X, path_matrix, blocks, specs, penalization,
                          nonzero, lambda, warning_non_convergence = TRUE, alpha)
{

  # X <- Data
  # specs <- s_satpls$model$specs
  # specs$modes <- c("B","B")
  # penalization <- "enet"
  # nonzero <- 40
  # lambda <- 1
  # warning_non_convergence = TRUE

  lambda <- rep(lambda,length(specs$modes))
  nonzero <- rep(nonzero,length(specs$modes))

  lvs = nrow(path_matrix)
  mvs = ncol(X)
  sdv = sqrt((nrow(X)-1) / nrow(X))   # std.dev factor correction
  blockinds = indexify(blocks)

  # outer design matrix 'ODM' and matrix of outer weights 'W'
  ODM = list_to_dummy(blocks)
  W = ODM %*% diag(1/(apply(X %*% ODM, 2, sd)*sdv), lvs, lvs)
  w_old = rowSums(W)
  iter = 1

  repeat
  {
    # external estimation of LVs 'Y'
    Y = X %*% W
    Y = scale(Y) * sdv
    # matrix of inner weights 'e'
    E <- switch(specs$scheme,
                "centroid" = sign(cor(Y) * (path_matrix + t(path_matrix))),
                "factorial" = cor(Y) * (path_matrix + t(path_matrix)),
                "path" = get_path_scheme(path_matrix, Y))
    # internal estimation of LVs 'Z'
    Z = Y %*% E
#    Z = Z %*% diag(1/(apply(Z,2,sd)*sdv), lvs, lvs)  # scaling Z
    # computing outer weights 'w'
    for (j in 1:lvs)
    {
      if (specs$modes[j] == "A")
        W[blockinds==j,j] <- (1/nrow(X)) * Z[,j] %*% X[,blockinds==j]
      if (specs$modes[j] == "B")
        #W[blockinds==j,j] <- solve.qr(qr(X[,blockinds==j]), Z[,j])

        #inject penalization models
        W[blockinds==j,j] <-
          switch(penalization,
                 "none" = t(solve(t(X[,blockinds==j]) %*% X[,blockinds==j])%*%
                              (t(X[,blockinds==j]) %*% Z[,j])),
                 "ust" = get_ust(X[,blockinds==j],Z[,j],nonzero[j]),
                 "enet" = calculateVectorEnet(X[,blockinds==j], Z[,j],
                                              lambda[j], nonzero[j]),
                 "enet2" = get_enet2(X[,blockinds==j], Z[,j],
                                     lambda[j], nonzero[j], alpha)
                 )

    }

    w_new = rowSums(W)
    w_dif = sum((abs(w_old) - abs(w_new))^2)
    if (w_dif < specs$tol || iter == specs$maxiter) break
    w_old = w_new
    iter = iter + 1
  } # end repeat

  # preparing results
  if (iter == specs$maxiter && warning_non_convergence) {
    print(paste("Iterative process did not converge with 'maxiter'=",
                specs$maxiter, ", 'tol'=", specs$tol, ", 'nonzero'=",
                nonzero, ", and 'lambda'=", lambda,
                ". Results might be suboptimal.",
                sep=""))
  }

  # preparing results
  W = W %*% diag(1/(apply(X %*% W, 2, sd)*sdv), lvs, lvs)
  w_new = rowSums(W)
  names(w_new) = colnames(X)
  dimnames(W) = list(colnames(X), rownames(path_matrix))
  results = list(w = w_new, W = W, ODM = ODM, iter = iter)

  # # preparing results
  # if (iter == specs$maxiter) {
  #   results = NULL
  # } else {
  #   W = W %*% diag(1/(apply(X %*% W, 2, sd)*sdv), lvs, lvs)
  #   w_new = rowSums(W)
  #   names(w_new) = colnames(X)
  #   dimnames(W) = list(colnames(X), rownames(path_matrix))
  #   results = list(w = w_new, W = W, ODM = ODM, iter = iter)
  # }
  # output
  results
}
