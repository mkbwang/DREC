

#' Deconvolve ecdfs of one sample with others
#'
#' @param y response vector
#' @param X predictor matrix
#' @param lambdas penalty parameters to choose from
#' @importFrom quadprog solve.QP
#' @returns weights, denoised value and chosen penalty parameter
deconv <- function(y, X, lambdas=c(1e-3, 1e-2, 1e-1)){

  n_predictor <- ncol(X)
  n_sample <- nrow(X)

  train_sample_id <- sample(n_sample, n_sample * 0.8)
  X_train <- X[train_sample_id, ]
  X_validation <- X[-train_sample_id, ]

  y_train <- y[train_sample_id]
  y_validation <- y[-train_sample_id]

  # quadratic optimization parameters
  dvec <- t(X_train) %*% y_train

  # linear constraint (one equality and others are inequality)
  Amat <-cbind(rep(1, n_predictor), diag(n_predictor))
  bvec <- c(1, rep(0, n_predictor))

  mses_validation <- rep(0, length(lambdas))
  weights <- matrix(0, ncol=length(lambdas), nrow=n_predictor)

  for (j in 1:length(lambdas)){

    penalty <- lambdas[j] * (0.25*n_sample) * diag(nrow=n_predictor)
    Dmat <- t(X_train) %*% X_train + penalty
    result <- solve.QP(Dmat, dvec, Amat, bvec, meq=1)
    weight <- result$solution
    weight[weight < 1e-8] <- 0

    weights[, j] <- weight
    y_validation_pred <- X_validation %*% weight
    mse_validation <- mean((y_validation - y_validation_pred)^2)
    mses_validation[j] <- mse_validation

  }

  best_choice <- which.min(mses_validation)
  best_lambda <- lambdas[best_choice]
  best_weight <- weights[, best_choice]

  y_denoise <- X %*% best_weight

  return(list(y_denoise=as.vector(y_denoise), lambda=best_lambda, weight=best_weight))

}




#' deconvolve every row of a matrix with ECDF values
#'
#' @param Q input matrix whose entries are ECDF values
#' @param lambdas penalty parameters to choose from
#'
#' @returns denoised matrix
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#'
#' @export
denoise <- function(Q, lambdas=c(1e-3, 1e-2, 1e-1)){

  nsamples <- nrow(Q)
  numCores <- detectCores() - 1
  cl <- makeCluster(numCores)
  registerDoParallel(cl)

  i <- 1
  denoised_Q <- foreach(i=1:nsamples, .combine=rbind,
                        .export="deconv",
                        .packages="quadprog") %dopar%{
    q_i <- Q[i, ]
    Q_it <- Q[-i, ] |> t()
    denoise_result <- deconv(y=q_i, X=Q_it, lambdas=lambdas)
    denoise_result$y_denoise
  }

  stopCluster(cl)
  denoised_Q
}

