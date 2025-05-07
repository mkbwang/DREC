

#' Deconvolve ecdfs of one sample with others (quadprog version)
#'
#' @param y response vector
#' @param X predictor matrix
#' @param lambda penalty parameter
#' @importFrom quadprog solve.QP
#' @returns weights, denoised value and chosen penalty parameter
deconv <- function(y, X, lambda){

  n_predictor <- ncol(X)
  n_sample <- nrow(X)

  # quadratic optimization parameters
  dvec <- t(X) %*% y

  # linear constraint (one equality and others are inequality)
  Amat <-cbind(rep(1, n_predictor), diag(n_predictor))
  bvec <- c(1, rep(0, n_predictor))


  penalty <- lambda * (0.25*n_sample) * diag(nrow=n_predictor)
  Dmat <- t(X) %*% X + penalty
  result <- solve.QP(Dmat, dvec, Amat, bvec, meq=1)
  weight <- result$solution
  weight[weight < 1e-8] <- 0

  y_denoise <- X %*% weight

  return(list(y_denoise=as.vector(y_denoise), weight=weight))

}

#' Deconvolve ecdfs of one sample with others (l1 loss)
#'
#' @param y response vector
#' @param X predictor matrix
#' @param lambda penalty parameter
#' @importFrom Matrix Matrix
#' @importFrom osqp osqp osqpSettings
#' @returns weights, denoised value and chosen penalty parameter
deconv.l1 <- function(y, X, lambda){

  n_predictor <- ncol(X)
  n_sample <- nrow(X)
  Pmat <- Matrix(diag(c(rep(1, n_predictor), rep(0, n_sample))) * lambda *2,
                 sparse=TRUE)
  qvec <- c(rep(0, n_predictor), rep(1, n_sample))

  # sum to one constraint
  A_constraint0 <- matrix(c(rep(1, n_predictor), rep(0, n_sample)), nrow=1)
  l_constraint0 <- 1
  u_constraint0 <- 1

  # absolute value constraint
  A_constraint1 <- cbind(X, diag(n_sample))
  l_constraint1 <- y
  u_constraint1 <- rep(Inf, n_sample)

  A_constraint2 <- cbind(-X, diag(n_sample))
  l_constraint2 <- -y
  u_constraint2 <- rep(Inf, n_sample)


  # nonegative constraint
  A_constraint3 <- cbind(diag(n_predictor), matrix(0, nrow=n_predictor, ncol=n_sample))
  l_constraint3 <- rep(0, n_predictor)
  u_constraint3 <- rep(Inf, n_predictor)

  Amat <- Matrix(rbind(A_constraint0, A_constraint1, A_constraint2, A_constraint3),
                 sparse=TRUE)
  l_vec <- c(l_constraint0, l_constraint1, l_constraint2, l_constraint3)
  u_vec <- c(u_constraint0, u_constraint1, u_constraint2, u_constraint3)

  settings <- osqpSettings(alpha = 1.0, max_iter=8000)
  model <- osqp(Pmat, qvec, Amat, l_vec, u_vec, settings)


  result <- model$Solve()
  weights <- result$x[1:n_predictor]
  weights[weights < 1/n_predictor/10] = 0
  weights <- weights / sum(weights)

  return(weights)

}



#' Cross validation of deconvolution
#'
#' @param y response vector
#' @param X predictor matrix
#' @param lambdas penalty parameters to try
#' @returns weights, mean error for the and chosen penalty parameter
#' @importFrom stats median
cv.deconv <- function(y, X, lambdas = c(1e-3, 1e-2, 1e-1)){

  n_sample <- nrow(X)
  n_feature <- ncol(X)
  cv_list <- split(sample(n_sample),
                   cut(seq(1, n_sample), breaks=10))

  validation_error <- rep(0, length(lambdas))
  weights_list <- list()

  for(j in 1:length(lambdas)){

    lambda <- lambdas[j]
    errors <- rep(0, 10)
    weights_mat <- matrix(0, nrow=10, ncol=n_feature)

    for (k in 1:10){
      y_validation <- y[cv_list[[k]]]
      X_validation <- X[cv_list[[k]], ]

      y_train <- y[-cv_list[[k]]]
      X_train <- X[-cv_list[[k]], ]

      fitted_result <- deconv(y_train, X_train, lambda=lambda)
      fitted_weights <- fitted_result$weight
      weights_mat[k, ] <- fitted_weights

      y_validation_pred <- X_validation %*% fitted_weights
      errors[k] <- median(abs(y_validation - y_validation_pred))
    }

    validation_error[j] <- mean(errors)
    weights_list[[j]] <- weights_mat

  }

  return(list(weights_list=weights_list,
              validation_error=validation_error,
              lambdas=lambdas))

}


#' deconvolve every row of a matrix with ECDF values
#'
#' @param Q input matrix whose entries are ECDF values
#' @param lambdas penalty parameters to choose from
#' @param ncores number of cores for parallel processing, by default 1
#'
#' @returns denoised matrix
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom parallelly availableCores
#' @importFrom foreach foreach %dopar%
#'
#' @export
denoise <- function(Q, lambdas=c(1e-3, 1e-2, 1e-1, 1, 10, 100), ncores=1){

  nsamples <- nrow(Q)
  numCores <- min(availableCores(), ncores)
  cl <- makeCluster(numCores)
  registerDoParallel(cl)

  i <- 1
  denoised_Q <- foreach(i=1:nsamples, .combine=rbind,
                        .export=c("deconv", "cv.deconv"),
                        .packages="quadprog") %dopar%{
    q_i <- Q[i, ]
    Q_it <- Q[-i, ] |> t()
    cv_result <- cv.deconv(y=q_i,
                           X=Q_it,
                           lambdas=lambdas)
    best_lambda <- lambdas[which.min(cv_result$validation_error)]
    cv_weights <- cv_result$weights_list[[which.min(cv_result$validation_error)]]
    selection_prob <- colMeans(cv_weights > 0)
    subset_sample <- which(selection_prob == 1)
    subset_Q_it <- Q_it[, subset_sample]

    denoise_result <- deconv(y=q_i, X=subset_Q_it, lambda=best_lambda)
    denoise_result$y_denoise
  }

  stopCluster(cl)
  denoised_Q
}

