

#' Deconvolve ecdfs of one sample with others (l2 loss)
#'
#' @param y response vector
#' @param X predictor matrix
#' @param lambda penalty parameter
#' @importFrom Matrix Matrix
#' @importFrom osqp osqp osqpSettings
#' @returns weights, denoised value and chosen penalty parameter
deconv.l2 <- function(y, X, lambda){

  n_predictor <- ncol(X)
  n_sample <- nrow(X)

  # quadratic optimization parameters
  qvec <- -t(X) %*% y
  penalty <- lambda * (0.25*n_sample) * diag(nrow=n_predictor)
  Pmat <- Matrix(t(X) %*% X + penalty, sparse=TRUE)

  # sum to one constraint
  A_constraint0 <- matrix(rep(1, n_predictor), nrow=1)
  l_constraint0 <- 1
  u_constraint0 <- 1

  # nonegative constraint
  A_constraint1 <- diag(n_predictor)
  l_constraint1 <- rep(0, n_predictor)
  u_constraint1 <- rep(Inf, n_predictor)

  Amat <- Matrix(rbind(A_constraint0, A_constraint1),
                 sparse=TRUE)
  l_vec <- c(l_constraint0, l_constraint1)
  u_vec <- c(u_constraint0, u_constraint1)

  settings <- osqpSettings(alpha = 1.0, max_iter=8000, verbose=FALSE)
  model <- osqp(Pmat, qvec, Amat, l_vec, u_vec, settings)

  result <- model$Solve()
  weights <- result$x[1:n_predictor]
  weights[weights < 1/n_predictor/10] = 0
  weights <- weights / sum(weights)

  return(weights)

}




#' Deconvolve ecdfs of one sample with others (l1 loss)
#'
#' @param y response vector
#' @param X predictor matrix
#' @param lambda penalty parameter
#' @importFrom Matrix Matrix
#' @importFrom osqp osqp osqpSettings
#' @returns weights, denoised value and chosen penalty parameter
#' @export
deconv.l1 <- function(y, X, lambda){

  n_predictor <- ncol(X)
  n_sample <- nrow(X)
  Pmat <- Matrix(diag(c(rep(1, n_predictor), rep(0, n_sample))) * (0.25*n_sample) * lambda,
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

  settings <- osqpSettings(alpha = 1.0, max_iter=8000, verbose=FALSE)
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
#' @param lambda penalty parameter
#' @param cv.fold number of folds for repeated fitting
#' @param type choose between absolute loss or square loss when fitting deconvolution
#' @returns weights, mean error for the and chosen penalty parameter
#' @importFrom stats median
#' @export
cv.deconv <- function(y, X, lambda = 0.1,
                      cv.fold=5, type=c("l1", "l2")){

  regression_type <- match.arg(type)
  print(regression_type)
  if (regression_type == "l1"){
    deconv_func <- deconv.l1
  } else{
    deconv_func <- deconv.l2
  }

  n_sample <- nrow(X)
  n_feature <- ncol(X)
  cv_list <- split(sample(n_sample),
                   cut(seq(1, n_sample), breaks=cv.fold))

  weights_mat <- matrix(0, nrow=cv.fold, ncol=n_feature)

  for (k in 1:cv.fold){
    y_validation <- y[cv_list[[k]]]
    X_validation <- X[cv_list[[k]], ]

    y_train <- y[-cv_list[[k]]]
    X_train <- X[-cv_list[[k]], ]

    fitted_weights <- deconv_func(y_train, X_train, lambda=lambda)
    weights_mat[k, ] <- fitted_weights

  }

  return(weights_mat)

}



#' deconvolve every row of a matrix with ECDF values
#'
#' @param Q input matrix whose entries are ECDF values
#' @param lambda penalty parameter
#' @param cv.fold number of folds for cross validation
#' @param type choose between absolute loss or square loss when fitting deconvolution
#' @param ncores number of cores for parallel processing, by default 1
#' @returns denoised matrix
#' @importFrom parallel detectCores makeCluster stopCluster clusterExport
#' @importFrom doParallel registerDoParallel
#' @importFrom parallelly availableCores
#' @importFrom foreach foreach %dopar%
#'
#' @export
denoise <- function(Q, lambda=0.1, cv.fold=10, type=c("l1", "l2"), ncores=1){

  type <- match.arg(type)
  nsamples <- nrow(Q)
  numCores <- min(availableCores(), ncores)
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  clusterExport(cl, varlist=c("deconv.l1", "deconv.l2", "cv.deconv"))
  i <- 1
  denoised_Q <- foreach(i=1:nsamples, .combine=rbind,
                        .export=c("deconv.l1", "deconv.l2", "cv.deconv"),
                        .packages=c("osqp", "Matrix")) %dopar% {
    q_i <- Q[i, ]
    Q_it <- Q[-i, ] |> t()

    cv_weights <- cv.deconv(y=q_i, X=Q_it, lambda=lambda, cv.fold=cv.fold, type=type)

    cv_mean_weights <- colMeans(cv_weights)

    y_denoise <- as.vector(Q_it %*% cv_mean_weights)
    y_denoise

  }

  stopCluster(cl)
  denoised_Q

}

