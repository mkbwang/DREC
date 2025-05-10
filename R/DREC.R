


#' Denoise metagenomics count table
#'
#' @param X raw count matrix, rows are samples and columns are features
#' @param lambda penalty parameter
#' @param cv.fold number of folds for cross validation
#' @param type l1 (absolute loss) or l2 (quadratic loss)
#' @param ncores number of cores for parallel computing
#' @returns denoised count matrix
#' @export
drec <- function(X, lambda=0.1, cv.fold=10, type=c("l1", "l2"), ncores=1){

    type <- match.arg(type)
    X_cdf <- count2ecdf(X)
    cat("denoising ECDFs...\n")
    denoised_cdf <- denoise(X_cdf, lambda=lambda, cv.fold=cv.fold, type=type,
                            ncores=ncores)
    cat("mapping ECDFs back to counts ... \n")
    denoised_count <- ecdf2count(Q_denoised=denoised_cdf,
                                 Q_observed=X_cdf,
                                 count_observed=X,
                                 ncores=ncores)

    return(denoised_count)

}
