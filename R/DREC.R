


#' Denoise metagenomics count table
#'
#' @param X raw count matrix, rows are samples and columns are features
#' @param lambdas penalty parameters to choose from
#' @param ncores number of cores for parallel computing
#'
#' @returns denoised count matrix
#' @export
drec <- function(X, lambdas=c(1e-3, 1e-2, 1e-1), ncores=1){

    X_cdf <- count2ecdf(X)
    cat("denoising ECDFs...\n")
    denoised_cdf <- denoise(X_cdf, lambdas, ncores=ncores)
    cat("mapping ECDFs back to counts ... \n")
    denoised_count <- ecdf2count(Q_denoised=denoised_cdf,
                                 Q_observed=X_cdf,
                                 count_observed=X,
                                 ncores=ncores)

    return(denoised_count)

}
