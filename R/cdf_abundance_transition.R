

#' @param x a vector of values
#' @importFrom stats ecdf
#' @returns the empirical CDF representing the rankings in x
calcecdf <- function(x){

    ecdf_func <- ecdf(x)
    ecdf_vals <- ecdf_func(x)

    return(ecdf_vals)
}


#' @param X a matrix of values. rows are samples and columns are features
#' @details
#' Calculate empirical CDF values for entries in each column
#'
count2ecdf <- function(X){
    ecdf_mat <- apply(X, 2, calcecdf)
}





#' @param Q_denoised quantiles after denoising
#' @param Q_observed the quantiles of the original observed values
#' @param values the original observed values
#' @returns the interpolated denoised values
calcval <- function(Q_denoised, Q_observed, values){

}




