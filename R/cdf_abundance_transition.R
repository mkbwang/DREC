

#' Calculate ECDF given a vector of counts
#'
#' @param x a vector of values
#' @importFrom stats ecdf
#' @returns the empirical CDF representing the rankings in x
calcecdf <- function(x){

    ecdf_func <- ecdf(x)
    ecdf_vals <- ecdf_func(x)

    return(ecdf_vals)
}



#' Calculate ECDF for every column of a count matrix
#'
#' @param X a matrix of values. rows are samples and columns are features
#' @export
count2ecdf <- function(X){
    ecdf_mat <- apply(X, 2, calcecdf)
    return(ecdf_mat)
}





#' Map denoised quantiles back to counts
#' @param Q_denoised quantiles after denoising
#' @param Q_observed the quantiles of the original observed values
#' @param values the original observed values
#' @returns the interpolated denoised values
calcval <- function(Q_denoised, Q_observed, values){

    nsample <- length(Q_denoised)

    unique_Q <- unique(Q_observed) |> sort()
    unique_values <- unique(values) |> sort()

    rep_unique_Q <- t(replicate(n=nsample, unique_Q))
    lindex <- rowSums(rep_unique_Q <= Q_denoised)
    uindex <- lindex + 1

    Q_lower <- unique_Q[lindex]
    Q_upper <- unique_Q[uindex]

    lower_bound <- unique_values[lindex]
    upper_bound <- unique_values[uindex]

    # linear interpolation
    denoised_values <- lower_bound  + (Q_denoised - Q_lower) / (Q_upper - Q_lower) *
        (upper_bound - lower_bound)

    return(denoised_values)

}


#' Map denoised quantiles back to counts for each column
#'
#' @param Q_denoised denoised empirical CDF table
#' @param Q_observed original empirical CDF table
#' @param count_observed observed counts table
#' @param ncores number of cores for parallel processing, by default 1
#'
#' @returns denoised counts
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom parallelly availableCores
#' @export
ecdf2count <- function(Q_denoised, Q_observed, count_observed,
                       ncores=1){

    nfeatures <- ncol(Q_denoised)
    numCores <- min(availableCores(), ncores)
    cl <- makeCluster(numCores)
    registerDoParallel(cl)

    i <- 1
    denoised_counts <- foreach(i=1:nfeatures, .combine=cbind,
                          .export="calcval") %dopar%{
                              denoised_count <- calcval(Q_denoised[, i],
                                                        Q_observed[, i],
                                                        count_observed[, i])
                              denoised_count
                          }

    stopCluster(cl)
    return(denoised_counts)

}


