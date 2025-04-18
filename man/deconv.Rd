% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/denoise.R
\name{deconv}
\alias{deconv}
\title{Deconvolve ecdfs of one sample with others (quadprog version)}
\usage{
deconv(y, X, lambda)
}
\arguments{
\item{y}{response vector}

\item{X}{predictor matrix}

\item{lambda}{penalty parameter}
}
\value{
weights, denoised value and chosen penalty parameter
}
\description{
Deconvolve ecdfs of one sample with others (quadprog version)
}
