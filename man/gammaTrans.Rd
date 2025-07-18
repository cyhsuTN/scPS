% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICC_gammaTrans.R
\name{gammaTrans}
\alias{gammaTrans}
\title{Parameter Transform in Gamma Distribution}
\usage{
gammaTrans(mean, q95)
}
\arguments{
\item{mean}{The mean of a gamma Distribution.}

\item{q95}{0.95 quantile of a gamma Distribution.}
}
\description{
A function to calculate shape and scale parameters in gamma distribution,
given its mean and 0.95 quantile.
}
\examples{
# gammaTrans(mean=0.01, q95=0.1)
}
