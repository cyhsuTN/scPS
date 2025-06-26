#' Parameter Transform in Gamma Distribution
#' @description A function to calculate shape and scale parameters in gamma distribution,
#' given its mean and 0.95 quantile.
#' @importFrom stats qgamma
#' @importFrom nleqslv nleqslv

#' @param mean The mean of a gamma Distribution.
#' @param q95 0.95 quantile of a gamma Distribution.

#' @examples # gammaTrans(mean=0.01, q95=0.1)

#' @export
gammaTrans <- function(mean, q95) {
  eq1 <- function(pb) {
    b <- exp(pb)
    a <- mean/b
    qgamma(0.95, shape=a, scale=b) - q95
  }
  root.solve <- nleqslv(x=log(0.05), fn=eq1, control=list(ftol=1e-6, maxit=15))
  b <- exp(root.solve$x)
  a <- mean/b
  round(c(shape=a, scale=b),3)
}





