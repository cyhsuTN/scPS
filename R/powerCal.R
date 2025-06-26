#' Power Calculation for DEA of a Single Cell Type between Two Independent Groups
#' @description A function to calculate powers for DEA of a single cell type between two independent groups.
#' @importFrom nleqslv nleqslv
#' @importFrom stats pt qt

#' @param ns Number of cells of interest per subject in control and experimental groups.
#' @param ms Numbers of subjects in the control and experimental groups, respectively.
#' @param vvmean1 Gene means of candidate genes in the control group.
#' @param FC Fold changes of candidate genes between the control and the experimental groups.
#' @param vvrho Cell-cell correlations of candidate genes within subjects.
#' @param hf A curve showing the relationship between gene means and gene standard deviations.
#' @param FDR False discovery rate.

#' @examples set.seed(12345)
#' @examples vvmean1 <- rep(1, 1000)
#' @examples FC <- c(rep(2, 50), rep(1, 950))
#' @examples ab <- gammaTrans(mean=0.01, q95=0.1) # Output shape and scale parameters
#' @examples vvrho <- rgamma(1000, shape=ab[1], scale=ab[2])
#' @examples hf <- function(x) sqrt(x*(1+3*x))
#' @examples powerCal(ns=c(30,50), ms=c(1,1)*10, vvmean1, FC, vvrho, hf, FDR=0.05)

#' @return \item{power}{Statistical power}
#' @return \item{alpha}{Type I error for each gene comparison}


#' @export
powerCal <- function(ns, ms, vvmean1, FC, vvrho, hf, FDR) {

  vvmean2 <- vvmean1 * FC

  idx <- which(vvmean1!=vvmean2)
  K1 <- length(idx)
  K0 <- length(vvmean1) - K1
  vmean1 <- vvmean1[idx]
  vmean2 <- vvmean2[idx]
  vrho <- vvrho[idx]

  alpha <- alpha.FDR(ePower=0.8, FDR, K0, K1)

  eq1 <- function(bb) {
    alpha1 <- exp(bb)/(1+exp(bb))

    falsePositives <- K0*alpha1
    mp <- multiple.power2(ns, ms, vmean1, vmean2, vrho,
                          hf=hf, alpha=alpha1)
    falsePositives/(falsePositives + mp[2]) - FDR
  }

  if(is.na(eq1(log(alpha/(1-alpha))))) {
    pw <- c(NA, NA)
    alpha <- NA
  } else {
    root.solve <- nleqslv(x=log(alpha/(1-alpha)), fn=eq1, control=list(ftol=1e-6, maxit=15))

    if(root.solve$fvec > 0.005) {
      alpha <- NA
      pw <- NA
    } else {
      alpha <- exp(root.solve$x)/(1+exp(root.solve$x))
      pw <- (K0*alpha * (1 - (root.solve$fvec+FDR)))/(root.solve$fvec+FDR)/K1
    }

  }

  c(power=pw, alpha=alpha)

}



multiple.power2 <- function(ns, ms, vmean1, vmean2, vrho,
                            hf = function(x) sqrt(x),
                            alpha=0.05) {

  paras <- cbind(vmean1, vmean2, vrho)
  vr <- apply(paras, 1, function(x) {
    suppressWarnings(vs <- single.power3(ns=ns, ms=ms, mu1=x[1], mu2=x[2], rho=x[3],
                                         hf=hf, alpha=alpha))
  })

  c(K1=nrow(paras), enRejects=sum(vr))
}



single.power3 <- function(ns, ms, mu1, mu2, rho, ## faster
                          hf = function(x) sqrt(x),
                          alpha=0.05) {

  mus <- c(mu1, mu2)
  sigmas <- hf(mus)

  csum.inv.Ms <- 1/(1-rho) * (1 - ns*(rho/(1-rho+ns*rho)) )

  DVD.list <- lapply(c(0,1), function(x) {

    n <- ns[x+1]
    csum.inv.M <- csum.inv.Ms[x+1]

    mu <- mus[x+1]
    sigma <- sigmas[x+1]

    csum.inv.V <- sigma^-2 * csum.inv.M

    v1 <- mu   * csum.inv.V
    v2 <- mu*x * csum.inv.V

    ms[x+1] * n * mu * rbind( c(v1, v1 * x),
                              c(v2, v2 * x))

  })

  list_sum <- Reduce("+", DVD.list)
  SIGMA <- solve(list_sum)

  effect.size <- (log(mu2)-log(mu1))/sqrt(SIGMA[2,2])
  q.alpha <- qt(1-alpha/2, df=sum(ms)-2)
  power1 <- pt(q.alpha - effect.size, df=sum(ms)-2, lower.tail = FALSE) +
    pt(-q.alpha - effect.size, df=sum(ms)-2)
  as.numeric(power1)

}



alpha.FDR <- function(ePower, FDR, K0, K1) {
  (ePower * K1)*FDR/(K0*(1-FDR))
}

