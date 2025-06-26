#' Power Calculation for DEA of a Single Cell Type between Paired Groups
#' @description A function to calculate powers for DEA of a single cell type between paired groups.
#' @importFrom nleqslv nleqslv
#' @importFrom stats pt qt

#' @param ns Numbers of cells of interest in two stages per subject.
#' @param m Number of subjects.
#' @param vvmean1 Gene means of candidate genes in the pre-stage.
#' @param FC Fold changes of candidate genes between two stages.
#' @param vvrho Cell-cell correlations of candidate genes within subjects.
#' @param hf A curve showing the relationship between gene means and gene standard deviations.
#' @param FDR False discovery rate.

#' @examples set.seed(12345)
#' @examples vvmean1 <- rgamma(1000, shape=2, scale=0.5)
#' @examples FC <- c(rep(2, 50), rep(1, 950))
#' @examples ab <- gammaTrans(mean=0.01, q95=0.1) # Output shape and scale parameters
#' @examples vvrho <- rgamma(1000, shape=ab[1], scale=ab[2])
#' @examples hf <- function(x) sqrt(x*(1+3*x))
#' @examples powerCal.BA(ns=c(1,1)*70, m=8, vvmean1, FC, vvrho, hf, FDR=0.05)

#' @return \item{power}{Statistical power}
#' @return \item{alpha}{Type I error for each gene comparison}


#' @export
powerCal.BA <- function(ns, m, vvmean1, FC, vvrho, hf, FDR) {

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
    mp <- multiple.power.BA(ns, m, vmean1, vmean2, vrho,
                            hf=hf, alpha=alpha1)
    falsePositives/(falsePositives + mp[2]) - FDR
  }

  if(is.na(eq1(log(alpha/(1-alpha))))) {
    pw <- c(NA, NA)
    alpha <- NA
  } else {
    (root.solve <- nleqslv(x=log(alpha/(1-alpha)), fn=eq1, control=list(ftol=1e-6, maxit=15)))

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


multiple.power.BA <- function(ns, m, vmean1, vmean2, vrho,
                              hf = function(x) sqrt(x),
                              alpha=0.05) {

  paras <- cbind(vmean1, vmean2, vrho)
  vr <- apply(paras, 1, function(x) {

    suppressWarnings(vs <- single.power.BA(ns=ns, m=m, mu1=x[1], mu2=x[2], rho=x[3],
                                           hf=hf, alpha=alpha))
  })

  c(K1=nrow(paras), enRejects=sum(vr))
}


single.power.BA <- function(ns, m, mu1, mu2, rho,
                            hf = function(x) sqrt(x),
                            alpha=0.05) {

  mus <- c(mu1, mu2)
  sigmas <- hf(mus)

  mu.by.sigma <- mus/sigmas

  n1 <- ns[1]
  n2 <- ns[2]

  tp1 <- 1/(1-rho)
  tp2 <- rho/(1-rho+(n1+n2)*rho)

  csum.inv.M1 <- tp1 * (1 - n1*tp2 )
  csum.inv.M2 <- tp1 * (  - n2*tp2 )
  csum.inv.M3 <- tp1 * (  - n1*tp2 )
  csum.inv.M4 <- tp1 * (1 - n2*tp2 )

  v11 <- sum(mu.by.sigma * c(csum.inv.M1, csum.inv.M2))
  v12 <- sum(mu.by.sigma * c(csum.inv.M3, csum.inv.M4))
  v21 <- sum(c(0,mu.by.sigma[2]) * c(csum.inv.M1, csum.inv.M2))
  v22 <- sum(c(0,mu.by.sigma[2]) * c(csum.inv.M3, csum.inv.M4))

  list_sum <- m * rbind(c(sum(ns * c(v11, v12) * mu.by.sigma), sum(ns * c(v11, v12) * c(0,mu.by.sigma[2]))),
                        c(sum(ns * c(v21, v22) * mu.by.sigma), sum(ns * c(v21, v22) * c(0,mu.by.sigma[2]))) )
  SIGMA <- solve(list_sum)

  effect.size <- (log(mu2)-log(mu1))/sqrt(SIGMA[2,2])
  q.alpha <- qt(1-alpha/2, df=m-2)
  power1 <- pt(q.alpha - effect.size, df=m-2, lower.tail = FALSE) +
    pt(-q.alpha - effect.size, df=m-2)
  as.numeric(power1)

}


alpha.FDR <- function(ePower, FDR, K0, K1) {
  (ePower * K1)*FDR/(K0*(1-FDR))
}
