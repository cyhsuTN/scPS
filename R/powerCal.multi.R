#' Power Calculation for DEA of Multiple Cell Types between Two Independent Groups
#' @description A function to calculate powers for DEA of multiple cell types
#' between two independent groups.
#' @importFrom nleqslv nleqslv

#' @param ns Number of total cells required per subject in control and experimental groups.
#' @param ms Numbers of subjects in control and experimental groups, respectively.
#' @param Genes.tested Object of estPreParas.multi.
#' @param FDR False discovery rate.

#' @examples # https://github.com/cyhsuTN/scPS
#' @examples # load(file = "DataForDemo/COVID19n.rda")
#' @examples # counts <- COVID19n$counts
#' @examples # cell.info <- COVID19n$cell.info
#' @examples # geneObject <- estPreParas.multi(counts, cell.info,
#' @examples #                     id="SampleId", x1="condition",
#' @examples #                     cells.interesting=c("T cells", "DC", "Prolif.T")[c(2,3)])
#' @examples # Genes.tested <- geneCandidate(geneObject)
#' @examples # powerCal.multi(ns=c(400,400), ms=c(1,1)*11, Genes.tested, FDR=0.05)

#' @export
powerCal.multi <- function(ns, ms, Genes.tested, FDR) {


  idx <- lapply(1:length(Genes.tested), function(cj)  which(Genes.tested[[cj]]$FC!=1) )
  vK1 <- lapply(1:length(Genes.tested), function(cj) length(idx[[cj]]))
  vK0 <- lapply(1:length(Genes.tested), function(cj) Genes.tested[[cj]]$K - vK1[[cj]])
  vmean1 <- lapply(1:length(Genes.tested), function(cj) Genes.tested[[cj]]$mean.control[idx[[cj]]])
  vmean2 <- lapply(1:length(Genes.tested), function(cj) vmean1[[cj]] * Genes.tested[[cj]]$FC[idx[[cj]]]   )
  vrho <- lapply(1:length(Genes.tested), function(cj) Genes.tested[[cj]]$icc[idx[[cj]]]   )
  vhf <- lapply(1:length(Genes.tested), function(cj) Genes.tested[[cj]]$hf)
  vnc <- lapply(1:length(Genes.tested), function(cj) round(Genes.tested[[cj]]$pc * ns) )

  K0 <- sum(unlist(vK0))
  K1 <- sum(unlist(vK1))

  alpha <- alpha.FDR(ePower=0.8, FDR, K0, K1)

  eq1 <- function(bb) {
    alpha1 <- exp(bb)/(1+exp(bb))

    falsePositives <- K0*alpha1
    vmp2 <- sapply(1:length(Genes.tested), function(cj) {
      mp <- multiple.power2(ns=vnc[[cj]], ms, vmean1[[cj]], vmean2[[cj]], vrho[[cj]],
                            hf=vhf[[cj]], alpha=alpha1)
      mp[2]
    })

    falsePositives/(falsePositives + sum(vmp2) ) - FDR
  }

  if(is.na(eq1(log(alpha/(1-alpha))))) {
    pw <- NA
    alpha <- NA
    vmp1 <- rep(NA, length(Genes.tested))
  } else {
    (root.solve <- nleqslv(x=log(alpha/(1-alpha)), fn=eq1, control=list(ftol=1e-6, maxit=15)))

    if(root.solve$fvec > 0.005) {
      alpha <- NA
      pw <- NA
      vmp1 <- rep(NA, length(Genes.tested))
    } else {
      alpha <- exp(root.solve$x)/(1+exp(root.solve$x))

      pw <- (K0*alpha * (1 - (root.solve$fvec+FDR)))/(root.solve$fvec+FDR)/K1

      vmp1 <- sapply(1:length(Genes.tested), function(cj) {
        mp <- multiple.power2(ns=vnc[[cj]], ms, vmean1[[cj]], vmean2[[cj]], vrho[[cj]],
                              hf=vhf[[cj]], alpha=alpha)
        mp[2]/mp[1]
      })

    }

  }

  out1 <- c(power=pw, alpha=alpha, sep.power=as.numeric(vmp1))
  names(out1)[-c(1,2)] <- paste0("power.",names(Genes.tested))
  out1

}
