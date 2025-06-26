#' Power Calculation for DEA of Multiple Cell Types between Paired Groups
#' @description A function to calculate powers for DEA of multiple cell types between paired groups.
#' @importFrom nleqslv nleqslv

#' @param ns Numbers of total cells in two stages per subject.
#' @param m Number of subjects.
#' @param Genes.tested Object of estPreParas.multi.
#' @param FDR False discovery rate.

#' @examples # https://github.com/cyhsuTN/scPS
#' @examples # load(file = "DataForDemo/GSE120575n.rda")
#' @examples # counts <- GSE120575n$counts
#' @examples # cell.info <- GSE120575n$cell.info
#' @examples # cell.info$TX <- factor(cell.info$TX, levels = c("Pre", "Post"))
#' @examples # geneObject <- estPreParas.multi(counts, cell.info, id="ptID", x1="TX",
#' @examples #                     cells.interesting=c("NK", "B"))
#' @examples # Genes.tested <- geneCandidate(geneObject)
#' @examples # powerCal.multi.BA(ns=c(1,1)*550, m=11, Genes.tested, FDR=0.05)

#' @export
powerCal.multi.BA <- function(ns=c(600,600), m=10, Genes.tested, FDR) {

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
      mp <- multiple.power.BA(ns=vnc[[cj]], m, vmean1[[cj]], vmean2[[cj]], vrho[[cj]],
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
        mp <- multiple.power.BA(ns=vnc[[cj]], m, vmean1[[cj]], vmean2[[cj]], vrho[[cj]],
                                hf=vhf[[cj]], alpha=alpha)
        mp[2]/mp[1]
      })

    }

  }

  out1 <- c(power=pw, alpha=alpha, sep.power=as.numeric(vmp1))
  names(out1)[-c(1,2)] <- paste0("power.",names(Genes.tested))
  out1

}
