#' Fit a Curve between Gene Means and Gene Standard Deviations
#' @description A function to fit a curve between gene means and gene standard deviations.
#' @import splines
#' @importFrom stats fitted lm pnorm predict loess
#' @importFrom graphics lines

#' @param preParas estPreParas object.
#' @param sp.method Spline methods: "bs" (cubic B-spline), "ns" (natural cubic spline),
#' and "loess" (LOESS smoothing). Default = "ns".
#' @param probs  Inter quantile that will be used as knots in "bs" and "ns" methods.
#' @param span A parameter which controls the degree of smoothing if sp.method = "loess".
#' @param cell.name cell name.

#' @examples # https://github.com/cyhsuTN/scPS
#' @examples # load(file = "DataForDemo/COVID19n.rda")
#' @examples # counts <- COVID19n$counts
#' @examples # cell.info <- COVID19n$cell.info
#' @examples # cells.interesting <- which(cell.info$cellcluster=="Prolif.T")
#' @examples # counts.small <- counts[,cells.interesting]
#' @examples # cell.info.small <- data.frame(id=cell.info[cells.interesting, "SampleId"],
#' @examples #        x1=ifelse(cell.info[cells.interesting, "condition"]=="Neg", 0, 1) )
#' @examples # aa <- estPreParas(counts=counts.small, cell.info=cell.info.small)
#'
#' @examples # ## Select 2000 candidate genes
#' @examples # idx.selected <- which(aa$nonZeroPs1>0.1 & aa$nonZeroPs2>0.1)
#' @examples # logFC <- log(aa$means2[idx.selected]/aa$means1[idx.selected])
#' @examples # idx.cg <- head(idx.selected[order(abs(logFC), decreasing = T)],2000)
#' @examples # aaP <- aa[idx.cg,]
#' @examples # curve.view <- fitSdmean(preParas=aaP, sp.method=c("bs","ns","loess")[2])
#' @examples # hf <- curve.view$hf.sigma


#' @export
fitSdmean <- function(preParas, sp.method="ns", probs=c(0.33, 0.67),
                      span=0.75, cell.name=NULL) {

  xx <- as.numeric(unlist(preParas[,grep("mean", colnames(preParas))]))
  yy <- sqrt(as.numeric(unlist(preParas[,grep("var", colnames(preParas))])))

  if(sp.method=="bs") {
    xx <- log(xx)
    fit.spline <- lm(log(yy) ~ bs(xx, knots = quantile(xx, probs )))
  } else if(sp.method=="loess") {
    xx <- log(xx)
    fit.spline <- loess(log(yy) ~ xx,
                        span = span)
  } else {
    xx <- log(xx)
    fit.spline <- lm(log(yy) ~ ns(xx, knots = quantile(xx, probs )))
  }


  plot(xx, log(yy), xlab = "log(gene mean)", ylab = "log(standard deviation)",
       main = paste0("SD-Mean curve - ", cell.name),
       cex.lab=1.5, cex.axis=1.5)
  lines(xx[order(xx)], fitted(fit.spline)[order(xx)], col="red", lwd=2)

  hf.sigma <- function(uu) {
    xx <- log(uu)
    exp(predict(fit.spline, data.frame(xx=xx)))
  }

  list(hf.sigma=hf.sigma, fit.spline=fit.spline)

}
