#' Select Candidate Genes and DEGs
#' @description A function to select candidate genes and DEGs for power and sample size calculation.

#' @param geneObject Object of estPreParas.multi.
#' @param cells.interesting Cell types of interest.
#' @param thres.nonZeroPs A threshold of non-zero proportion to exclude non-expressed genes.
#' @param nGenesCandidate Number of candidate genes (pre-specified). Genes with large observed
#' fold-changes are selected as candidate genes.
#' @param propDEGs A proportion of DEGs among candidate genes (pre-specified).
#' @param FC.way FC.way = "median" (default) assumes same logFC = \emph{d}*sign(logFC)
#' for each DEG, where \emph{d} is the median of |logFC|'s among DEGs;
#' FC.way = "individual" assumes logFC of each DEG = its observed log ratio of means between groups.
#' @param cg.name Names of candidate genes if users want to select candidate genes by themselves.
#' The default = NULL, meaning genes with large observed fold-changes are selected as candidate genes.
#' @param ... Arguments in fitSdmean: \emph{sp.method}, \emph{probs}, and \emph{span}.

#' @return \item{logFC.cg}{Observed logFCs for candidate genes}
#' @return \item{alogFC.deg}{Median of |logFC|'s among DEGs}
#' @return \item{FC}{Assumed FCs for power and sample size calculation}
#' @return \item{icc}{ICCs for candidate genes}
#' @return \item{hf}{SD-Mean relationship function}
#' @return \item{idx.cg}{Indices of candidate genes}
#' @return \item{idx.deg}{Indices of DEGs}

#' @examples # https://github.com/cyhsuTN/scPS
#' @examples # load(file = "DataForDemo/COVID19n.rda")
#' @examples # counts <- COVID19n$counts
#' @examples # cell.info <- COVID19n$cell.info
#' @examples # geneObject <- estPreParas.multi(counts, cell.info,
#' @examples #                     id="SampleId", x1="condition",
#' @examples #                     cells.interesting=c("T cells", "DC", "Prolif.T")[c(2,3)])
#' @examples # Genes.tested <- geneCandidate(geneObject, FC.way="median")

#' @examples # https://github.com/cyhsuTN/scPS
#' @examples # load(file = "DataForDemo/GSE120575n.rda")
#' @examples # counts <- GSE120575n$counts
#' @examples # cell.info <- GSE120575n$cell.info
#' @examples # cell.info$TX <- factor(cell.info$TX, levels = c("Pre", "Post"))
#' @examples # geneObject <- estPreParas.multi(counts, cell.info, id="ptID", x1="TX",
#' @examples #                  cells.interesting=c("CD8", "CD4", "Macrophage", "NK", "B")[4:5])
#' @examples # Genes.tested <- geneCandidate(geneObject, FC.way="median")


#' @export
geneCandidate <- function(geneObject, cells.interesting=NULL,
                          thres.nonZeroPs = 0.1,
                          nGenesCandidate = 2000,
                          propDEGs = 0.01,
                          FC.way = c("median", "individual")[1],
                          cg.name=NULL,
                          ...
                          ) {

  counts <- geneObject$counts
  cell.info <- geneObject$cell.info

  if(is.null(cells.interesting)) cells.interesting <- names(geneObject$geneInfo)

  idx.id <- grep("id", colnames(geneObject$cell.info))
  idx.x1 <- grep("x1", colnames(geneObject$cell.info))

  BA <- any(apply(table(cell.info[,idx.id], cell.info[,idx.x1]), 1, function(x) all(x>0, na.rm = T)))
  if(BA) print("Paired-group comparison") else print("Independent two-group comparison")

  Genes.tested <- lapply(cells.interesting, function(cc) {

    counts.small <- counts[,cell.info$cellcluster==cc]
    cell.info.small <- cell.info[cell.info$cellcluster==cc, c(idx.id, idx.x1)] #c("id", "x1")

    aa <- geneObject$geneInfo[[cc]]

    if(is.null(cg.name)) {
      idx.selected <- which(aa$nonZeroPs1>thres.nonZeroPs & aa$nonZeroPs2>thres.nonZeroPs)

      logFC <- log(aa$means2[idx.selected]/aa$means1[idx.selected])

      idx.cg <- head(idx.selected[order(abs(logFC), decreasing = T)], nGenesCandidate)
    } else {
      idx.cg <- which(rownames(counts) %in% cg.name)
    }

    aaP <- aa[idx.cg,]

    curve.view <- fitSdmean(preParas=aaP, cell.name=cc, ...)
    hf <- curve.view$hf.sigma

    if(BA) {
      geeout <- geeDEA.BA(counts.test=counts.small[idx.cg,], cell.info=cell.info.small, hf=hf, BCC=1)
    } else {
      geeout <- geeDEA(counts.test=counts.small[idx.cg,], cell.info=cell.info.small, hf=hf, BCC=1)
    }

    mean.control <- aaP$means1
    logFC <- log(aaP$means2/aaP$means1)
    idx.DEG <- head(order(geeout$unadj.p), propDEGs*length(mean.control))
    s.logFC <- summary(abs(logFC[idx.DEG]))

    FC <- rep(1, length(mean.control))
    if(FC.way=="individual") {
      FC[idx.DEG] <- exp(logFC[idx.DEG])
    } else {
      FC[idx.DEG] <- ifelse(logFC[idx.DEG]>0, round(exp(s.logFC[3]),1), 1/round(exp(s.logFC[3]),1))
    }

    list(pc=round(ncol(counts.small)/ncol(counts),3),
         K=length(mean.control),
         mean.control=mean.control,
         logFC.cg=logFC,  # Observed logFCs for candidate genes
         alogFC.deg=s.logFC, # Median of |logFC|'s among DEGs
         FC=FC,        # Assumed FCs for power and sample size calculation
         icc=aaP$icc,
         hf=hf,
         idx.cg=idx.cg,
         idx.deg=idx.cg[idx.DEG])
  })

  names(Genes.tested) <- cells.interesting
  Genes.tested

}
