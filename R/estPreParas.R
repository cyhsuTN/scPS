#' Estimate Required Parameters from a Pilot Data (a Single Cell Type)
#' @description A function to estimate required parameters from a pilot data for scPS.
#' @import Matrix
#' @importFrom stats sd

#' @param counts Normalized count matrix.
#' @param cell.info Cell information.

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

#' @export
estPreParas <- function(counts, cell.info) {

  cell.info$id <- as.character(cell.info$id)

  oout <- t(apply(counts, 1, info, id=cell.info$id, x1=cell.info$x1))

  return(data.frame(oout))
}



info <- function(y, id, x1) {

  means <- tapply(y, x1, mean)
  sigmas <- tapply(y, x1, sd)

  nonZeroPs <- tapply(y, x1, function(y) mean(y!=0))

  if(length(means)==1) means <- rep(means,2)
  if(length(sigmas)==1) sigmas <- rep(sigmas,2)

  x1keep <- !(x1 %in% (which(sigmas==0) - 1))
  y  <-  y[x1keep]
  id <- id[x1keep]
  x1 <- x1[x1keep]

  if(all(sigmas==0)) {
    icc <- NA
  } else {
    if(sigmas[1]==0 & sigmas[2]!=0) {
      resi <- (y-means[2])/sigmas[2]
    } else if(sigmas[1]!=0 & sigmas[2]==0) {
      resi <- (y-means[1])/sigmas[1]
    } else {
      resi <- ifelse(x1==0, (y-means[1])/sigmas[1], (y-means[2])/sigmas[2])
    }

    np <- length(unique(x1))
    e_sum <- sum(tapply(resi, id, function(e) sum(e)^2 - sum(e^2)))
    n_sum <- sum(tapply(resi, id, function(e) length(e)*(length(e)-1))) - 2 * np

    icc <- (e_sum/n_sum)
  }

  if(length(means)==1) {
    c(icc=icc, mean=as.numeric(means)[1], var=as.numeric((sigmas^2))[1], nonZeroPs=as.numeric(nonZeroPs)[1])
  } else {
    c(icc=icc, means=as.numeric(means), vars=as.numeric(sigmas^2), nonZeroPs=as.numeric(nonZeroPs))
  }

}
