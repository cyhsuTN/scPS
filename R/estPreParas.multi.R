#' Estimate Required Parameters from a Pilot Data (Multiple Cell Types)
#' @description A function to estimate required parameters from a pilot data for scPS.

#' @param counts Normalized count matrix.
#' @param cell.info Cell information.
#' @param id Sample/subject ID.
#' @param x1 Groups/conditions/stages.
#' @param cellcluster cell types.
#' @param cells.interesting Cell types of interest.

#' @examples # https://github.com/cyhsuTN/scPS
#' @examples # load(file = "DataForDemo/COVID19n.rda")
#' @examples # counts <- COVID19n$counts
#' @examples # cell.info <- COVID19n$cell.info
#' @examples # geneObject <- estPreParas.multi(counts, cell.info,
#' @examples #                     id="SampleId", x1="condition",
#' @examples #                     cells.interesting=c("T cells", "DC", "Prolif.T")[c(2,3)])

#' @examples # https://github.com/cyhsuTN/scPS
#' @examples # load(file = "DataForDemo/GSE120575n.rda")
#' @examples # counts <- GSE120575n$counts
#' @examples # cell.info <- GSE120575n$cell.info
#' @examples # cell.info$TX <- factor(cell.info$TX, levels = c("Pre", "Post"))
#' @examples # geneObject <- estPreParas.multi(counts, cell.info, id="ptID", x1="TX",
#' @examples #                  cells.interesting=c("CD8", "CD4", "Macrophage", "NK", "B")[4:5])
#'


#' @export
estPreParas.multi <- function(counts, cell.info,
                              id="SampleId",
                              x1="condition",
                              cellcluster="cellcluster",
                              cells.interesting=NULL) {


  if(is.factor(cell.info[,x1])) {
    groups <- levels(cell.info[,x1])
  } else {
    stop("Please set x1 of cell.info as factor with a reference group")
  }

  cell.info$id <- cell.info[,id]
  cell.info$x1 <- ifelse(cell.info[,x1]==groups[1], 0, 1)
  cell.info$cellcluster <- cell.info[,cellcluster]

  if(is.null(cells.interesting)) cells.interesting <-
    names(table(cell.info$cellcluster))[order(table(cell.info$cellcluster), decreasing = T)]

  if (length(groups)==1) {
    print("One group")
  } else {
    BA <- any(apply(table(cell.info$id, cell.info$x1), 1, function(x) all(x>0, na.rm = T)))
    if(BA) print("Paired groups") else print("Independent two groups")
  }

  geneInfo <- lapply(cells.interesting, function(cc) {
    counts.small <- counts[,cell.info$cellcluster==cc]
    cell.info.small <- cell.info[cell.info$cellcluster==cc, c("id", "x1")]
    aa <- estPreParas(counts=counts.small, cell.info=cell.info.small)
    aa
  })

  names(geneInfo) <- cells.interesting
  list(counts=counts,
       cell.info=cell.info,
       geneInfo=geneInfo)

}


