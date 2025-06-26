#' Sample Size Calculation for DEA of Multiple Cell Types between Paired Groups
#' @description A function to calculate sample sizes and numbers of total cells
#' between paired groups.
#' @import ggplot2
#' @importFrom nleqslv nleqslv
#' @importFrom stats pt qt


#' @param low.up.m Lower and upper bounds of subjects required.
#' @param low.up.n Lower and upper bounds of total cells required in the pre-stage per subject.
#' @param ePower  An expected overall power (an expected proportion of DEGs detected).
#' @param FDR False discovery rate.
#' @param grid.m An integer by which the \emph{low.up.m} is divided.
#' @param grid.n An integer by which the \emph{low.up.n} is divided.
#' @param rc A ratio of no. of cells in post-stage to no. of cells in pre-stage.
#' @param Genes.tested Object of estPreParas.multi.

#' @examples # https://github.com/cyhsuTN/scPS
#' @examples # load(file = "DataForDemo/GSE120575n.rda")
#' @examples # counts <- GSE120575n$counts
#' @examples # cell.info <- GSE120575n$cell.info
#' @examples # cell.info$TX <- factor(cell.info$TX, levels = c("Pre", "Post"))
#' @examples # geneObject <- estPreParas.multi(counts, cell.info, id="ptID", x1="TX",
#' @examples #                     cells.interesting=c("NK", "B"))
#' @examples # Genes.tested <- geneCandidate(geneObject)
#' @examples # m.n.power <- sizeCal.multi.BA(low.up.m=c(10,14), low.up.n=c(400,700),
#' @examples #     ePower=0.8, FDR=0.05, grid.m=1, grid.n=50, rc=1, Genes.tested)
#' @examples # m.n.power$fig
#'

#' @export
sizeCal.multi.BA <- function(low.up.m=c(2,10), low.up.n=c(10,100), ePower=0.8, FDR=0.05,
                             grid.m=2, grid.n=10, rc=1,
                             Genes.tested) {

  a <- seq(low.up.m[1], low.up.m[2], grid.m)
  b <- seq(low.up.n[1], low.up.n[2], grid.n)

  comb.m.n <- expand.grid(a,b)

  z <- apply(comb.m.n, 1, function(x) {
    m <- x[1]
    ns <- c(x[2], round(x[2]*rc) )

    powerCal.multi.BA(ns, m, Genes.tested, FDR)
  })

  dat2 <- data.frame(cbind(m=comb.m.n[,1], n=comb.m.n[,2], t(z)))


  color1 <- ifelse(dat2$power > ePower, "blue", "red")
  value1 <- c("blue", "red")
  if(all(dat2$power < ePower)) {
    color1 <- "red"; value1 <- "red"
  }

  fig <- ggplot(data=NULL, aes(x=dat2$m, y=dat2$n, fill=dat2$power)) +
    geom_point(size = 10, shape=21, colour = "transparent") +
    geom_text(aes(label = round(dat2$power, 2), color = color1, fontface=2),
              size = 3.2, show.legend = FALSE) +
    scale_color_manual(values = value1) +
    scale_fill_gradient(low = "yellow", high = "green") +
    scale_x_continuous(breaks = dat2$m) +
    scale_y_continuous(breaks = dat2$n) +
    theme_minimal() +
    theme(axis.text=element_text(size=12),
          axis.title.y = element_text(size = 12),
          axis.title.x = element_text(size = 14)) +
    xlab("No. of subjects") +
    ylab("No. of total cells per subject in pre-stage") +
    ggtitle(paste0("Cells Tested (", paste0(names(Genes.tested), collapse = ", "), ": ",
                   round(sum(sapply(Genes.tested, function(x) x$pc))*100,1),"%)"))

  fig$labels$fill <- "power"

  pc <- sapply(Genes.tested, function(x) x$pc)
  K  <- sapply(Genes.tested, function(x) x$K)

  list(m.n.power=dat2, fig=fig, pc=pc, K=K, BA=TRUE, rc=rc)

}
