#' Figure Showing Numbers of Cells Required by Cell Types of Interest
#' @description A function to plot numbers of cells required by cell types of interest.
#' @import ggplot2
#' @importFrom ggpubr ggarrange

#' @param m.n.power Object of sizeCal.multi or sizeCal.multi.BA.
#' @param ct.names Names of cell types.
#' @param ePower  Expected overall power (expected proportion of DEGs detected).
#' @param nrow  Number of rows to show figures.
#' @param ncol  Number of columns to show figures.

#' @examples # https://github.com/cyhsuTN/scPS
#' @examples # load(file = "DataForDemo/COVID19n.rda")
#' @examples # counts <- COVID19n$counts
#' @examples # cell.info <- COVID19n$cell.info
#' @examples # geneObject <- estPreParas.multi(counts, cell.info,
#' @examples #                     id="SampleId", x1="condition",
#' @examples #                     cells.interesting=c("T cells", "DC", "Prolif.T")[c(2,3)])
#' @examples # Genes.tested <- geneCandidate(geneObject)
#' @examples # m.n.power <- sizeCal.multi(low.up.m=c(10,14), low.up.n=c(200,500),
#' @examples #     ePower=0.8, FDR=0.05, grid.m=1, grid.n=50, r=1, rc=1, total=NULL, Genes.tested)
#' @examples # plotPower.sep(m.n.power)

#' @export
plotPower.sep <- function(m.n.power, ct.names=NULL, ePower=0.8,
                          nrow = NULL, ncol = NULL) {

  if(is.null(ct.names)) ct.names <- names(m.n.power$pc)

  plotlist <- lapply(ct.names, function(ct.name) {
    pc <- m.n.power$pc
    K  <- m.n.power$K
    dat2 <- m.n.power$m.n.power

    ct.index <- which(names(pc)==ct.name)
    dat2$nc <- round(dat2$n * pc[ct.index])
    dat2$power <- dat2[,grep("power.", colnames(dat2))[ct.index]]

    color1 <- ifelse(dat2$power > ePower, "blue", "red")
    value1 <- c("blue", "red")
    if(all(dat2$power < ePower, na.rm = T)) {
      color1 <- "red"; value1 <- "red"
    }

    if(m.n.power$BA) {
      xxlab <- "No. of subjects"
      yylab <- "No. of cells of interest per subject in pre-stage"
    } else {
      xxlab <- "No. of subjects in control"
      yylab <- "No. of cells of interest per subject in control"
    }


    fig <- ggplot(data=NULL, aes(x=dat2$m, y=dat2$nc, fill=dat2$power)) +
      geom_point(size=10, shape=21, colour = "transparent") +
      geom_text(aes(label = round(dat2$power, 2), color = color1, fontface=2),
                size = 3.2, show.legend = FALSE) +
      scale_color_manual(values = value1) +
      scale_fill_gradient(low = "yellow", high = "green") +
      scale_x_continuous(breaks = dat2$m) +
      scale_y_continuous(breaks = dat2$nc) +
      theme_minimal() +
      theme(axis.text=element_text(size=12),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_text(size = 14)) +
      xlab(xxlab) +
      ylab(yylab) +
      ggtitle(paste0(names(pc)[ct.index],
                     " (prop. = ",round(pc[ct.index]*100,1),"%), K = ",
                     K[ct.index]))

    fig$labels$fill <- "power"

    fig
  })

  if(is.null(nrow)) nrow <- 1
  if(is.null(ncol)) ncol <- length(ct.names)
  ggarrange(plotlist=plotlist, nrow = nrow,  ncol = ncol)

}
