#' Optimal Combination of m and n
#' @description A function to show optimal m and n combination given a cost function.

#' @param m.n.power Object of sizeCal/sizeCal.multi or sizeCal.BA/sizeCal.multi.BA.
#' @param costfun Cost function. For example, costfun = function(m, n) m*n.
#' @param ePower  Expected overall power.
#' @param budget  Budget. The default = NULL.
#' @param prop.this.cell Proportion of this cell type of interest. The default = NULL.

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
#' @examples # head(optimalCost(m.n.power))

#' @export
optimalCost <- function(m.n.power, costfun=function(m, n) 0.24*m*n,
                        ePower=0.8, budget=NULL, prop.this.cell=NULL) {

  m.n <- m.n.power$m.n.power
  rc  <- m.n.power$rc
  pw  <- m.n.power$m.n.power$power
  tol <- m.n.power$total

  if(is.null(prop.this.cell)) prop.this.cell <- 1

  rm.idx <- which(colnames(m.n) %in% c("m", "n"))

  if(m.n.power$BA) {

    m <- m.n$m
    n1 <- m.n$n
    n2 <- round(m.n$n * rc)

    cost <- costfun(m, n1+n2)/prop.this.cell
    m.n <- cbind(cost, m, n1, n2, m.n[,-rm.idx, drop=F])

  } else {

    n1 <- m.n$n
    n2 <- round(m.n$n * rc)

    m1 <- m.n$m

    if (is.null(tol)) {
      r   <- m.n.power$r
      m2 <- round(m.n$m * r)
    } else {
      m2 <- tol - m.n$m
    }

    cost <- (costfun(m1, n1) + costfun(m2, n2))/prop.this.cell
    m.n <- cbind(cost, m1, m2, n1, n2, m.n[,-rm.idx, drop=F])
  }

  Rank <- 1:nrow(m.n)
  if(is.null(budget)) {
    tf.pw <- as.numeric(1*(pw>=ePower))
    out <- cbind(Rank, m.n[order(tf.pw, cost, pw, decreasing = c(T, F, T)),])
    out <- out[out$power>=ePower,]
  } else {
    tf.pw <- as.numeric(1*(m.n$cost<=budget))
    out <- cbind(Rank, m.n[order(tf.pw, pw, cost, decreasing = c(T, T, F)),])
    out <- out[out$cost<=budget,]
  }

  out

}
