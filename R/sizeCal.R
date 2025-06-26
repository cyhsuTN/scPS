#' Sample Size Calculation for DEA of a Single cell Type between Two Independent Groups
#' @description A function to calculate sample sizes and numbers of cells of interest
#' between two independent groups.
#' @import ggplot2
#' @importFrom nleqslv nleqslv
#' @importFrom stats pt qt


#' @param low.up.m Lower and upper bounds of subjects in control.
#' @param low.up.n Lower and upper bounds of cells of interest in each subject of control.
#' @param ePower  An expected overall power (an expected proportion of DEGs detected).
#' @param FDR False discovery rate.
#' @param grid.m An integer by which the \emph{low.up.m} is divided.
#' @param grid.n An integer by which the \emph{low.up.n} is divided.
#' @param r An allocation ratio (a ratio of no. of experimental to no. of control).
#' @param rc A ratio of no. of cells in experimental to in control.
#' @param total A total of sample size.
#' @param vvmean1 Gene means of candidate genes in the control group.
#' @param FC Fold changes of candidate genes between the control and the experimental groups.
#' @param vvrho Cell-cell correlations of candidate genes within subjects.
#' @param hf A curve showing the relationship between gene means and gene deviations.

#' @examples set.seed(12345)
#' @examples vvmean1 <- rep(1, 1000)
#' @examples FC <- c(rep(2, 50), rep(1, 950))
#' @examples ab <- gammaTrans(mean=0.01, q95=0.1) # Output shape and scale parameters
#' @examples vvrho <- rgamma(1000, shape=ab[1], scale=ab[2])
#' @examples hf <- function(x) sqrt(x*(1+3*x))
#' @examples sizeCal(low.up.m=c(7,11), low.up.n=c(30,60), ePower=0.8, FDR=0.05,
#' @examples        grid.m=1, grid.n=5, r=1, rc=1, total=NULL, vvmean1, FC, vvrho, hf)

#' @return \item{m.n.power}{A dataframe showing powers}
#' @return \item{fig}{A heat map}


#' @export
sizeCal <- function(low.up.m=c(6,10), low.up.n=c(10,100), ePower=0.8, FDR=0.05,
                     grid.m=1, grid.n=10, r=1, rc=1, total=NULL,
                     vvmean1, FC, vvrho,
                     hf=function(x) sqrt(x)) {

  if(!is.null(total)) low.up.m[2] <- ifelse(low.up.m[2]>total, total - 1, low.up.m[2])

  a <- seq(low.up.m[1], low.up.m[2], grid.m)
  b <- seq(low.up.n[1], low.up.n[2], grid.n)

  comb.m.n <- expand.grid(a,b)

  z <- apply(comb.m.n, 1, function(x) {
    ms <- c(x[1], ifelse(!is.null(total), total-x[1], round(x[1]*r) ) )
    ns <- c(x[2], round(x[2]*rc) )
    powerCal(ns=ns, ms=ms, vvmean1, FC, vvrho, hf, FDR)[1]
  })

  dat2 <- data.frame(m=comb.m.n[,1], n=comb.m.n[,2], power=z)

  color1 <- ifelse(dat2$power > ePower, "blue", "red")
  value1 <- c("blue", "red")
  if(all(dat2$power < ePower, na.rm = T)) {
    color1 <- "red"; value1 <- "red"
  }

  fig <- ggplot(data=NULL, aes(x=dat2$m, y=dat2$n, fill=dat2$power)) +
    geom_point(size=10, shape=21, colour = "transparent") +
    geom_text(aes(label = round(dat2$power, 2), color = color1, fontface=2),
              size = 3.2, show.legend = FALSE) +
    scale_color_manual(values = value1) +
    scale_fill_gradient(low = "yellow", high = "green") +
    scale_x_continuous(breaks = dat2$m) +
    scale_y_continuous(breaks = dat2$n) +
    theme_minimal() +
    theme(axis.text=element_text(size=12),
          axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14)) +
    xlab("No. of subjects in control") +
    ylab("No. of cells of interest per subject in control")

  fig$labels$fill <- "power"

  K <- length(vvmean1)
  list(m.n.power=dat2, fig=fig, pc=1, K=K, BA=FALSE, r=r, rc=rc, total=total)

}



