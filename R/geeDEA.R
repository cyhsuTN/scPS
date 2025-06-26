#' DEA Using GEE for Independent Two-Group Comparison
#' @description A function to perform differential expression analysis using GEE
#' for independent two-group comparison.
#' @importFrom stats pt qt

#' @param counts.test Count matrix of tested genes.
#' @param cell.info Cell information (patients' id and x1 are required).
#' @param hf A curve showing the relationship between gene means and gene deviations.
#' The default = NULL, assuming there is no relationship between
#' gene means and gene deviations, and gene variances of group 1 are the same as
#' those of group 2.
#' @param BCC Bias covariance correction at small sample sizes. = 0, no correction;
#' = 1 (Default), t distribution with a degree of freedom = (total size) - 2; = 2,
#' t distribution and MD's bias covariance correction.

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
#' @examples # ## Select 2000 genes that will be tested
#' @examples # idx.selected <- which(aa$nonZeroPs1>0.1 & aa$nonZeroPs2>0.1)
#' @examples # logFC <- log(aa$means2[idx.selected]/aa$means1[idx.selected])
#' @examples # idx.cg <- head(idx.selected[order(abs(logFC), decreasing = T)],2000)
#' @examples # aaP <- aa[idx.cg,]
#' @examples # curve.view <- fitSdmean(preParas=aaP, sp.method = "ns")
#' @examples # hf <- curve.view$hf.sigma
#'
#' @examples # ## Test the selected genes
#' @examples # counts.test <- counts.small[idx.cg,]
#' @examples # geeout <- geeDEA(counts.test, cell.info=cell.info.small, hf=hf, BCC=1)
#' @examples # head(geeout)


#' @export
geeDEA <- function(counts.test, cell.info, hf=NULL, BCC=1) {

  if(ncol(counts.test)!=nrow(cell.info)) {
    stop("ncol(counts.test)!=nrow(cell.info)")
  }

  if(is.null(hf)) hf <- function(x) rep(1, length(x))

  res.muti <- data.frame(t(apply(counts.test, 1, function(x) {
    cell.info$y <- x
    res <- gee.output(cell.info, hf=hf, BCC=BCC)
    c(logFC=res$coefficients[2,1], unadj.p=res$coefficients[2,4], icc=res$icc,
      est.mu=res$est.mus,
      raw.mu=res$raw.mus)
  })))

  return(res.muti)
}



gee.onestep.est <- function(dat, beta0, beta1, rho, hf) {

  id <- unique(dat$id)

  mus <- c(exp(beta0), exp(beta0 + beta1))
  sigmas <- hf(mus)

  DVD.list <- lapply(id, function(x) {

    ni <- sum(dat$id==x)
    tx <- dat$x1[dat$id==x][1]
    yi <- dat$y[dat$id==x]

    mu <- mus[tx+1]
    sigma <- sigmas[tx+1]

    csum.inv.V <- sigma^-2 * 1/(1-rho) * (1 - ni*(rho/(1-rho+ni*rho)) )
    v1 <- mu    * csum.inv.V
    v2 <- mu*tx * csum.inv.V

    Ui <- rbind(v1 * sum(yi - mu), v2 * sum(yi - mu))
    Wi <- ni * mu * rbind( c(v1, v1 * tx),
                           c(v2, v2 * tx))

    list(Wi=Wi, Ui=Ui)

  })

  list_W <- lapply(DVD.list, function(x) {
    x$Wi
  })

  list_U <- lapply(DVD.list, function(x) {
    x$Ui
  })

  W_sum <- Reduce("+", list_W)
  inv.W <- solve(W_sum)

  U_sum <- Reduce("+", list_U)

  as.numeric(c(beta0, beta1) + inv.W %*% U_sum)

}




gee.icc.est <- function(dat, beta0, beta1, hf) {

  mus <- c(exp(beta0), exp(beta0 + beta1))
  sigmas <- hf(mus)

  dat$resi <- ifelse(dat$x1==0, (dat$y-mus[1])/sigmas[1], (dat$y-mus[2])/sigmas[2])

  np <- length(unique(dat$x1))
  e_sum <- sum(tapply(dat$resi, dat$id, function(e) sum(e)^2 - sum(e^2)))
  n_sum <- sum(tapply(dat$resi, dat$id, function(e) length(e)*(length(e)-1))) - 2 * np
  e2_sum <- sum(tapply(dat$resi, dat$id, function(e) sum(e^2)))

  scale.para <- e2_sum/(nrow(dat) - np)  ## Here is the inverse of the scale parameter defined in Liang and Zeger (1986)

  icc <- (e_sum/n_sum)/scale.para

  c(icc=icc, scale.para=scale.para)

}


gee.cov.est <- function(dat, beta0, beta1, rho, hf) {

  id <- unique(dat$id)

  mus <- c(exp(beta0), exp(beta0 + beta1))
  sigmas <- hf(mus)

  DVD.list <- lapply(id, function(x) {
    #x <- id[1]
    ni <- sum(dat$id==x)
    tx <- dat$x1[dat$id==x][1]
    yi <- dat$y[dat$id==x]

    mu <- mus[tx+1]
    sigma <- sigmas[tx+1]

    csum.inv.V <- sigma^-2 * 1/(1-rho) * (1 - ni*(rho/(1-rho+ni*rho)) )
    v1 <- mu    * csum.inv.V
    v2 <- mu*tx * csum.inv.V

    Ui <- rbind(v1 * sum(yi - mu), v2 * sum(yi - mu))
    Wi <- ni * mu * rbind( c(v1, v1 * tx),
                           c(v2, v2 * tx))

    list(Wi=Wi, UiUi=Ui %*% t(Ui))

  })

  list_W <- lapply(DVD.list, function(x) {
    x$Wi
  })

  list_U <- lapply(DVD.list, function(x) {
    x$UiUi
  })

  W_sum <- Reduce("+", list_W)
  inv.W <- solve(W_sum)

  U_sum <- Reduce("+", list_U)

  inv.W %*% U_sum %*% inv.W

}


gee.cov.est.MD <- function(dat, beta0, beta1, rho, hf) {

  id <- unique(dat$id)

  mus <- c(exp(beta0), exp(beta0 + beta1))
  sigmas <- hf(mus)

  DVD.list <- lapply(id, function(x) {

    ni <- sum(dat$id==x)
    tx <- dat$x1[dat$id==x][1]
    yi <- dat$y[dat$id==x]

    mu <- mus[tx+1]
    sigma <- sigmas[tx+1]

    csum.inv.V <- sigma^-2 * 1/(1-rho) * (1 - ni*(rho/(1-rho+ni*rho)) )
    v1 <- mu    * csum.inv.V
    v2 <- mu*tx * csum.inv.V

    Wi <- ni * mu * rbind( c(v1, v1 * tx),
                           c(v2, v2 * tx))

    list(Wi=Wi)

  })

  list_W <- lapply(DVD.list, function(x) {
    x$Wi
  })
  W_sum <- Reduce("+", list_W)
  inv.W <- solve(W_sum)


  DVD.list2 <- lapply(id, function(x) {

    ni <- sum(dat$id==x)
    tx <- dat$x1[dat$id==x][1]
    yi <- dat$y[dat$id==x]

    mu <- mus[tx+1]
    sigma <- sigmas[tx+1]

    csum.inv.V <- sigma^-2 * 1/(1-rho) * (1 - ni*(rho/(1-rho+ni*rho)) )
    v1 <- mu    * csum.inv.V
    v2 <- mu*tx * csum.inv.V

    hii <- sum(c(t(c(mu, mu*tx)) %*% inv.W) * c(v1, v2))
    inv.1.minus.H.e <- (yi - mu) + hii/(1-ni*hii) * sum(yi - mu)

    Ui <- rbind(v1 * sum(inv.1.minus.H.e), v2 * sum(inv.1.minus.H.e))

    list(UiUi=Ui %*% t(Ui))

  })

  list_U <- lapply(DVD.list2, function(x) {
    x$UiUi
  })

  U_sum <- Reduce("+", list_U)

  inv.W %*% U_sum %*% inv.W

}



gee.indep.est <- function(dat, beta0, beta1, hf=function(x) x) {

  id <- unique(dat$id)

  mus <- c(exp(beta0), exp(beta0 + beta1))
  sigmas <- hf(mus)

  estfun <- function(betas) {
    U.list <- lapply(id, function(x) {

      ni <- sum(dat$id==x)
      tx <- dat$x1[dat$id==x][1]
      yi <- dat$y[dat$id==x]

      mu <- mus[tx+1]
      sigma <- sigmas[tx+1]

      v1 <- mu    * sigma^-2
      v2 <- mu*tx * sigma^-2

      rbind(v1 * sum(yi - mu), v2 * sum(yi - mu))

    })

    U_sum <- Reduce("+", U.list)

    as.numeric(U_sum)
  }
  beta00 <- diff(c(0,log(tapply(dat$y, dat$x1, mean))))
  root.solve <- nleqslv(x=beta00, fn=estfun, control=list(ftol=1e-6, maxit=10))
  list(root=root.solve$x, fvec=root.solve$fvec)
}



gee.output <- function(dat, hf, BCC) {

  mus <- mus.pre <- tapply(dat$y, dat$x1, mean)

  if(all(mus==0)) {
    betas.one <- rep(NA,2)

    cov.betas <- NA

    robust.se <- rep(NA,2)
    z.value <- rep(NA,2)
    p.value <- rep(NA,2)

    list(coefficients=data.frame(est=betas.one, robust.se, z.value, p.value),
         cov.betas=cov.betas,
         icc=NA,
         scale.para=NA,
         est.mus=rep(NA,2),
         raw.mus=as.numeric(mus.pre))
  } else {
    if(any(mus==0)) mus <- mus + 1E-8
    betas.initial <- diff(c(0,log(mus)))

    icc.scale <- gee.icc.est(dat, beta0=betas.initial[1], beta1=betas.initial[2], hf)

    betas.one <- gee.onestep.est(dat, beta0=betas.initial[1], beta1=betas.initial[2], rho=icc.scale["icc"], hf)

    if(BCC==2) {
      cov.betas <- tryCatch({
        cov.betas <- gee.cov.est.MD(dat, beta0=betas.one[1], beta1=betas.one[2], rho=icc.scale["icc"], hf)
      }, error = function(e) {
        return(gee.cov.est.MD(dat, beta0=betas.initial[1], betas.initial[2], rho=icc.scale["icc"], hf) )
      })
      var1 <- diag(cov.betas)
      var1 <- ifelse(var1<0, 1E-16, var1)
      robust.se <- sqrt(var1)
      z.value <- betas.one/robust.se
      p.value <- 2 * pt(abs(z.value), df=length(unique(dat$id))-2, lower.tail = FALSE)
    } else if(BCC==1) {
      cov.betas <- tryCatch({
        cov.betas <- gee.cov.est(dat, beta0=betas.one[1], beta1=betas.one[2], rho=icc.scale["icc"], hf)
      }, error = function(e) {
        return(gee.cov.est(dat, beta0=betas.initial[1], betas.initial[2], rho=icc.scale["icc"], hf) )
      })
      var1 <- diag(cov.betas)
      var1 <- ifelse(var1<0, 1E-16, var1)
      robust.se <- sqrt(var1)
      z.value <- betas.one/robust.se
      p.value <- 2 * pt(abs(z.value), df=length(unique(dat$id))-2, lower.tail = FALSE)
    } else {
      cov.betas <- tryCatch({
        cov.betas <- gee.cov.est(dat, beta0=betas.one[1], beta1=betas.one[2], rho=icc.scale["icc"], hf)
      }, error = function(e) {
        return(gee.cov.est(dat, beta0=betas.initial[1], betas.initial[2], rho=icc.scale["icc"], hf) )
      })
      var1 <- diag(cov.betas)
      var1 <- ifelse(var1<0, 1E-16, var1)
      robust.se <- sqrt(var1)
      z.value <- betas.one/robust.se
      p.value <- 2 * pnorm(abs(z.value), lower.tail = FALSE)
    }


    list(coefficients=data.frame(est=betas.one, robust.se, z.value, p.value),
         cov.betas=cov.betas,
         icc=as.numeric(icc.scale["icc"]),
         scale.para=as.numeric(icc.scale["scale.para"]),
         est.mus=as.numeric( c(exp(betas.one[1]), exp(betas.one[1]+betas.one[2])) ),
         raw.mus=as.numeric(mus.pre))
  }

}
