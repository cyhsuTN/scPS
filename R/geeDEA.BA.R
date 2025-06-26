#' DEA Using GEE for Paired-Group Comparison
#' @description A function to perform differential expression analysis using GEE
#' for paired-group comparison.
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
#' @examples # load(file = "DataForDemo/GSE120575n.rda")
#' @examples # counts <- GSE120575n$counts
#' @examples # cell.info <- GSE120575n$cell.info
#' @examples # cells.interesting <- which(cell.info$cellcluster=="NK")
#' @examples # counts.small <- counts[,cells.interesting]
#' @examples # cell.info.small <- data.frame(id=cell.info[cells.interesting, "ptID"],
#' @examples #        x1=ifelse(cell.info[cells.interesting, "TX"]=="Pre", 0, 1) )
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
#' @examples # geeout <- geeDEA.BA(counts.test, cell.info=cell.info.small, hf=hf, BCC=1)
#' @examples # head(geeout)


#' @export
geeDEA.BA <- function(counts.test, cell.info, hf=NULL, BCC=1) {

  if(ncol(counts.test)!=nrow(cell.info)) {
    stop("ncol(counts.test)!=nrow(cell.info)")
  }

  if(is.null(hf)) hf <- function(x) rep(1, length(x))

  tab <- table(cell.info$id, cell.info$x1)
  if(any(apply(tab, 1, function(x) any(x==0)))) {
    rm.id <- rownames(tab)[apply(tab, 1, function(x) any(x==0))]
    counts.test <- counts.test[,!(cell.info$id %in% rm.id)]
    cell.info <- cell.info[!(cell.info$id %in% rm.id),]
    print(paste0("Warning: id ", paste0(rm.id, collapse = ", "),
                 " are excluded from calculation because empty in one of paired groups"))
  }


  res.muti <- data.frame(t(apply(counts.test, 1, function(x) {
    cell.info$y <- x
    res <- gee.output.BA(cell.info, hf=hf, BCC=BCC)

    c(logFC=res$coefficients[2,1], unadj.p=res$coefficients[2,4], icc=res$icc,
      est.mu=res$est.mus,
      raw.mu=res$raw.mus)
  })))

  return(res.muti)
}



gee.onestep.est.BA <- function(dat, beta0, beta1, rho, hf) {

  id <- unique(dat$id)

  mus <- c(exp(beta0), exp(beta0 + beta1))
  sigmas <- hf(mus)
  mu.by.sigma <- mus/sigmas

  DVD.list <- lapply(id, function(x) {

    ns <- as.numeric(table(dat$x1[dat$id==x]))
    tx <- dat$x1[dat$id==x]
    yi <- dat$y[dat$id==x]

    y.minus.mu.by.sigma <- as.numeric(tapply(yi-mus[tx+1], tx, sum)/sigmas)

    n1 <- ns[1]
    n2 <- ns[2]

    tp1 <- 1/(1-rho)
    tp2 <- rho/(1-rho+(n1+n2)*rho)

    csum.inv.M1 <- tp1 * (1 - n1*tp2 )
    csum.inv.M2 <- tp1 * (  - n2*tp2 )
    csum.inv.M3 <- tp1 * (  - n1*tp2 )
    csum.inv.M4 <- tp1 * (1 - n2*tp2 )

    v11 <- sum(mu.by.sigma * c(csum.inv.M1, csum.inv.M2))
    v12 <- sum(mu.by.sigma * c(csum.inv.M3, csum.inv.M4))
    v21 <- sum(c(0,mu.by.sigma[2]) * c(csum.inv.M1, csum.inv.M2))
    v22 <- sum(c(0,mu.by.sigma[2]) * c(csum.inv.M3, csum.inv.M4))

    Ui <- rbind(sum(c(v11, v12) * y.minus.mu.by.sigma),
                sum(c(v21, v22) * y.minus.mu.by.sigma))

    Wi <- rbind(c(sum(ns * c(v11, v12) * mu.by.sigma), sum(ns * c(v11, v12) * c(0,mu.by.sigma[2]))),
                c(sum(ns * c(v21, v22) * mu.by.sigma), sum(ns * c(v21, v22) * c(0,mu.by.sigma[2]))) )

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




gee.icc.est.BA <- function(dat, beta0, beta1, hf) {

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


gee.cov.est.BA <- function(dat, beta0, beta1, rho, hf) {

  id <- unique(dat$id)

  mus <- c(exp(beta0), exp(beta0 + beta1))
  sigmas <- hf(mus)
  mu.by.sigma <- mus/sigmas

  DVD.list <- lapply(id, function(x) {
    ns <- as.numeric(table(dat$x1[dat$id==x]))
    tx <- dat$x1[dat$id==x]
    yi <- dat$y[dat$id==x]

    y.minus.mu.by.sigma <- as.numeric(tapply(yi-mus[tx+1], tx, sum)/sigmas)

    n1 <- ns[1]
    n2 <- ns[2]

    tp1 <- 1/(1-rho)
    tp2 <- rho/(1-rho+(n1+n2)*rho)

    csum.inv.M1 <- tp1 * (1 - n1*tp2 )
    csum.inv.M2 <- tp1 * (  - n2*tp2 )
    csum.inv.M3 <- tp1 * (  - n1*tp2 )
    csum.inv.M4 <- tp1 * (1 - n2*tp2 )

    v11 <- sum(mu.by.sigma * c(csum.inv.M1, csum.inv.M2))
    v12 <- sum(mu.by.sigma * c(csum.inv.M3, csum.inv.M4))
    v21 <- sum(c(0,mu.by.sigma[2]) * c(csum.inv.M1, csum.inv.M2))
    v22 <- sum(c(0,mu.by.sigma[2]) * c(csum.inv.M3, csum.inv.M4))

    Ui <- rbind(sum(c(v11, v12) * y.minus.mu.by.sigma),
                sum(c(v21, v22) * y.minus.mu.by.sigma))

    Wi <- rbind(c(sum(ns * c(v11, v12) * mu.by.sigma), sum(ns * c(v11, v12) * c(0,mu.by.sigma[2]))),
                c(sum(ns * c(v21, v22) * mu.by.sigma), sum(ns * c(v21, v22) * c(0,mu.by.sigma[2]))) )

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


gee.cov.est.BA.MD <- function(dat, beta0, beta1, rho, hf) {

  id <- unique(dat$id)

  mus <- c(exp(beta0), exp(beta0 + beta1))
  sigmas <- hf(mus)
  mu.by.sigma <- mus/sigmas

  DVD.list <- lapply(id, function(x) {
    ns <- as.numeric(table(dat$x1[dat$id==x]))
    tx <- dat$x1[dat$id==x]
    yi <- dat$y[dat$id==x]

    n1 <- ns[1]
    n2 <- ns[2]

    tp1 <- 1/(1-rho)
    tp2 <- rho/(1-rho+(n1+n2)*rho)

    csum.inv.M1 <- tp1 * (1 - n1*tp2 )
    csum.inv.M2 <- tp1 * (  - n2*tp2 )
    csum.inv.M3 <- tp1 * (  - n1*tp2 )
    csum.inv.M4 <- tp1 * (1 - n2*tp2 )

    v11 <- sum(mu.by.sigma * c(csum.inv.M1, csum.inv.M2))
    v12 <- sum(mu.by.sigma * c(csum.inv.M3, csum.inv.M4))
    v21 <- sum(c(0,mu.by.sigma[2]) * c(csum.inv.M1, csum.inv.M2))
    v22 <- sum(c(0,mu.by.sigma[2]) * c(csum.inv.M3, csum.inv.M4))

    Wi <- rbind(c(sum(ns * c(v11, v12) * mu.by.sigma), sum(ns * c(v11, v12) * c(0,mu.by.sigma[2]))),
                c(sum(ns * c(v21, v22) * mu.by.sigma), sum(ns * c(v21, v22) * c(0,mu.by.sigma[2]))) )

    list(Wi=Wi)

  })

  list_W <- lapply(DVD.list, function(x) {
    x$Wi
  })
  W_sum <- Reduce("+", list_W)
  inv.W <- solve(W_sum)


  DVD.list2 <- lapply(id, function(x) {
    ns <- as.numeric(table(dat$x1[dat$id==x]))
    tx <- dat$x1[dat$id==x]
    yi <- dat$y[dat$id==x]

    n1 <- ns[1]
    n2 <- ns[2]

    tp1 <- 1/(1-rho)
    tp2 <- rho/(1-rho+(n1+n2)*rho)

    csum.inv.M1 <- tp1 * (1 - n1*tp2 )
    csum.inv.M2 <- tp1 * (  - n2*tp2 )
    csum.inv.M3 <- tp1 * (  - n1*tp2 )
    csum.inv.M4 <- tp1 * (1 - n2*tp2 )

    v11 <- sum(mu.by.sigma * c(csum.inv.M1, csum.inv.M2))
    v12 <- sum(mu.by.sigma * c(csum.inv.M3, csum.inv.M4))
    v21 <- sum(c(0,mu.by.sigma[2]) * c(csum.inv.M1, csum.inv.M2))
    v22 <- sum(c(0,mu.by.sigma[2]) * c(csum.inv.M3, csum.inv.M4))

    DW1 <- t(c(mu.by.sigma[1],0))              %*% inv.W
    DW2 <- t(c(mu.by.sigma[2],mu.by.sigma[2])) %*% inv.W

    a <- sum(c(DW1) * c(v11, v21))
    b <- sum(c(DW1) * c(v12, v22))
    c <- sum(c(DW2) * c(v11, v21))
    d <- sum(c(DW2) * c(v12, v22))

    one.minus.Hii <- diag(1, (n1+n2)) - rbind(cbind(matrix(a, n1, n1), matrix(b,n1,n2)),
                                              cbind(matrix(c, n2, n1), matrix(d,n2,n2)))

    k1 <- c*b*n1 + c*b*a*n1^2/(1-n1*a) + d
    k2 <- (1 + a*n1/(1-n1*a)) * (1 + k1*n2/(1-n2*k1))
    k3 <- a/(1-n1*a) + c*b*k2*n1*(1+a*n1/(1-n1*a))

    y.minus.mu.by.sigma <- as.numeric(tapply(yi-mus[tx+1], tx, sum)/sigmas)

    y.minus.mu.by.sigma <- y.minus.mu.by.sigma + c(n1 * sum(c(k3, b*k2) * y.minus.mu.by.sigma),
                                                   n2 * sum(c(c*k2, k1/(1-n2*k1)) * y.minus.mu.by.sigma))

    Ui <- rbind(sum(c(v11, v12) * y.minus.mu.by.sigma),
                sum(c(v21, v22) * y.minus.mu.by.sigma))

    list(UiUi=Ui %*% t(Ui))

  })

  list_U <- lapply(DVD.list2, function(x) {
    x$UiUi
  })

  U_sum <- Reduce("+", list_U)

  inv.W %*% U_sum %*% inv.W

}





gee.output.BA <- function(dat, hf, BCC) {

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

    icc.scale <- gee.icc.est.BA(dat, beta0=betas.initial[1], beta1=betas.initial[2], hf)

    betas.one <- gee.onestep.est.BA(dat, beta0=betas.initial[1], beta1=betas.initial[2], rho=icc.scale["icc"], hf)

    if(BCC==2) {
      cov.betas <- tryCatch({
        cov.betas <- gee.cov.est.BA.MD(dat, beta0=betas.one[1], beta1=betas.one[2], rho=icc.scale["icc"], hf)
      }, error = function(e) {
        return(gee.cov.est.BA.MD(dat, beta0=betas.initial[1], betas.initial[2], rho=icc.scale["icc"], hf) )
      })
      var1 <- diag(cov.betas)
      var1 <- ifelse(var1<0, 1E-16, var1)
      robust.se <- sqrt(var1)
      z.value <- betas.one/robust.se
      p.value <- 2 * pt(abs(z.value), df=length(unique(dat$id))-2, lower.tail = FALSE)
    } else if(BCC==1) {
      cov.betas <- tryCatch({
        cov.betas <- gee.cov.est.BA(dat, beta0=betas.one[1], beta1=betas.one[2], rho=icc.scale["icc"], hf)
      }, error = function(e) {
        return(gee.cov.est.BA(dat, beta0=betas.initial[1], betas.initial[2], rho=icc.scale["icc"], hf) )
      })
      var1 <- diag(cov.betas)
      var1 <- ifelse(var1<0, 1E-16, var1)
      robust.se <- sqrt(var1)
      z.value <- betas.one/robust.se
      p.value <- 2 * pt(abs(z.value), df=length(unique(dat$id))-2, lower.tail = FALSE)
    } else {
      cov.betas <- tryCatch({
        cov.betas <- gee.cov.est.BA(dat, beta0=betas.one[1], beta1=betas.one[2], rho=icc.scale["icc"], hf)
      }, error = function(e) {
        return(gee.cov.est.BA(dat, beta0=betas.initial[1], betas.initial[2], rho=icc.scale["icc"], hf) )
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

