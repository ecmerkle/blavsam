gibbs_step2 <- function(fit, step1samps, blpt, psinblk, psiblkse) {
  ## fit: sam object with do.fit = FALSE
  ## blpt: blavaan partable (with prior column)
  stopifnot(lavInspect(fit, "ngroups") == 1)

  dat <- lavInspect(fit, "data")
  step1mod <- lavaan(subset(parTable(fit), step == 1), data = dat)
  
  niter <- nrow(step1samps)

  ## info from full model
  mats <- lavInspect(fit, "est")
  frmats <- lavInspect(fit, "free")
  p <- nrow(mats$lambda)
  m <- ncol(mats$lambda)
  std.lv <- lavInspect(fit, "options")$std.lv

  ## indexing of free parms in alpha/beta + priors
  if (!("beta" %in% names(frmats))) frmats$beta <- array(0, dim=0)
  lvmat <- t(with(frmats, cbind(alpha, beta)))
  fridx2 <- NULL
  alphidx <- NULL
  betidx <- NULL
  if (any(lvmat > 0)) {
    fridx2 <- which(lvmat > 0)
    alphidx <- which(which(lvmat > 0, arr.ind = TRUE)[,'row'] == 1)
    betidx <- which(t(mats$beta) != 0)
  }
  fixalphidx <- which(frmats$alpha == 0)

  parnums <- t(as.numeric(lvmat))[t(as.numeric(lvmat)) > 0]
  pris <- blavaan:::dist2r(blpt$prior[match(parnums, blpt$free)], target = "stan")
  primn <- as.numeric(sapply(pris, function(x) x[2]))
  priprec <- as.numeric(sapply(pris, function(x) {
    pow <- -2
    if (length(x) == 4) {
      if (grepl("prec", x[4])) pow <- 1
      if (grepl("var", x[4])) pow <- -1
    }
    x[3]^pow
  }))
  priprec <- diag(priprec)
  
  ## indexing of free parm in psi (if any)
  psiidx <- which(diag(frmats$psi) > 0)
  psi1idx <- which(diag(mats$psi) == 1)
  corrfac <- rep(1, m)

  parnums <- diag(frmats$psi)[diag(frmats$psi) > 0]
  pris <- blavaan:::dist2r(blpt$prior[match(parnums, blpt$free)], target = "stan")
  psi_prior_shape <- as.numeric(sapply(pris, function(x) x[2]))
  psi_prior_shape <- diag(psi_prior_shape)
  psi_prior_rate <- as.numeric(sapply(pris, function(x) x[3]))
  psi_prior_rate <- diag(psi_prior_rate)
  
  ## initialize alpha, Beta, psi using final step 1 draw
  lavmod <- blavaan:::fill_params(colMeans(step1samps), step1mod@Model, parTable(step1mod))
  mats <- lavmod@GLIST ## TODO Ng > 1
  Lambda <- mats$Lambda
  Theta <- mats$Theta
  Nu <- mats$Nu
  etainit <- sample_eta(fit, dat, Nu, matrix(0, m), Lambda, matrix(0, m, m), Theta, diag(m), corrfac)
  locinit <- sample_loc(matrix(0, m), matrix(0, m, m), diag(m), etainit, fridx2, alphidx, betidx, primn, priprec)
  Alpha <- locinit$Alpha
  Beta <- locinit$Beta
  scinit <- sample_scale(eta, Alpha, Beta, psinblk, psiblkse, psi_prior_shape, psi_prior_rate, std.lv, psi1idx)
  Psi <- scinit$Psi
  corrfac <- scinit$corrfac

  samps <- vector("list", niter)
  for (i in 1:niter) {
    ## use blavaan:::fill_params to fill in step 1 model with step1samps;
    ## get filled matrices from GLIST
    lavmod <- blavaan:::fill_params(step1samps[i,], step1mod@Model, parTable(step1mod))
    mats <- lavmod@GLIST ## TODO Ng > 1
    Lambda <- mats$Lambda
    Theta <- mats$Theta
    Nu <- mats$Nu

    ## sample factor scores
    eta <- sample_eta(fit, dat, Nu, Alpha, Lambda, Beta, Theta, Psi, corrfac)

    ## sample alpha, beta
    locparms <- sample_loc(Alpha, Beta, Psi, eta, alphidx, betidx, primn, priprec)
    Alpha <- locparms$Alpha
    Beta <- locparms$Beta

    ## sample psi
    scparms <- sample_scale(eta, Alpha, Beta, psinblk, psiblkse, psi_prior_shape, psi_prior_rate, std.lv, psi1idx)
    Psi <- scparms$Psi
    corrfac <- scparms$corrfac

    ## save samples TODO save covs separately from vars
    samps[[i]] <- c(locparms$alpha, locparms$beta, diag(Psi), as.numeric(t(eta)))
  }

  samps
}


sample_eta <- function(fit, dat, Nu, Alpha, Lambda, Beta, Theta, Psi, corrfac) {
  dummy.ov.x.idx <- unlist(fit@Model@ov.x.dummy.ov.idx)
  dummy.lv.x.idx <- unlist(fit@Model@ov.x.dummy.lv.idx)
  dummy.ov.idx <- unlist(c(fit@Model@ov.x.dummy.ov.idx, fit@Model@ov.y.dummy.ov.idx))
  dummy.lv.idx <- unlist(c(fit@Model@ov.x.dummy.lv.idx, fit@Model@ov.y.dummy.lv.idx))
  m <- ncol(Beta)

  IBinv <- solve(diag(m) - Beta)
  IBinv[dummy.lv.x.idx, ] <- 0; IBinv[, dummy.lv.x.idx] <- 0; diag(IBinv)[dummy.lv.x.idx] <- 1
  Psi0 <- IBinv %*% Psi %*% t(IBinv)
  LamtThetinv <- t(Lambda) %*% solve(Theta)
  D <- solve(LamtThetinv %*% Lambda + solve(Psi0))
  d <- apply(dat, 1, function(y) LamtThetinv %*% (y - Nu) + solve(Psi0) %*% IBinv %*% (Alpha + Beta[, dummy.lv.x.idx, drop = FALSE] %*% y[dummy.ov.x.idx]))
  if (m == 1) d <- t(d)
  fs <- t(D %*% d)

  for (j in 1:ntot) {
    eta[j,] <- rmnorm(1, fs[j,], D)
  }
  for (j in 1:m) {
    eta[,j] <- eta[,j] / corrfac[j]
  }
  if (length(dummy.lv.idx) > 0) {
    eta[, dummy.lv.idx] <- dat[, dummy.ov.idx]
  }

  eta
}


sample_loc <- function(Alpha, Beta, Psi, eta, fridx2, alphidx, betidx, primn, priprec) {
  m <- ncol(Psi)
  ntot <- nrow(eta)
  F_i <- lapply(1:ntot, function(ii) diag(m) %x% matrix(c(1, eta[ii,]), 1))
  F_i <- lapply(F_i, function(F) F[, fridx2])

  Vinv <- solve(Psi)

  FVF <- lapply(F_i, function(F) t(F) %*% Vinv %*% F)
  ## TODO in first part of FVF, sum equality-constrained rows and columns, only save the first
  FVFsum <- Reduce("+", FVF) + priprec
  D <- solve(FVFsum)

  fixsubl <- rep(0, m)
  FVz <- lapply(1:ntot, function(i) {
    fixsubl[fixalphidx] <- Alpha[fixalphidx]
    t(F_i[[i]]) %*% Vinv %*% (eta[i, ] - fixsubl)
      })
  d <- Reduce("+", FVz) + priprec %*% primn

  ## TODO sum entries of equality-constrained parameters, only save the first one
    
  ## samples free parameters
  lvpars <- rmnorm(1, D %*% d, D)

  Alpha[frmats$alpha > 0] <- lvpars[alphidx]
  bet <- lvpars
  if (length(alphidx) > 0) {
    bet <- bet[-alphidx]
  }
  tbeta <- t(Beta)    
  tbeta[betidx] <- bet
  Beta <- t(tbeta)

  list(Alpha = Alpha, Beta = Beta, alpha = lvpars[alpidx], beta = bet)
}


sample_scale <- function(eta, Alpha, Beta, psinblk, psiblkse, psi_prior_shape, psi_prior_rate, std.lv, psi1idx) {
  ntot <- nrow(eta)
  m <- ncol(Beta)

  residcp <- tcrossprod(sapply(1:ntot, function(k) (eta[k,] - Alpha - Beta %*% eta[k,])^2))

  ## TODO handle psiorder, psirevord as in stanmarg?
  Psiblk <- matrix(0, m, m)
  for (k in 1:psinblk) {
    srow <- psiblkse[k,1]
    erow <- psiblkse[k,2]

    if (erow > srow) {
      Psiblk[srow:erow, srow:erow] <- MCMCpack::riwish(1, ntot + psi_prior_shape[srow, srow], residcp[srow:erow, srow:erow] + psi_prior_rate(srow:erow, srow:erow))
    } else {
      Psiblk[srow, srow] <- 1/rgamma(1, .5 * ntot + psi_prior_shape[k, k], rate = .5 * residcp[k,k] + psi_prior_rate[k, k])
    }
  }

  if (std.lv) {
    for (k in 1:m) {
      if (k %in% psi1idx) {
        corrfac[k] <- sqrt(1/rgamma(1, .5 * ntot + psi_prior_shape[k, k], .5 * residcp[k,k] + psi_prior_rate[k, k]))
      }
    }
  }
      
  list(Psi = Psiblk, corrfac = corrfac)
}
