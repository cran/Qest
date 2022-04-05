# QCoxph covariance matrix

Qcox.covar <- function(theta, z, y, d, X, w, knots, tau, type){

  # Gradient and Jacobian
  ee <- (if(type == "ct") QCox.ee.ct else QCox.ee.c)
  G <- ee(theta, 0, z, y, d, X, w, knots, tau, J = TRUE)
  g.i <- G$g.i
  J <- G$J

  Omega <- chol2inv(chol(t(w*g.i)%*%g.i))
  V <- t(J)%*%Omega%*%J
  chol2inv(chol(V))
}

# Loss of Cox model.
# If the data are not censored and not truncated, there is an analytical expression:
# A <- .rowSums(t(t(BBy$dWmins)*BBy$c1) + BBy$dRmins, n, k + 1)
# AA1 <- .rowSums(t(t(BBy$diWknots_dWknots1)*BBy$c1) + BBy$diRknots_dRknots1, n, k + 1)
# L <- sum(w*(y*(tau$Wfun(Fy) - tau$iWfun(1)) + AA1 - A))
# However, I assume that cox model will (almost) always be used with censored or truncated data
# (which is why I did not bother to implement QCox.ee.u), and I keep the code simple by
# directly computing the numeric integral that is needed with censored or truncated data.
coxLoss <- function(theta, z, y, d, X, w, knots, tau, type, Fy, Fz){

  Ltau <- function(opt, tau){opt$tau <- tau; opt}
  n <- nrow(X)
  q <- ncol(X)
  k <- length(knots) - 2
  beta <- (theta[1:q]) # Cox regression coefficients
  egamma <- exp(theta[(q + 1):length(theta)]) # coefficients of H0
  HR1 <- exp(-c(X%*%beta))

  # Quantile function evaluated on a grid

  taugrid <- c(1e-6, (1:99)/100, 1 - 1e-6); dtau <- 0.01; ntau <- length(taugrid)
  TAU <- t(matrix(taugrid, ntau, n))
  log1_tau <- log(1 - taugrid)
  u <- -HR1*rep(log1_tau, each = n) # vector of length n*ntau

  knotsinv <- c(plfcox(knots, knots = knots)%*%egamma) # knots of the inverse of b
  knotsinv_long <- unique(c(0, knotsinv, -max(HR1)*log1_tau[ntau])) # plfcox can project outside range, approx cannot.
  binv <- plfcox(knotsinv_long, knotsinv)

  Qcoxtemp <- c(binv%*%(1/egamma))
  Qcox <- approx(knotsinv_long, Qcoxtemp, xout = u, method = "linear")$y + knots[1]
  qq <- matrix(Qcox, n, ntau) # A matrix n*ntau.

  # Pseudo-Loss for censored or truncated data. With "u" data, omega = omegay.

  delta <- y - qq
  omegay <- (delta <= 0)
  omega <- omegay + (1 - d)*omegay*(TAU - 1)/(1 - pmin(Fy, 1 - 1e-10))
  if(type == "ct"){
    deltaz <- z - qq
    omegaz <- (deltaz <= 0)
    omega <- omega - omegaz - omegaz*(TAU - 1)/(1 - pmin(Fz, 1 - 1e-10)) + TAU
  }

  wtau <- do.call(tau$wfun, Ltau(tau$opt, taugrid))
  out <- delta*(TAU - omega)
  out <- out %*% (dtau*wtau)
  sum(w*out)
}
