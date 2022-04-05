# Qest covariance matrix
Qest.covar <- function(fit, eps, w){

  g.i <- t(t(fit$ee.i)/eps)
  J <- fit$jacobian*nrow(g.i) # J was divided by n in Qest.newton

  Omega <- chol2inv(chol(t(w*g.i)%*%g.i))
  V <- t(J)%*%Omega%*%J
  chol2inv(chol(V))
}

# Compute omega = I(y <= Q(tau | X)).
# If the data are censored or truncated, returns \hat\omega.
# Only used in Loss
omega <- function(d, tau, type, Fy, Fz){

  tau <- tau$TAU

  Sy <- Sz <- 1
  if(type != "u"){Sy <- 1 - Fy$Fy}
  if(type == "ct"){Sz <- 1 - Fz$Fy}

  deltay <- Fy$delta
  omegay <- (deltay <= 0)
  if(type == "u"){omega <- omegay}
  else{omega <- omegay + (1 - d)*omegay*(tau - 1)/Sy}
  if(type == "ct"){
    deltaz <- Fz$delta
    omegaz <- (deltaz <= 0)
    omega <- omega - omegaz - omegaz*(tau - 1)/Sz + tau
  }

  list(deltay = deltay, omega = omega)
}

# Loss function. If the data are neither censored nor truncated, there is an "analytical"
# expression that uses A and AA. However, its evaluation is not going to be particularly fast,
# and I just use numerical integration in all cases.
Loss <- function(w, d, tau, type, Fy, Fz){

  q <- omega(d, tau, type, Fy, Fz)
  out <- q$deltay*(tau$TAU - q$omega)
  out <- out %*% (tau$dtau*tau$wtau)
  2*sum(w*out) # the dtau was divided by 2
}
