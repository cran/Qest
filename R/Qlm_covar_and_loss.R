# Qlm covariance matrix
Qlm.covar <- function(g.i, w, H){
  omega <- crossprod(g.i*sqrt(w))
  Hinv <- chol2inv(chol(H))
  covar <- Hinv%*%omega%*%Hinv
  covar
}

# Qlm loss function
qlmLoss <- function(theta, y, X, w, bfun){

  # building blocks

  sigma <- exp(theta[1]); beta <- theta[-1]
  delta <- c(y - X%*%beta)
  p.i <- pmax(pmin(pnorm(delta,0,sigma), 1-1e-6), 1e-6)
  deltaW <- bfun$W(p.i) - bfun$Wbar
  deltaZW <- bfun$ZW(p.i) - bfun$ZWbar

  Loss.i <- delta*(deltaW) - sigma*(deltaZW)
  Loss <- sum(w*Loss.i)
  return(list(Loss = Loss, p.i = p.i, delta = delta, deltaW = deltaW, deltaZW = deltaZW))
}
