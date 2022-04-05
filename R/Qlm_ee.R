# Qlm estimating equation

Qlm.ee.u <- function(theta, X, w, bfun, EE, J = FALSE){

  sigma <- exp(theta[1]); beta <- theta[-1]
  p.i <- EE$p.i; delta <- EE$delta; deltaW <- EE$deltaW; deltaZW <- EE$deltaZW

  # gradient

  g.sigma.i <- -sigma*deltaZW
  g.sigma <- sum(w*g.sigma.i)
  g.beta.i <- -(X*deltaW)
  g.beta <- c(w%*%g.beta.i)
  g <- c(g.sigma, g.beta)
  g.i <- cbind(g.sigma.i, g.beta.i)

  H <- NULL
  if(J){
    # hessian
    f <- dnorm(delta,0,sigma)
    wf <- w*f
    zi <- bfun$z(p.i)
    wi <- do.call(bfun$w, Ltau(bfun$opt, p.i))
    zwi <- zi*wi

    a <- wf*X*wi
    Hbeta <- t(a)%*%X
    Hsigma <- sigma^2*sum(wf*zwi*zi) + g[1]
    Hbetasigma <- sigma*rep(1,nrow(X))%*%(a*zi)
    H <- cbind(c(Hsigma, Hbetasigma), rbind(Hbetasigma, Hbeta))
  }

  # output

  list(g.i = g.i, g = g, H = H)
}
