# Gradient
myg <- function(theta, f, f0, eps, ...){
  g <- NULL
  for(j in 1:length(theta)){
    th <- theta; th[j] <- th[j] - eps[j]
    g[j] <- (f0 - f(th, ...))/eps[j]
  }
  g
}

# Gradient search algorithm (used for starting points in Qest).
# The first "choose_eps" must work, otherwise there is no way to proceed.
# The others may not work due to too large theta step. That's why I use "try".
gs <- function(theta0, f, ..., tol = 1e-4, maxit = 100){

  A <- list(...)
  # eps <- choose_eps(A$Q, theta0, A$y, A$data, eps0 = rep(1e-6, length(theta0)), obj = 0.01)
  eps <- rep(.01, length(theta0))
  f0 <- f(theta0, ...)
  g0 <- myg(theta0, f, f0, eps, ...)

  alpha <- 0.001
  conv <- FALSE
  for(i in 1:maxit){

    cond <- FALSE
    while(!cond){
      step <- alpha*g0

      if(max(abs(step)) < tol){conv <- TRUE; break}
      theta1 <- theta0 - step

      f1 <- f(theta1, ...)
      cond <- (!is.na(f1) && is.finite(f1) && f1 < f0)

      if(cond){
        if(i < 5 | round(i/10) == i/10){
          eps_new <- try(choose_eps(A$Q, theta1, A$y, A$data, eps0 = sign(eps)*pmax(abs(eps)/20, 1e-6), obj = 0.01), silent = TRUE)
          cond <- !(inherits(eps_new, "try-error"))
        }
        if(cond){
          # eps <- eps_new
          g1 <- myg(theta1, f, f1, eps, ...)
          cond <- (cond & all(g1 != 0)) # (g = 0) means that the distribution is degenerated
        }
      }

      alpha <- alpha*0.5
    }
    if(conv){break}
    f0 <- f1
    g0 <- g1
    theta0 <- theta1
    alpha <- alpha*5
  }

  # eps <- choose_eps(A$Q, theta0, A$y, A$data, eps0 = sign(eps)*pmax(abs(eps)/30, 1e-6), obj = 0.01)
  list(theta = theta0, f = f0, n.it = i, conv = conv, g = g0, eps = eps, alpha = alpha)
}

# Choose a "meaningful" epsilon for numerical derivatives.
# The starting eps0 must be a vector, can be negative.
choose_eps <- function(Q, theta, y, data, eps0, obj = 0.01){
  tau <- invQ(Q, theta, y, data, n.it = 13)

  npar <- length(theta)
  EPS <- NULL
  for(j in 1:npar){
    cond <- count <- FALSE
    eps <- eps0[j]
    while(!cond){
      if(abs(eps) > 1000){
        if(!count){eps <- -eps0[j]}
        else{stop(paste(
        "\n
        The parameter in position", j, "does not appear to be identified.
        Or, most likely, the starting points are too bad."))}
      }

      eps <- eps*2
      thetastar <- theta
      thetastar[j] <- theta[j] - eps
      taustar <- invQ(Q, thetastar, y, data, n.it = 13)

      if(any(is.na(taustar))){stop("The provided starting points are too bad. The algorithm cannot move from here.")}
      delta_tau <- abs(tau - taustar)
      cond <- if(all(delta_tau == 0)) FALSE else (mean(delta_tau[delta_tau != 0]) > obj)
    }
    EPS[j] <- eps
  }
  EPS
}

# Derivatives of Q(tau | theta) w.r.t. theta
# Each element is a matrix n*ntau
derQtheta <- function(theta, eps, Q, Q1, data, tau, ind){

  if(missing(ind)){ind <- 1:length(theta)}
  Qtheta <- list()
  tau <- tau$TAU
  for(j in ind){
    e <- eps[j]
    thetaL <- theta
    thetaL[j] <- thetaL[j] - e
    der <- (Q1 - callQ(Q, thetaL, tau, data))
    Qtheta <- c(Qtheta, list(der))
  }
  Qtheta # OBS: eps is not divided by here!
}

# Second derivatives of Q(tau | theta) w.r.t. theta
# Each element is a matrix n*ntau
der2Qtheta <- function(theta, eps, Q, Qtheta1, data, tau){

  Q2theta <- list()
  npar <- length(theta)

  for(j in 1:npar){
    ind <- j:npar
    e <- eps[j]
    thetaL <- theta
    thetaL[j] <- thetaL[j] - e
    Q1 <- callQ(Q, thetaL, tau$TAU, data)
    der <- dlist(Qtheta1[ind], derQtheta(thetaL, eps, Q, Q1, data, tau, ind))
    Q2theta <- c(Q2theta, der)
  }
  Q2theta # OBS: eps is not divided by here!
}


# Given a list of matrices n*ntau, calculate their integral functions, by row (another list of matrices n*ntau).
# If ifun = FALSE, only return the integral between 0 and 1 (used for AA(1) only in Qest.ee.u)
# NOTE: A*tau$dtau_long is identical to (t(t(A)*dtau)), but much faster. The same idea is used a lot in Qest.ee.
intA <- function(A, tau, ifun = TRUE){

  intA.internal <- function(A, tau){
    ntau <- tau$ntau
    n <- nrow(A)
    A <- cbind(A[,1], A[,2:ntau] + A[,1:(ntau - 1)])
    if(ifun){out <- rowCumsums(A*tau$dtau_long)}
    else{out <- .rowSums(A*tau$dtau_long, n, ntau)}
    out
  }
  applyfun <- (if(ifun) lapply else sapply)
  applyfun(A, intA.internal, tau = tau)
}


# Given a list "A" in which each element is a matrix of size n*ntau, extrapolates A(p).
# The function assumes that A(0) = 0, which is true when A is an integral, but is not true
# if A corresponds to Q_theta(p). For this reason, when I call findAp(Fy, tau, Qtheta) to compute Qtheta(Fy),
# I preliminary replace v = pmax(v,2), and Fy = pmax(Fy, tau[2]), avoiding any "outL". To make this accurate,
# tau[2] is VERY CLOSE to tau[1]. This is achieved by tau[1] = omicron and tau[2] = omicron*1.1.
findAp <- function(p, tau, A){

  tau1v <- p$tau1v
  dtau1v <- p$dtau1v
  indL <- p$indL
  indR <- p$indR
  outL <- p$outL
  p <- p$Fy

  n <- length(p)

  delta <- p - tau1v
  delta[outL] <- 0 # fix p = tauL

  findAp.internal <- function(A, delta, dtau1v, indL, indR){
    A <- cbind(0, A, 100) # 100 is just a dummy value, it gets multiplied by zero when p = tauR
    AL <- A[indL]; AR <- A[indR]
    dA <- AR - AL
    delta*(dA/dtau1v) + AL
  }

  # A matrix of size n*npar
  sapply(A, findAp.internal, delta = delta, dtau1v = dtau1v, indL = indL, indR = indR)
}


# Given a matrix "Q1" of size n*ntau, that contains vales of Q(p) on a grid, extrapolates F(y).
# This is faster than using invQ bisection, provided that Q1 must be computed anyway (as in Qest.ee etc).
# For how the function is designed, a fix is needed for v = 1 that corresponds to F(y) = tauL.
# Quantile crossing will result in a disaster.

# If "F" is provided by the user, this function computes anyway the quantities used in findAp,
# but returns the value of Fy obtained by applying F(y).

findp <- function(y, tau, Q1){

  n <- length(y)
  ntau1 <- tau$ntau1

  Q1 <- cbind(-Inf, Q1, Inf)
  delta <- (y - Q1)
  v <- .rowSums(delta >= 0, n, ntau1) # between 1 and (ntau1 - 1)

  outL <- which(v == 1)
  outR <- which(v == (ntau1 - 1))
  indL <- cbind(1:n, v)
  indR <- cbind(1:n, v + 1)

  tau1v <- tau$tau1[v]
  dtau1v <- tau$dtau1[v]

  dQ1 <- Q1[indR] - Q1[indL]
  Fy <- delta[indL]*pmin((dtau1v/dQ1), 1e+200) + tau1v
  # if(any(is.infinite(1/dQ1))) Fy[is.infinite(1/dQ1)] <- tau1v[is.infinite(1/dQ1)]
  Fy[outL] <- tau$tauL # fix p = tauL

  list(Fy = Fy, delta = delta[,2:(ntau1 - 1)], tau1v = tau1v, dtau1v = dtau1v,
       indL = indL, indR = indR, outL = outL, outR = outR)
}

# Given Q, computes F. This is faster than using findp, if you can avoid computing the whole Q(p)
# on a n*ntau grid. This is the case in choose_eps, which is where invQ is used.
# Quantile crossing will result in a disaster.
invQ <- function(Q, theta, y, data, n.it = 17){

  n <- length(y)
  tau <- rep.int(0.5,n)
  for(i in 2:n.it){
    omega <- sign(y - callQ(Q, theta, tau, data))
    tau <- tau + 1/(2^i)*omega
  }
  tau
}


# Control parameters
Qest.control <- function(tol = 1e-8, maxit, safeit, alpha0, display = FALSE, restart = FALSE){
  if(missing(maxit)){maxit <- NULL}
  if(missing(alpha0)){alpha0 <- NULL}
  if(missing(safeit)){safeit <- NULL}
  list(tol = tol, maxit = maxit, alpha0 = alpha0, safeit = safeit, display = display, restart = restart)
}


