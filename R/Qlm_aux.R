
# numerical integral
num.fun <- function(dx,fx){
  n <- length(dx) + 1
  fL <- fx[1:(n-1)]
  fR <- fx[2:n]
  cumsum(c(0, 0.5*dx*(fL + fR)))
}

Qlm.bfun <- function(wtau, ...){

  r <- 4999
  p <- pbeta(seq.int(qbeta(1e-6,2,2), qbeta(1 - 1e-6,2,2), length.out = r),2,2)
  dp <- p[-1] - p[-r]

  z <- qnorm(p)
  w <- wtau(p, ...)
  zw <- z*w
  W <- num.fun(dp,w)
  ZW <- num.fun(dp,zw)
  Wbar <- num.fun(dp,W)[r]
  ZWbar <- num.fun(dp,ZW)[r]

  z <- splinefun(p,z, method = "hyman") # faster than qnorm, no trouble at tau = 0 and tau = 1
  list(
    z = z,
    w = wtau,
    W = splinefun(p,W, method = "hyman"),
    ZW = approxfun(p,ZW, rule = 2), # may not be monotone, use approxfun
    Wbar = Wbar, ZWbar = ZWbar,
    opt = list(...)
  )
}

# Transform the scale of y and X in a linear model.
scalevars.Qlm <- function(X,y){

  sX <- apply(X, 2, sd)
  int <- which(sX == 0)
  vars <- which(sX > 0)
  is.int <- (length(int) > 0)
  mX <- colMeans(X)*is.int # zero if the model has no intercept

  X[,vars] <- scale(X[,vars], center = mX[vars], scale = sX[vars])

  my <- min(y); My <- max(y)
  y <- (y - my*is.int)/(My - my)*10 # y goes in (0,10). If no intercept, y is only scaled, not shifted
  Stats <- list(mX = mX, sX = sX, int = int, is.int = is.int, vars = vars, my = my, My = My)

  list(X = X, y = y, Stats = Stats)
}

start.Qlm <- function(x, y, w, start, ok, Stats){
  # if(is.null(start)) {
    # Qy <- quantile(y, probs = c(.05, .95))
    # use <- (y > Qy[1] & y < Qy[2])
    # temp <- glm.fit(x[use, ], y[use], w[use], family=gaussian())
    # temp$dispersion <- sum((w[use] * temp$residuals^2)[w[use] > 0]) / temp$df.residual
    # names(temp$dispersion) <- "log(sigma)"
    # names(temp$coefficients) <- colnames(x)
    # theta <- c(log(temp$dispersion), temp$coefficients)
  p1 <- .7; p2 <- .3
  temp <- rq.fit.br2(x, y, tau = .5)
  temp$dispersion <- (quantile(temp$residuals, p1) - quantile(temp$residuals, p2))/(qnorm(p1) - qnorm(p2))
  names(temp$dispersion) <- "log(sigma)"
  names(temp$coefficients) <- colnames(x)
  theta <- c(log(temp$dispersion), temp$coefficients)
  # }
  # else{
  #   npar <- length(ok) + 1
  #   if(length(start) != npar) stop("Wrong size of 'start'")
  #   # if(any(is.na(start))) stop("NAs are not allowed in 'start'")
  #   start[is.na(start)] <- 0.
  #
  #   K <- (Stats$My - Stats$my) / 10
  #   theta <- double(sum(ok) + 1)
  #   theta[1] <- start[1] - log(K) #ok
  #   theta[-1] <- start[-1][ok]
  #
  #   if(length(Stats$vars) == 0){theta[2] <- (theta[2] - Stats$my) / K}
  #   else{
  #     beta <- theta[-1]
  #     beta[Stats$vars] <- (beta * Stats$sX)[Stats$vars] / K
  #     if(Stats$is.int) beta[Stats$int] <- (beta[Stats$int]  - Stats$my) / K + sum((beta * Stats$mX / Stats$sX)[Stats$vars])
  #     theta[-1] <- beta
  #   }

    nms <- c("log(sigma)", colnames(x))
    names(theta) <- nms
  # }
  theta
}


# Transforms back the vector of coefficients in a linear model
descalecoef.Qlm <- function(theta, Stats){

  beta <- theta[-1]
  log.sigma <- theta[1]

  beta[Stats$int] <- beta[Stats$int] - sum((beta*Stats$mX/Stats$sX)[Stats$vars])
  beta[Stats$vars] <- (beta/Stats$sX)[Stats$vars]

  beta <- beta * (Stats$My - Stats$my) / 10
  beta[Stats$int] <- beta[Stats$int] + Stats$my*Stats$is.int

  log.sigma <- log.sigma + log((Stats$My - Stats$my)/10)
  c(log.sigma, beta)
}



########### NEWTON-RAPHSON ########################################

invJ <- function(J, type){
  if(type == "u") {
    Jinv <- try(chol(J), silent = TRUE)
    err <- (inherits(Jinv, "try-error"))
  }
  else{
    Jinv <- qr(J)
    r <- Jinv$rank
    err <- (r != ncol(J))
  }
  list(Jinv = Jinv, err = err)
}

Qlm.sgs.internal <- function(theta, type, tol, maxit, alpha0, ee, display, n.it, y, X, w, bfun){

  n <- length(y)
  L0 <- qlmLoss(theta, y, X, w, bfun)
  G <- ee(theta, X, w, bfun, L0, J = FALSE)
  g <- G$g/n
  nmg1 <- nmg <- sqrt(sum(g^2))
  alpha <- alpha0
  conv <- FALSE

  if(display) cat("\n Stochastic gradient-based iterations: \n")
  if(display & n.it == 0) cat("Iter =   0",
                              "  alpha = ", formatC(alpha, digits = 9, width = 10, format = "f"),
                              " --- ||g||_2 = ", formatC(nmg, digits = 15, width = 16, format = "f"),
                              " --- max|g| = ", formatC(max(abs(g)), digits = 9, width = 10, format = "f"),"\n", sep = "")

  w0 <- w

  hist_gs <- 0
  for(i in 1:maxit){

    hist_gs <- hist_gs + g^2
    step <- g * alpha / sqrt(tol + hist_gs)
    new.theta <- theta - step

    sampleid <- sample(n, round(n*.5))
    w <- w0; w[sampleid] <- 0

    L1 <- qlmLoss(new.theta, y, X, w, bfun)
    G1 <- ee(new.theta, X, w, bfun, L1, J = FALSE)
    g1 <- G1$g/n
    nmg1 <- sqrt(sum(g1^2)); if(is.na(nmg1)){nmg1 <- Inf}

    g <- g1; L0 <- L1; G <- G1; theta <- new.theta; nmg <- nmg1

    if(display) cat("Iter =",formatC(i + n.it, digits = 0, width = 4, format = "f"),
                    "  alpha = ", formatC(alpha, digits = 9, width = 10, format = "f"),
                    " --- ||g||_2 = ", formatC(nmg, digits = 15, width = 16, format = "f"),
                    " --- max|g| = ", formatC(max(abs(g)), digits = 9, width = 10, format = "f"),"\n", sep = "")

  }

  L0 <- qlmLoss(theta, y, X, w0, bfun)
  G <- ee(theta, X, w, bfun, L0, J = FALSE)
  g <- G$g/n
  nmg <- sqrt(sum(g^2))

  list(theta = theta, conv = TRUE, n.it = i, nmg = nmg, g = g, G = G, alpha = alpha)
}

Qlm.gs.internal <- function(theta, type, tol, maxit, alpha0, ee, display, n.it, y, X, w, bfun){

  n <- length(y)
  L0 <- qlmLoss(theta, y, X, w, bfun)
  G <- ee(theta, X, w, bfun, L0, J = FALSE)
  g <- G$g/n
  nmg1 <- nmg <- sqrt(sum(g^2))
  alpha <- alpha0
  conv <- FALSE

  if(display) cat("\n Gradient-based iterations: \n")
  if(display & n.it == 0) cat("Iter =   0",
                              "  alpha = ", formatC(alpha, digits = 9, width = 10, format = "f"),
                              " --- ||g||_2 = ", formatC(nmg, digits = 15, width = 16, format = "f"),
                              " --- max|g| = ", formatC(max(abs(g)), digits = 9, width = 10, format = "f"),"\n", sep = "")

  for(i in 1:maxit){
    if(i > 1 & (conv | nmg < tol)){break}
    cond <- FALSE

    while(!cond){
      step <- g * alpha
      if((max(abs(step)) < tol & i > (maxit/2)) | alpha < 1e-12){conv <- TRUE; break}
      new.theta <- theta - step
      L1 <- qlmLoss(new.theta, y, X, w, bfun)
      G1 <- ee(new.theta, X, w, bfun, L1, J = FALSE)
      g1 <- G1$g/n
      nmg1 <- sqrt(sum(g1^2)); if(is.na(nmg1)){nmg1 <- Inf}
      cond <- (nmg1 < nmg)
      alpha <- alpha*0.5
    }

    if(conv){break}
    alpha.print <- alpha*2

    g <- g1; L0 <- L1; G <- G1; theta <- new.theta; nmg <- nmg1
    alpha <- min(alpha*3, 1)

    if(display) cat("Iter =",formatC(i + n.it, digits = 0, width = 4, format = "f"),
                    "  alpha = ", formatC(alpha.print, digits = 9, width = 10, format = "f"),
                    " --- ||g||_2 = ", formatC(nmg, digits = 15, width = 16, format = "f"),
                    " --- max|g| = ", formatC(max(abs(g)), digits = 9, width = 10, format = "f"),"\n", sep = "")

  }

  list(theta = theta, conv = conv, n.it = i, nmg = nmg, g = g, L0 = L0, G = G, alpha = alpha*2)
}

Qlm.gs <- function(theta, type, tol, maxit, alpha0, ee, display, y, X, w, bfun){

  n <- length(y)
  model.type <- "Qlm"

  # Use gs, check that the result corresponds to an invertible jacobian
  # if(display) cat(" Preliminary gradient-based iterations: \n")
  if(display) cat(" Preliminary iterations: \n")
  fit0.ok <- n.it <- FALSE
  count <- 0
  while(!fit0.ok & count <= 5){
    count <- count + 1
    fit0 <- Qlm.sgs.internal(theta = theta, type = type, tol = tol, maxit = maxit,
                            alpha0 = alpha0, ee = ee, display = display, n.it = n.it,
                            y = y, X = X, w = w, bfun = bfun)
    fit0 <- Qlm.gs.internal(theta = fit0$theta, type = type, tol = tol, maxit = maxit,
                             alpha0 = alpha0, ee = ee, display = display, n.it = n.it,
                             y = y, X = X, w = w, bfun = bfun)

    # Compute J. If in "Qest", use a new eps (consistently re-calculate g, and not only J, with the new eps).
    G <- ee(fit0$theta, X, w, bfun, J = TRUE, EE = fit0$L0)
    H <- G$H/n

    # Check that J is def+
    Hinv <- invJ(H, type)
    fit0.ok <- !Hinv$err
    alpha0 <- min(alpha0 * 5, 1)
    # theta <- fit0$theta
    # alpha0 <- min(fit0$alpha, 1)
    n.it <- n.it + fit0$n.it - 1
  }
  if(!fit0.ok && fit0$conv){stop("Could not find theta s.t. J(theta) is def+")} # The user will not see this

  fit0$H <- H
  fit0$Hinv <- Hinv$Jinv
  fit0$G <- G
  fit0
}

Qlm.newton <- function(theta, type = "u", tol, maxit, safeit, alpha0, display, y, X, w, bfun){

  n <- length(y)
  model.type <- "Qlm"
  ee <- Qlm.ee.u

  # Preliminary gradient search
  fit0 <- Qlm.gs(theta = theta, type = type, tol = tol * 100, maxit = safeit,
                  alpha0 = alpha0, ee = ee, display = display, y = y, X = X,
                  w = w, bfun = bfun)

  new.theta <- theta <- fit0$theta
  g1 <- g <- fit0$g; H <- fit0$H; Hinv <- fit0$Hinv; G1 <- G <- fit0$G; L1 <- L0 <- fit0$L0
  nmg1 <- nmg <- fit0$nmg
  alpha <- 1.0 #min(fit0$alpha*10, 1)
  fit.ok <- TRUE # Just a reminder that NR starts at a value of theta where J is def+

  conv <- FALSE
  if(display) cat("\n Newton-Raphson: \n")
  for(i in 1:maxit){

    if(conv | nmg < tol){break}

    delta <- chol2inv(Hinv) %*% g
    cond <- FALSE
    while(!cond){
      step <- c(alpha*delta)
      if(max(abs(step)) < tol){conv <- TRUE; break}
      new.theta <- theta - step
      L1 <- qlmLoss(new.theta, y, X, w, bfun)
      G1 <- ee(new.theta, X, w, bfun, L1, J = FALSE)
      g1 <- G1$g/n; nmg1 <- sqrt(sum(g1^2)); if(is.na(nmg1)){nmg1 <- Inf}
      cond <- (nmg1 < nmg)
      alpha <- alpha*0.5
    }

    alpha.print <- alpha*2

    G1 <- ee(new.theta, X, w, bfun, L1, J = TRUE); H1 <- G1$H/n
    H1inv <- invJ(H1, type)
    fit.ok <- !H1inv$err

    if(conv){break}
    if(fit.ok){
      g <- g1; H <- H1; Hinv <- H1inv$Jinv; L0 <- L1; G <- G1; theta <- new.theta; nmg <- nmg1
      alpha <- min(alpha*3, 1)
    }
    else{alpha <- alpha*0.5} # In practice, this will repeat the iteration with alpha_t = 0.25*alpha_t-1

    if(display) cat("Iter =",formatC(i, digits = 0, width = 4, format = "f"),
                    "  alpha = ", formatC(alpha.print, digits = 9, width = 10, format = "f"),
                    " --- ||g||_2 = ", formatC(nmg, digits = 15, width = 16, format = "f"),
                    " --- max|g| = ", formatC(max(abs(g)), digits = 9, width = 10, format = "f"),"\n",sep = "")
  }

  conv <- (i < maxit)
  if(display && conv){cat("Converged.\n\n")}

  list(theta = new.theta, converged = conv, n.it = i - conv,
       ee = g1, ee.i = G1$g.i, L0 = L1, hessian = H1)
}




