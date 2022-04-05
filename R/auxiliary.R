expit <- function(x){1/(1 + exp(-x))}
logit <- function(x){log(x) - log(1 - x)}
pmax0 <- function (x){(x + abs(x))/2}
dlist <- function(x1,x2){for(j in 1:length(x1)){x1[[j]] <- x1[[j]] - x2[[j]]}; x1}
Ltau <- function(opt, tau){opt$tau <- tau; opt}
format.perc <- function (probs, digits) paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")
minabs <- function(x1,x2){
  i <- which(abs(x1) < abs(x2))
  out <- x2; out[i] <- x1[i]
  out
}
rq.fit.br2 <- function (x, y, tau = 0.5) {
  tol <- .Machine$double.eps^(2/3)
  eps <- tol
  big <- .Machine$double.xmax
  x <- as.matrix(x)
  p <- ncol(x)
  n <- nrow(x)
  ny <- NCOL(y)
  nsol <- 2
  ndsol <- 2
  lci1 <- FALSE
  qn <- rep(0, p)
  cutoff <- 0

  z <- .Fortran("rqbr", as.integer(n), as.integer(p), as.integer(n + 5), as.integer(p + 3), as.integer(p + 4),
                as.double(x), as.double(y), as.double(tau), as.double(tol), flag = as.integer(1),
                coef = double(p), resid = double(n), integer(n), double((n + 5) * (p + 4)), double(n),
                as.integer(nsol), as.integer(ndsol), sol = double((p + 3) * nsol), dsol = double(n * ndsol),
                lsol = as.integer(0), h = integer(p * nsol), qn = as.double(qn),
                cutoff = as.double(cutoff), ci = double(4 * p), tnmat = double(4 * p), as.double(big), as.logical(lci1))

  coef <- z$coef
  dual <- z$dsol[1:n]
  names(coef) <- dimnames(x)[[2]]
  return(list(coefficients = coef, x = x, y = y, residuals = y - x %*% z$coef, dual = dual))
}

# 1) if delta.right = 0: zero weight for higher quantiles
# 2) if delta.left = 0: zero weight for lower quantiles
# 3) if delta.right != 0 & delta.left != 0: zero weight for lower and higher quantiles
# 4) if smooth = TRUE, all the previous version are smoothed
wtrunc <- function(tau, delta.left = 0.05, delta.right = 0.05, smooth = FALSE, sigma = 0.01){
  w1 <- if(delta.right !=0){if(smooth) 1 - pnorm((tau - (1 - delta.right))/sigma) else as.numeric(tau <= 1 - delta.right)} else 1
  w2 <- if(delta.left != 0){if(smooth) pnorm((tau - delta.left)/sigma) else as.numeric(tau >= delta.left)} else 1
  w1*w2
}


########### Qest ##################################################

callQ <- function(Q, theta, tau, data){

  Qoptions <- attr(Q, "Qoptions")
  Qoptions$theta <- theta
  Qoptions$tau <- tau
  if(is.null(attr(Q, "X"))) {
    Qoptions$data <- data
    do.call(Q, Qoptions)
  }
  else {
    Qoptions$X <- attr(Q, "X")
    do.call(Q, Qoptions)$Q
  }
}

callwtau <- function(wtau, tau, opt){
  opt$tau <- tau
  do.call(wtau, opt)
}

buildTau <- function(ntau, wtau = NULL, nobs, wtauoptions = NULL) {
  # tau and wtau. To see the reason behind omicron*1.1, see the comment to findAp
  omicron <- 0.0005; tauL <- omicron; tauR <- 1 - omicron
  tau <- pbeta(seq.int(qbeta(omicron*1.1,2,2), qbeta((1 - omicron*1.1),2,2), length.out = ntau - 2), 2, 2)
  tau <- c(tauL, tau, tauR)
  dtau <- 0.5*(c(tau[1], tau[2:ntau] - tau[1:(ntau - 1)])) # The 0.5 is for integral with trapezoids

  wfun <- wtau
  wfit <- (!is.null(wfun)) # is there a weight w(tau)?
  if(!wfit){wfun <- function(tau, ...){rep(1,length(tau))}}
  wtau <- callwtau(wfun, tau, wtauoptions)
  w01 <- callwtau(wfun, 0:1, wtauoptions)

  if(wfit){
    if(any(is.na(wtau)) || any(wtau < 0) | any(!is.finite(c(wtau, w01))))
    {stop("'wtau(tau)' must return a non-negative, finite value for 0 <= tau <= 1")}
    wtau <- pmax(wtau, 1e-12)
  }

  TAU <- t(matrix(tau, ntau, nobs))
  tau <- list(
    tau = tau, dtau = dtau, ntau = ntau, TAU = TAU,
    tau1 = c(0, tau, 1), dtau1 = c(dtau, 1 - tau[ntau]), ntau1 = ntau + 2,
    dtau_long = rep(dtau, rep.int(nobs, ntau)), wtau_long = rep(wtau, rep.int(nobs, ntau)),
    wtau = wtau, wfun = wfun,
    tauL = tauL, tauR = tauR,
    opt = wtauoptions)
  tau
}

start.Qest <- function(z, y, d, x, w, tau, Q, opt, start, data, type, control) {
  nobs <- length(y)

  #### check singularities
  QR <- qr(x)
  nok <- (abs(diag(qr.R(QR))) < 1e-6)
  ok <- !nok
  x <- x[, ok, drop = FALSE]
  if(sum(nok) > 0) warning("\nSome variables have been dropped out\n")

  argQ <- formalArgs(Q)
  Qoptions <- opt[names(opt) %in% argQ[-(1:3)]]

  # Q function, starting values
  if(any(argQ[1:3] != c("theta", "tau", "data")))
  {stop("Q must have argument names 'theta', 'tau', 'data', in this order, followed by other optional arguments")}
  attr(Q, "Qoptions") <- Qoptions

  if(missing(start) || is.null(start))
  {stop(
    "\n
      Please provide starting values!
      Note that length(start) is used to determine the number of parameters.
      If you supply start = c(NA,NA, ...), Qest will use start = c(0,0, ...)."
  )}
  start <- start[ok]
  npar <- length(start)
  nostart <- any(is.na(start))
  start[is.na(start)] <- 0

  testQ <- callQ(Q, start, t(matrix((1:9)/10, 9, nobs)), data)
  if(!inherits(testQ, "matrix") || any(dim(testQ) != c(nobs, 9)))
  {stop("Q(theta, tau, data) must return a matrix when tau is a matrix")}
  if(any(is.na(testQ) | !is.finite(testQ)))
  {stop("Q(start, tau, data) does not return a finite value")}

  #######################################################################################################
  # Starting points #####################################################################################
  #######################################################################################################

  # We compute a "good" starting points, using the provided "start" (if available, otherwise 0,0,0...)
  # as starting points. If a "start" is provided, and contro$restart = FALSE, we believe your starting values and skip all this.
  if(nostart | control$restart){
    pch.fit <- getFromNamespace("pch.fit.ct", ns = "pch")
    predF.pch <- getFromNamespace("predF.pch", ns = "pch")
    m0 <- pch.fit(z, y, d, x, w)
    Fy <- 1 - predF.pch(m0, x, y)[,"Surv"]

    obj <- function(theta, tau, data, y, w, Q) 1/length(y)/var(y)*sum(w*(y - callQ(Q, theta, tau, data))^2)
    use <- (Fy > 0.1 & Fy < 0.9)

    # I use an amazing trick: if "start" is very bad, the distribution may be degenerated.
    # However, if gs finds an eps that "moves" some quantiles, the same eps can be used to move "start"
    # in a direction that may not be "right" (the eps has nothing to do with the gradient) but at least
    # yields a non-degenerated distribution.
    fit0 <- list(theta = start, n.it = 1, eps = 0)
    temp.fit <- suppressWarnings(nlm(obj, fit0$theta, tau = Fy[use], data = data[use,], y = y[use],
                                     w = w[use], Q = Q))
    temp.fit2 <- gs(fit0$theta - fit0$eps, obj, tau = Fy[use], data = data[use,],
                    y = y[use], w = w[use], Q = Q, maxit = 20 + 5*npar, tol = 1e-4)
    fit0$theta <- if(temp.fit$minimum > temp.fit2$f) temp.fit2$theta else temp.fit$estimate
    fit0 <- gs(fit0$theta - fit0$eps, obj, tau = Fy[use], data = data[use,],
               y = y[use], w = w[use], Q = Q, maxit = 20 + 5*npar, tol = 1e-4)
    # while(fit0$n.it == 1){
    # fit0 <- gs(fit0$theta - fit0$eps, obj, tau = Fy[use], data = data[use,],
    #            y = y[use], w = w[use], Q = Q, maxit = 20 + 5*npar, tol = 1e-4)
    # }
    theta <- fit0$theta
    eps <- fit0$eps
    alpha <- fit0$alpha/10
  }
  else{
    theta <- start
    eps <- choose_eps(Q, theta, y, data, eps0 = rep(0.001, npar), obj = 0.01)
    alpha <- min(eps) / npar / 100
  }
  list(z = z, y = y, x = x, tau = tau, Q = Q, theta = theta, eps = eps, alpha = alpha)
}

start.Qest.family <- function(z, y, d, x, w, tau, wtau, wtauoptions, Q, opt, start, data, type, control) {
  if(Q$family$family == "gaussian" & type == "u")
    warning("To fit linear models with Normal homoskedastic errors, please use Qlm")

  if(Q$family$family %in% c("uniform", "poisson") & type %in% c("c", "ct"))
    stop("Censored and truncated data not supported")

  nobs <- length(y)

  #### offset
  offset <- NULL
  if(Q$family$family == "poisson") {
    offset <- Q$offset
    if (!is.null(offset)) {
      if(!is.null(id <- attr(data, "na.action"))) offset <- offset[-id]
      if (length(offset) != nobs)
        stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                      length(offset), nobs), domain = NA)
    }
    if (is.null(offset)) offset <- rep.int(0, nobs)
    tau <- tau.pois(tau)
  }

  #### check singularities
  nms <- colnames(x)
  QR <- qr(x)
  nok <- (abs(diag(qr.R(QR))) < 1e-6)
  ok <- !nok
  names(ok) <- nms
  x <- x[, ok, drop = FALSE]
  if(sum(nok) > 0) warning("\nSome variables have been dropped out\n")
  attr(x, "offset") <- offset
  attr(x, "min") <- Q$min

  temp <- Q$scale(x, y, z)
  colnames(temp$Xc) <- colnames(temp$Stats$X) <- nms[ok]
  temp$nms <- nms
  x <- temp$Xc; y <- temp$yc; z <- temp$zc
  attr(x, "offset") <- offset

  theta <- Q$initialize(x, z, y, d, w, Q, start, tau, data, ok, temp$Stats)
  npar <- length(theta)

  Q0 <- Q$Q
  attr(Q0, "scale") <- Q$scale
  attr(Q0, "descale") <- Q$descale
  attr(Q0, "glmFamily") <- Q$family
  if(!is.null(Q$fix.Fy)) attr(Q0, "fix.Fy") <- Q$fix.Fy
  if(!is.null(Q$bfun)) tau$bfun <- Q$bfun(wtau, wtauoptions)
  attr(Q0, "X") <- x
  attr(Q0, "temp") <- temp
  attr(Q0, "ok") <- ok
  class(Q0) <- class(Q)
  Q <- Q0

  eps <- 1 #choose_eps(Q, theta, y, data, eps0 = rep(0.001, npar), obj = 0.01)
  alpha <- 0.01 / npar

  list(z = z, y = y, x = x, tau = tau, Q = Q, theta = theta, eps = eps, alpha = alpha)
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

Qest.sgs.internal <- function(theta, type, tol, maxit, alpha0, ee, display, eps, n.it, ...){

  temp <- list(...)
  n <- length(temp$y)
  # G <- ee(theta, ..., J = FALSE, eps = eps)

  w0 <- temp$w
  sampleid <- sample(n, round(n*.5))
  w <- w0; w[sampleid] <- 0

  G <- ee(theta, eps = eps, temp[[1]], temp[[2]], temp[[3]], temp[[5]], w, temp[[6]], tau = temp$tau, J = FALSE)
  g <- G$g/n
  nmg1 <- nmg <- sqrt(sum(g^2))
  alpha <- alpha0
  conv <- FALSE

  if(display) cat("\n Stochastic gradient-based iterations: \n")
  if(display & n.it == 0) cat("Iter =   0",
                              "  alpha = ", formatC(alpha, digits = 9, width = 10, format = "f"),
                              " --- ||g||_2 = ", formatC(nmg, digits = 15, width = 16, format = "f"),
                              " --- max|g| = ", formatC(max(abs(g)), digits = 9, width = 10, format = "f"),"\n", sep = "")

  hist_gs <- 0
  for(i in 1:maxit){

    hist_gs <- hist_gs + g^2
    step <- g * alpha / sqrt(tol + hist_gs)
    new.theta <- theta - step

    G1 <- ee(new.theta, temp[[1]], temp[[2]], temp[[3]], temp[[5]], w, temp[[6]], tau = temp$tau, J = FALSE, eps = eps)
    g1 <- G1$g/n
    nmg1 <- sqrt(sum(g1^2)); if(is.na(nmg1)){nmg1 <- Inf}

    g <- g1; G <- G1; theta <- new.theta; nmg <- nmg1

    # print(c(nmg, nmg1))
    # print(rbind(new.theta, g, alpha / sqrt(tol + hist_gs)))

    if(display) cat("Iter =",formatC(i + n.it, digits = 0, width = 4, format = "f"),
                    "  alpha = ", formatC(alpha, digits = 9, width = 10, format = "f"),
                    " --- ||g||_2 = ", formatC(nmg, digits = 15, width = 16, format = "f"),
                    " --- max|g| = ", formatC(max(abs(g)), digits = 9, width = 10, format = "f"),"\n", sep = "")

    w0 <- temp$w
    sampleid <- sample(n, round(n*.5))
    w <- w0; w[sampleid] <- 0
    G <- ee(theta, eps = eps, temp[[1]], temp[[2]], temp[[3]], temp[[5]], w, temp[[6]], tau = temp$tau, J = FALSE)
    g <- G$g/n
    nmg <- sqrt(sum(g^2))
  }

  G <- ee(theta, temp[[1]], temp[[2]], temp[[3]], temp[[5]], temp[[4]], temp[[6]], tau = temp$tau, J = FALSE, eps = eps)
  g <- G$g/n
  nmg <- sqrt(sum(g^2))

  list(theta = theta, conv = TRUE, n.it = i, nmg = nmg, g = g, G = G, alpha = alpha)
}

Qest.gs.internal <- function(theta, type, tol, maxit, alpha0, ee, display, eps, n.it, ...){

  n <- length(list(...)$y)
  G <- ee(theta, ..., J = FALSE, eps = eps)
  g <- G$g/n
  nmg1 <- nmg <- sqrt(sum(g^2))
  alpha <- min(alpha0, 1)
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
      # if((max(abs(step)) < tol & i > (maxit/2)) | alpha < 1e-12){conv <- TRUE; break}
      if((nmg1 < tol) | alpha < 1e-12){conv <- TRUE; break}
      new.theta <- theta - step
      G1 <- ee(new.theta, ..., J = FALSE, eps = eps)
      g1 <- G1$g/n
      nmg1 <- sqrt(sum(g1^2)); if(is.na(nmg1)){nmg1 <- Inf}
      cond <- (nmg1 <= nmg)
      alpha <- alpha*0.5
    }

    if(conv){break}
    alpha.print <- alpha*2

    g <- g1; G <- G1; theta <- new.theta; nmg <- nmg1
    alpha <- min(alpha*3, 1)

    if(display) cat("Iter =",formatC(i + n.it, digits = 0, width = 4, format = "f"),
                    "  alpha = ", formatC(alpha.print, digits = 9, width = 10, format = "f"),
                    " --- ||g||_2 = ", formatC(nmg, digits = 15, width = 16, format = "f"),
                    " --- max|g| = ", formatC(max(abs(g)), digits = 9, width = 10, format = "f"),"\n", sep = "")

  }

  list(theta = theta, conv = conv, n.it = i, nmg = nmg, g = g, G = G, alpha = alpha*2)
}

Qest.gs <- function(theta, type, tol, maxit, alpha0, ee, display, eps, ...){

  n <- length(list(...)$y)
  model.type <- attributes(type)$model
  A <- list(...)

  # Use gs, check that the result corresponds to an invertible jacobian

  # print(theta)

  # if(display) cat(" Preliminary gradient-based iterations: \n")
  if(display) cat(" Preliminary iterations: \n")
  fit0.ok <- n.it <- FALSE
  count <- 0L
  # max.count <- 5L
  # alpha00 <- alpha0
  while(!fit0.ok & count <= 5){
    count <- count + 1L
    fit0 <- Qest.sgs.internal(theta = theta, type = type, tol = tol, maxit = maxit,
                              alpha0 = alpha0, ee = ee, display = display, eps = eps, n.it = n.it, ...)

    # print(fit0$theta)

    fit0 <- Qest.gs.internal(theta = fit0$theta, type = type, tol = tol, maxit = maxit,
                             alpha0 = alpha0, ee = ee, display = display, eps = eps, n.it = n.it, ...)

    # fit0 <- Qest.gs.internal(theta = theta, type = type, tol = tol, maxit = maxit,
    #                          alpha0 = alpha0, ee = ee, display = display, eps = eps, n.it = n.it, ...)

# browser()
    # pippoX <- seq(-2, 2, l=40)
    # pippoY <- seq(-5, 5, l=50)
    # pippo <- expand.grid(pippoX, pippoY)
    # sampleid <- sample(n, round(.5*n))
    # pippoZ <- sapply(seq.int(nrow(pippo)), function(i) ee(unlist(pippo[i,]), ..., J = FALSE, eps = .1)$g/n)
    # persp(pippoX, pippoY, matrix(pippoZ[1,], nrow=length(pippoX), ncol=length(pippoY)), ticktype = "detailed", theta = 45, phi = 0, shade = 0.5, xlab = "log(shape)", ylab = "scale", zlab = "ee(shape)")
    # persp(pippoX, pippoY, matrix(pippoZ[2,], nrow=length(pippoX), ncol=length(pippoY)), ticktype = "detailed", theta = 45, phi = 0, shade = 0.5, xlab = "log(shape)", ylab = "scale", zlab = "ee(scale)")
    #
    # temp <- matrix(pippoZ[1,], nrow=length(pippoX), ncol=length(pippoY), dimnames = list(round(pippoX, 2), round(pippoY, 2)))
    # heatmap(temp, scale = "none", Rowv = NA, Colv = NA, labRow = round(pippoX, 2), labCol = round(pippoY, 2), main = "log(shape)")
    # temp2 <- matrix(pippoZ[2,], nrow=length(pippoX), ncol=length(pippoY), dimnames = list(round(pippoX, 2), round(pippoY, 2)))
    # heatmap(temp2, scale = "none", Rowv = NA, Colv = NA, labRow = round(pippoX, 2), labCol = round(pippoY, 2), main = "scale")

    # pippo <- seq(-2, 10, l=100)
    # grad <- sapply(1:100, function(i) {ee(pippo[i], ..., J = FALSE, eps = .1)$g/n})
    # plot(pippo, grad)

    # print(fit0$theta)

    # Compute J. If in "Qest", use a new eps (consistently re-calculate g, and not only J, with the new eps).
    if(model.type == "Qest" & is.null(attr(list(...)$Q, "glmFamily")$family)){
      # if(!any(attr(list(...)$Q, "glmFamily")$family %in% c("Gamma", "poisson", "uniform", "gaussian")))
      eps <- choose_eps(A$Q, fit0$theta, A$y, A$data, eps0 = sign(eps)*pmax(abs(eps)/10, 1e-6), obj = 0.01)
      G <- ee(fit0$theta, ..., J = TRUE, eps = eps)
      J <- G$J/n; fit0$g <- G$g/n
    }
    else{
      G <- ee(fit0$theta, ..., J = TRUE, EE = fit0$G, eps = eps)
      J <- G$J/n
    }

    # Check that J is def+
    Jinv <- invJ(J, type)
    fit0.ok <- !Jinv$err # && ifelse(count >= max.count, TRUE, (fit0$alpha >= alpha00))
    # if(alpha0 == 1) break
    # alpha0 <- min(alpha0 * 5, 1)
    alpha0 <- alpha0 * ifelse(alpha0 >= 1, 2, 5)
    # if(!fit0.ok && fit0$conv){stop("Could not find theta s.t. J(theta) is def+")} # The user will not see this # && count >= max.count
    # theta <- fit0$theta
    # alpha0 <- min(fit0$alpha, 1)
    # if(alpha0 < tol) alpha0 <- alpha00
    # alpha0 <- max(fit0$alpha, alpha0)
    n.it <- n.it + fit0$n.it - 1
  }
  if(!fit0.ok && fit0$conv){stop("Could not find theta s.t. J(theta) is def+")} # The user will not see this # && count >= max.count
  if(!fit0.ok && !fit0$conv){stop("The gradient search algorithm didn't converge.
  Moreover, could not find theta s.t. J(theta) is def+")} # The user will not see this

  fit0$J <- J
  fit0$Jinv <- Jinv$Jinv
  fit0$eps <- eps
  fit0$G <- G
  fit0
}

Qest.newton <- function(theta, type, tol, maxit, safeit, alpha0, display, eps, ...){

  n <- length(list(...)$y)
  model.type <- attributes(type)$model

  tempF <- attr(list(...)$Q, "glmFamily")$family
  if(model.type == "Qest"){
    if(type == "u"){
      ee <- Qest.ee.u
      if(!is.null(tempF)) ee <- switch(tempF, "Gamma" = QestGamma.ee.u, "poisson" = QestPois.ee.u,
                                       "uniform" = QestUnif.ee.u, "gaussian" = QestNorm.ee.u)
    }
    else if(type == "c"){
      ee <- Qest.ee.c
      if(!is.null(tempF)) ee <- switch(tempF, "Gamma" = QestGamma.ee.c, "gaussian" = QestNorm.ee.c)
    }
    else{
      ee <- Qest.ee.ct
      if(!is.null(tempF)) ee <- switch(tempF, "Gamma" = QestGamma.ee.ct, "gaussian" = QestNorm.ee.ct)
    }
  }
  else{
    if(type != "ct"){ee <- QCox.ee.c}
    else{ee <- QCox.ee.ct}
  }

  # Preliminary gradient search
  fit0 <- Qest.gs(theta = theta, type = type, tol = tol * 100, maxit = safeit,
    alpha0 = alpha0, ee = ee, display = display, eps = eps, ...)

  new.theta <- theta <- fit0$theta
  g1 <- g <- fit0$g; J <- fit0$J; Jinv <- fit0$Jinv; G1 <- G <- fit0$G
  nmg1 <- nmg <- fit0$nmg
  alpha <- 1.0 #min(fit0$alpha*10, 1);
  eps <- fit0$eps
  fit.ok <- TRUE # Just a reminder that NR starts at a value of theta where J is def+

  conv <- FALSE
  if(display) cat("\n Newton-Raphson: \n")
  for(i in 1:maxit){

    if(conv | nmg < tol){break}

    delta <- (if(type == "u") chol2inv(Jinv) %*% g else qr.solve(Jinv) %*% g)
    cond <- FALSE
    while(!cond){
      step <- c(alpha*delta)
      if(max(abs(step)) < tol){conv <- TRUE; break}
      new.theta <- theta - step
      G1 <- ee(new.theta, ..., J = FALSE, eps = eps)
      g1 <- G1$g/n; nmg1 <- sqrt(sum(g1^2)); if(is.na(nmg1)){nmg1 <- Inf}
      cond <- (nmg1 <= nmg)
      alpha <- alpha*0.5
    }

    alpha.print <- alpha*2

    G1 <- ee(new.theta, ..., J = TRUE, EE = G1, eps = eps); J1 <- G1$J/n
    J1inv <- invJ(J1, type)
    fit.ok <- !J1inv$err

    if(conv){break}
    if(fit.ok){
      g <- g1; J <- J1; Jinv <- J1inv$Jinv; G <- G1; theta <- new.theta; nmg <- nmg1
      alpha <- min(alpha*3, 1)
    }
    else{alpha <- alpha*0.5} # In practice, this will repeat the iteration with alpha_t = 0.25*alpha_t-1

    if(display) cat("Iter =",formatC(i, digits = 0, width = 4, format = "f"),
                    "  alpha = ", formatC(alpha.print, digits = 9, width = 10, format = "f"),
                    " --- ||g||_2 = ", formatC(nmg, digits = 15, width = 16, format = "f"),
                    " --- max|g| = ", formatC(max(abs(g)), digits = 9, width = 10, format = "f"),"\n",sep = "")
  }

  Q <- list(...)$Q
  if(inherits(Q, "Qfamily")) {
    temp <- attr(Q, "temp")$Stats
    attr(Q, "X") <- temp$X
    theta <- attr(Q, "descale")(theta, temp)
    G1 <- ee(theta, eps, temp$z, temp$y, list(...)$d, Q, list(...)$w, list(...)$data, list(...)$tau, J = TRUE)
    J <- G1$J/n; Jinv <- invJ(J, type); fit.ok <- !Jinv$err; g <- G1$g/n; G <- G1; nmg <- sum(g^2)

    if(attr(Q, "glmFamily")$family %in% c("Gamma", "uniform", "gaussian")) {
      G$Fy <- attr(Q, "fix.Fy")(list(Fy=G$Fy), theta, list(...)$tau, temp$y, temp$X, Q)
      if(type == "ct") G$Fz <- attr(Q, "fix.Fy")(list(Fy=G$Fz), theta, list(...)$tau, temp$y, temp$X, Q)
    }
  }

  conv <- (i < maxit)
  if(display && conv){cat("Converged.\n\n")}

  list(theta = theta, converged = conv, n.it = i - conv,
       ee = g, ee.i = G$g.i, jacobian = J, eps = eps,
       Fy = G$Fy, fy = G$fy, Fz = G$Fz, fz = G$fz)
}




