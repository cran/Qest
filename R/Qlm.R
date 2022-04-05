# start = beta only, sigma is not included.

Qlm <- function(formula, data, subset, weights, na.action, start = NULL, contrasts = NULL,
                wtau = NULL, control = Qest.control(), ...){

  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf)

  if(any((w <- model.weights(mf)) < 0)){stop("negative 'weights'")}
  if (!is.null(w) && !is.numeric(w)) stop("'weights' must be a numeric vector")
  if(is.null(w)){w <- rep.int(1, nrow(mf)); alarm <- FALSE}
  else{
    alarm <- (w == 0)
    sel <- which(!alarm)
    mf <- mf[sel,]
    y <- y[sel]
    w <- w[sel]
    w <- w / mean(w)
  }
  if(any(alarm)){warning("observations with null weight will be dropped")}
  if((n <- nrow(mf)) == 0){stop("zero non-NA cases")}

  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = if(is.matrix(y)) matrix(NA_real_, 0, ncol(y)) else numeric(),
              residuals = y, fitted.values = 0 * y, weights = w, rank = 0L,
              df.residual = if(!is.null(w)) sum(w != 0) else if (is.matrix(y)) nrow(y) else length(y))
    if(!is.null(offset)) {
      z$fitted.values <- rep(0, n)
      z$residuals <- y
    }
  }
  else{
    x <- model.matrix(mt, mf, contrasts)
    z <- Qlm.fit(y, x, w, start, wtau, control, ...)
  }
  z$na.action <- attr(mf, "na.action")
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  z$model <- mf
  class(z) <- c("Qlm", "lm", "Qest")
  z
}

Qlm.fit <- function(y, X, w = rep(1, nobs), start = NULL, wtau = NULL,
                    control = Qest.control(), ...){

  if(is.null(wtau)) wtau <- function(tau, ...) rep(1, length(tau))

  #### check singularities
  QR <- qr(X)
  nok <- (abs(diag(qr.R(QR))) < 1e-6)
  ok <- !nok
  nobs <- nrow(X); nvars <- ncol(X)
  rank <- QR$rank
  assign <- attr(X, "assign")
  QR <- structure(QR[c("qr", "qraux", "pivot", "rank")], class = "qr")

  xnames <- dimnames(X)[[2L]]
  ynames <- names(y)
  X <- X[, ok, drop = FALSE]

  #### scale vars
  X0 <- X; y0 <- y
  scaleVars <- scalevars.Qlm(X, y)
  X <- scaleVars$X; y <- scaleVars$y

  # starting points
  theta0 <- start.Qlm(X, y, w, start, ok, scaleVars$Stats)
  if(!is.null(start)){
    if(length(start) != ncol(X)){stop("length(start) does not match the number of predictors")}
    if(any(is.na(start))){stop("NAs are not allowed in 'start'")}

    ss <- scaleVars$Stats
    start[ss$int] <- start[ss$int] - ss$my*ss$is.int
    start <- start*10/(ss$My - ss$my)
    start[ss$vars] <- (start*ss$sX)[ss$vars]
    start[ss$int] <- start[ss$int] + sum((start*ss$mX/ss$sX)[ss$vars])
    theta0[-1] <- start
  }
  # m0 <- lm.wfit(X, y, w)
  # theta0 <- c(log(sd(m0$residuals)), m0$coefficients)
  # if(!is.null(start)){
  #   if(length(start) != ncol(X)){stop("length(start) does not match the number of predictors")}
  #   if(any(is.na(start))){stop("NAs are not allowed in 'start'")}
  #
  #   ss <- scaleVars$Stats
  #   start[ss$int] <- start[ss$int] - ss$my*ss$is.int
  #   start <- start*10/(ss$My - ss$my)
  #   start[ss$vars] <- (start*ss$sX)[ss$vars]
  #   start[ss$int] <- start[ss$int] + sum((start*ss$mX/ss$sX)[ss$vars])
  #   theta0[-1] <- start
  # }

  # assign arguments in control list
  tol <- control$tol
  # c_1 <- 1e-4
  # fac <- .5
  epsilon <- tol
  maxit <- control$maxit; if(is.null(maxit)){maxit <- 50 + 25*rank}
  safeit <- control$safeit; if(is.null(safeit)){safeit <- 10 + 5*rank}
  display <- control$display
  alpha0 <- control$alpha0; if(is.null(alpha0)){alpha0 <- 0.01  / length(theta0)} # For gradient step
  bfun <- Qlm.bfun(wtau, ...)
  fit.ok <- count <- FALSE

  while(!fit.ok){
    if(count > 0 & display){cat("\n\n Non-invertible Jacobian matrix. Restarting... \n\n")}
    fit <- try(Qlm.newton(theta = theta0, type = "u", tol = tol, maxit = maxit, safeit = safeit,
                          alpha0 = alpha0, display = display, y = y, X = X, w = w, bfun = bfun),
               silent = FALSE)

    fit.ok <- (!inherits(fit, "try-error"))
    if(fit.ok){
      theta1 <- descalecoef.Qlm(fit$theta, scaleVars$Stats)
      L1<- qlmLoss(theta1, y0, X0, w, bfun)
      temp <- Qlm.ee.u(theta1, X0, w, bfun, L1, J = TRUE)
      g <- temp$g / nobs
      H <- temp$H
      covar <- try(Qlm.covar(temp$g.i, w, H), silent = TRUE)
      fit.ok <- !(inherits(covar, "try-error"))
    }
    else{
      count <- count + 1
      if(count == 5){stop("Could not fit the model (non-invertible Jacobian matrix).
        Is the model identifiable? Try with better starting points.
        Please, do not use Qest to fit discrete distributions.")}
      alpha0 <- alpha0 / 5
      maxit <- round(maxit * 1.25)
      safeit <- round(safeit * 1.5)
      theta <- fit$theta - 1e-3
    }
  }

  # L0 <- qlmLoss(theta0, y, X, w, bfun)
  # temp <- Qlm.ee.u(theta0, X, w, bfun, L0, J = FALSE)
  # g <- temp$g
  #
  # alpha0 <-  control$alpha0; if(is.null(alpha0)){alpha0 <- 0.05}
  # alpha <- alpha0
  # for(i in 1:safeit){
  #   Step <- -g
  #   dtg <- sum(Step * g)
  #   cond <- TRUE
  #   while(cond){
  #     theta1 <- theta0 + alpha*Step
  #     L1 <- qlmLoss(theta1, y, X, w, bfun)
  #     cond <- (L1$Loss > L0$Loss + c_1*alpha*dtg)
  #     alpha <- fac * alpha
  #     if(alpha <= tol) break
  #   }
  #   temp <- Qlm.ee.u(theta1, X, w, bfun, L1, J = FALSE)
  #   g <- temp$g
  #   nmg <- max(abs(g))
  #
  #   if(display & i == 1) cat(" Preliminary gradient-based iterations: \n")
  #   if(display) cat("Iter = ", formatC(i, digits = 0, width = 3, format = "d"),
  #                   " alpha = ", formatC(alpha*2, digits = 9, width = 10, format = "f"),
  #                   " --- ||g||^2 = ", formatC(sum(g^2), digits = 15, width = 16, format = "f"),
  #                   " --- max |g| = ", formatC(nmg, digits = 9, width = 10, format = "f"),"\n",sep = "")
  #
  #   if(nmg < tol){break}
  #   alpha <- alpha0
  #   L0 <- L1; theta0 <- theta1
  # }
  #
  # temp <- Qlm.ee.u(theta1, X, w, bfun, L1, J = TRUE)
  # g <- temp$g; H <- temp$H
  #
  # alpha0 <-  alpha <- 1
  # for(i in 1:maxit){
  #   Step <- -(g %*% solve(H + .0001*diag(ncol(H))))
  #   dtg <- sum(Step * g)
  #   cond <- TRUE
  #   while(cond){
  #     theta1 <- theta0 + alpha*Step
  #     L1 <- qlmLoss(theta1, y, X, w, bfun)
  #     cond <- (L1$Loss > L0$Loss + c_1*alpha*dtg)
  #     alpha <- fac * alpha
  #     if(alpha <= tol) break
  #   }
  #   temp <- Qlm.ee.u(theta1, X, w, bfun, L1, J = TRUE)
  #   g <- temp$g; H <- temp$H
  #   nmg <- max(abs(g))
  #
  #   if(display & i == 1) cat("\n Newton-Raphson:: \n")
  #   if(display) cat("Iter = ", formatC(i, digits = 0, width = 3, format = "d"),
  #                   " alpha = ", formatC(alpha*2, digits = 9, width = 10, format = "f"),
  #                   " --- ||g||^2 = ", formatC(sum(g^2), digits = 15, width = 16, format = "f"),
  #                   " --- max |g| = ", formatC(nmg, digits = 9, width = 10, format = "f"),"\n",sep = "")
  #
  #   if(nmg < tol){break}
  #   alpha <- alpha0
  #   L0 <- L1; theta0 <- theta1
  # }
#
#   theta1 <- descalecoef.Qlm(theta1, scaleVars$Stats)
#   L1 <- qlmLoss(theta1, y0, X0, w, bfun)
#   temp <- Qlm.ee.u(theta1, X0, w, bfun, L1, J = TRUE)
#   g <- temp$g
#   H <- temp$H

  theta <- rep(NA, nvars+1); theta[c(TRUE, ok)] <- c(theta1)
  # covar <- Qlm.covar(temp$g.i, w, H)
  gradient <- g[-1]
  covarPhi <- covar[1, 1]; covar <- covar[-1, -1, drop = FALSE]
  H <- H[-1, -1, drop = FALSE]

  coefficients <- theta[-1]
  dispersion <- exp(theta[1])^2
  attr(dispersion, "std.err") <- sqrt(covarPhi)

  se <- sqrt(diag(covar))
  std.errs <- coefficients; std.errs[ok] <- se

  fitted.values <- drop(X0 %*% theta1[-1])
  residuals <- drop(y0 - fitted.values)
  df.residual <- nobs - rank

  colnames(covar) <- rownames(covar) <- names(gradient) <- xnames[ok]
  names(coefficients) <- names(std.errs) <- xnames
  names(fitted.values) <- names(residuals) <- ynames

  list(coefficients = coefficients, std.errs = std.errs, covar = covar, dispersion = dispersion, residuals = residuals,
       rank = rank, fitted.values = fitted.values, assign = assign, qr = QR, df.residual = df.residual,
       obj.function = L1$Loss, gradient = gradient, hessian = H, CDF = L1$p.i, convergence = fit$converged,
       n.it = fit$n.it, control = control)
}

