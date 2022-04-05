# Note: start is for beta only
Qcoxph <- function(formula, weights, start, data, knots, wtau = NULL, control = Qcoxph.control(), ...){

  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "weights", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.response(mf)

  if(any((w <- model.weights(mf)) < 0)){stop("negative 'weights'")}
  if (!is.null(w) && !is.numeric(w)) stop("'weights' must be a numeric vector")
  if(is.null(w)){w <- rep.int(1, nrow(mf)); alarm <- FALSE}
  else{
    alarm <- (w == 0)
    sel <- which(!alarm)
    mf <- mf[sel,]
    Y <- Y[sel, ]
    w <- w[sel]
    w <- w / mean(w)
  }
  if(any(alarm)){warning("observations with null weight will be dropped")}
  if((n <- nrow(mf)) == 0){stop("zero non-NA cases")}

  Terms <- attr(mf, "terms")
  X <- model.matrix(Terms, mf, contrasts.arg = NULL)
  Xatt <- attributes(X)
  xdrop <- Xatt$assign %in% 0
  X <- X[, !xdrop, drop = FALSE]
  attr(X, "assign") <- Xatt$assign[!xdrop]
  attr(X, "contrasts") <- Xatt$contrasts
  nms <- colnames(X)
  if(is.null(nms)) nms <- paste0("Theta", seq_len(ncol(X)))

  # Response variable
  type <- attributes(Y)$type
  zyd <- cbind(Y)

  if(is.null(type)){
    y <- zyd[, 1]
    z <- rep.int(-Inf, n)
    d <- rep.int(1, n)
    type <- "u"
  }
  else if (type == "right"){
    y <- zyd[, 1]
    z <- rep.int(-Inf, n)
    d <- zyd[, 2]
    type <- (if (any(d == 0)) "c" else "u")
  }
  else if(type == "counting"){
    z <- zyd[, 1]
    y <- zyd[, 2]
    d <- zyd[, 3]
    type <- "ct"
    if(all(z < min(y))){type <- (if(any(d == 0)) "c" else "u")}
  }
  if(!(any(d == 1))) {stop("all data are censored")}
  attributes(type)$model <- "Cox"

  # tau, wtau, etc. Note that r, R, and iR are functions of taustar = log(-log(1 - tau)), for better accuracy

  tauL <- 1e-10; tauR <- 1 - 1e-10
  taustarL <- log(-log(1 - tauL)); taustarR <- log(-log(1 - tauR))
  wfit <- (!is.null(wtau)) # is there a weight w(tau)?

  if(!wfit){
    wfun <- function(tau){tau*0 + 1}
    Wfun <- I
    iWfun <- function(tau){0.5*tau^2}
    w1fun <- function(tau){0*tau}

    rfun <- function(taustar, ...){exp(taustar)} # -log(1 - tau)
    Rfun <- function(taustar){
      etaustar <- exp(taustar) # -log(1 - tau)
      tau1 <- exp(-etaustar) # 1 - tau
      out <- 1 - tau1*(1 + etaustar) # tau + (1 - tau)*log(1 - tau)
      pmax0(out) # should be always > 0, but it is not due to numerical issues
    }
    iRfun <- function(taustar){
      etaustar <- exp(taustar) # -log(1 - tau)
      tau1 <- exp(-etaustar) # 1 - tau
      tau <- 1 - tau1
      out <- 0.75*tau^2 - 0.5*(tau - (tau1^2)*etaustar) # 0.75*tau^2 - 0.5*(tau + (1 - tau)^2*log(1 - tau))
      pmax0(out) # should be always > 0, but it is not due to numerical issues
    }
    r1fun <- function(tau, ...){1/(1 - tau)}
  }
  else{
    # Notes:
    # w(tau) and r(tau) are not replaced by spline functions, for max precision.
    # w1(tau) and r1(tau) may not be monotone: I use approxfun, which is much safer (and faster) than splinefun.
    # For monotone functions, such as W(tau) and R(tau), splinefun(method = "hyman") is the best option.
    # w1(tau) and r1(tau) are computed using both left- and right-derivative, and choosing the side with min abs value.
      # This protects the estimates of the standard errors in case w(tau) is not continuous.

    ntau <- 4999
    tau <- seq(tauL, tauR, length = ntau); dtau <- (tau[2] - tau[1])*0.5
    taustar <- seq(taustarL, taustarR, length = ntau); dtaustar <- (taustar[2] - taustar[1])*0.5
    h <- exp(taustar - exp(taustar)) #

    #######################

    wfun <- wtau
    wtau <- wfun(tau, ...)
    if(any(is.na(wtau)) || any(wtau < 0) | any(!is.finite(c(wtau, wfun(0, ...), wfun(1, ...)))))
      {stop("'wtau(tau)' must return a non-negative, finite value for 0 <= tau <= 1")}

    #######################

    Wtau <- cumsum(wtau[1:(ntau - 1)] + wtau[2:ntau])*dtau
    iWtau <- cumsum(Wtau[1:(ntau - 2)] + Wtau[2:(ntau - 1)])*dtau
    w1tau <- (wtau[-1] - wtau[-ntau])/(2*dtau) # dtau was divided by 2
    Wfun <- splinefun(c(0,tau[-1]), c(0, Wtau), method = "hyman")
    iWfun <- splinefun(c(0,tau[-(1:2)]), c(0, iWtau), method = "hyman")
    w1fun <- approxfun(tau[-1], minabs(w1tau, c(w1tau[-1],Inf)), rule = 2)

    #######################
    #######################

    rfun <- function(taustar, ...){
      etaustar <- exp(taustar) # -log(1 - tau)
      tau1 <- exp(-etaustar) # 1 - tau
      etaustar*wfun(1 - tau1, ...) # -log(1 - tau)*wfun(tau)
    }
    rtau <- rfun(taustar, ...)

    #######################

    rtauh <- rtau*h
    Rtau <- cumsum(rtauh[1:(ntau - 1)] + rtauh[2:ntau])*dtaustar
    Rtauh <- Rtau*h[-1]
    iRtau <- cumsum(Rtauh[1:(ntau - 2)] + Rtauh[2:(ntau - 1)])*dtaustar
    rtaubis <- rfun(log(-log(1 - tau)), ...); r1tau <- (rtaubis[-1] - rtaubis[-ntau])/(2*dtau) # dtau was divided by 2

    Rfun <- splinefun(c(taustarL,taustar[-1]), c(0, Rtau), method = "hyman")
    iRfun <- splinefun(c(taustarL,taustar[-(1:2)]), c(0, iRtau), method = "hyman")
    r1fun <- function(tau, ...){tau <- pmin(tau, 1-1e-10); wfun(tau, ...)/(1 - tau) - log(1 - tau)*w1fun(tau)}
  }

  tau <- list(wfun = wfun, Wfun = Wfun, iWfun = iWfun, w1fun = w1fun,
              rfun = rfun, Rfun = Rfun, iRfun = iRfun, r1fun = r1fun,
              tauL = tauL, tauR = tauR, taustarL = taustarL, taustarR = taustarR,
              opt = list(...))

  # knots.
  # If 'knots' is missing, by default max(1, min(floor(sum(d)/30), 5)) internal knots are used.
  # If 'knots' is a scalar, its value is used to determine the number of internal knots
  # (knots = 0 is allowed, and fit an exponential model).
  # In both cases, the knots are positioned at the empirical quantiles of the observed events.

  # If 'knots' is a vector, it is used to identify the exact position of *all* knots,
  # including boundaries (so 'knots' should be a vector of at least two elements).

  ry <- range(y)
  if(missing(knots)){
    k <- max(1, min(floor(sum(d)/30), 3)) # n. of internal knots
    knots <- quantile(y[d == 1], (1:k)/(k + 1))
    knots <- c(ry[1] - 1e-6, knots, ry[2] + 1e-6)
    attr(knots, "user-defined") <- FALSE
  }
  else if(length(knots) == 1){
    knots <- (if(knots > 0) quantile(y[d == 1], (1:knots)/(knots + 1)) else NULL)
    knots <- c(ry[1] - 1e-6, knots, ry[2] + 1e-6)
    attr(knots, "user-defined") <- FALSE
  }
  else{
    knots <- sort(unique(knots))
    if(any(!is.finite(knots))){stop("All knots must be finite")}
    if(knots[1] < min(y) | knots[length(knots)] > max(y)){stop("All y values must be between the first and last knot")}
    # knots <- c(ry[1] - 1e-6, knots, ry[2] + 1e-6)
    test <- cut(y[d == 1], knots, include.lowest = TRUE)
    if(any(table(test) == 0)){stop("There must be some event between any two knots")}
    attr(knots, "user-defined") <- TRUE
  }

  z <- pmax(z, knots[1]) # This is fundamental! For z < knots[1], plfcox(z,knots) would be negative, causing H0(z) < 0.

  #######################################################################################################
  # Starting points #####################################################################################
  #######################################################################################################

  X0 <- X; z0 <- z; y0 <- y; knots0 <- knots; Y0 <- Y
  scaleVars <- scalevars.Qcoxph(X,z,y,knots)
  X <- scaleVars$X; z <- scaleVars$z; y <- scaleVars$y; knots <- scaleVars$knots; Y <- data.frame(y,d)

  #### check singularities
  temp <- check.singularities(X, scaleVars)
  X <- temp$X; scaleVars <- temp$scaleVars; nok <- temp$nok; ok <- temp$ok; rank <- temp$rank

  #### starting points
  theta <- starting.points.Qcox(X, Y, n, w, mf, knots)
  if(!attr(knots, "user-defined")) {
    knots <- attr(theta, "knots")
    knots0 <- knots*scaleVars$Stats$My/10
  }
  if(!missing(start)){
    if(length(start) != ncol(X0)){stop("length(start) does not match the number of predictors")}
    if(any(is.na(start))){stop("NAs are not allowed in 'start'")}
    start <- start[ok]
    theta[1:rank] <- start*scaleVars$Stats$sX
  }

  # Fit the model #########################################################################################

  # Control arguments
  npar <- length(theta)
  tol <- control$tol
  safeit <- control$safeit; if(is.null(safeit)){safeit <- 10 + 5*npar}
  maxit <- control$maxit; if(is.null(maxit)){maxit <- 50 + 25*npar}

  # safeit <- 10 + 2*npar
  alpha0 <- control$alpha0; if(is.null(alpha0)){alpha0 <- 0.01 / npar}
  display <- control$display

  fit.ok <- count <- FALSE
  while(!fit.ok){

    if(count > 0 & display){cat("Non-invertible Jacobian matrix. Restarting... \n\n")}
    fit <- try(Qest.newton(theta = theta, type = type, tol = tol, maxit = maxit, safeit = safeit,
      alpha0 = alpha0, display = display, eps = NA,
      z = z, y = y, d = d, w = w, X = X, knots = knots, tau = tau), silent = TRUE)

    # Sometimes this is ok, but covar is not
    fit.ok <- (!inherits(fit, "try-error"))
    if(fit.ok){
      theta0 <- descalecoef.Qcoxph(drop(fit$theta), scaleVars$Stats)
      covar <- try(Qcox.covar(theta0, z0, y0, d, X0[,ok,drop = FALSE], w, knots0, tau, type), silent = TRUE)
      fit.ok <- !(inherits(covar, "try-error"))
    }

    count <- count + 1
    if(count == 10){
      stop("Could not fit the model (non-invertible Jacobian matrix).
      Try changing the number or position of the knots.
      Verify if there is a sufficient number of events between any two knots.")
    }
    alpha0 <- alpha0/5
    maxit <- round(maxit * 1.25)
    safeit <- round(safeit * 1.5)
    theta[1:rank] <- theta[1:rank]*0.95 # Just shake it a little bit
  }

  # Finish ##################################################################################

  theta <- theta0; X <- X0; z <- z0; y <- y0; knots <- knots0; Y <- Y0
  fit$fy <- fit$fy/scaleVars$Stats$My
  fit$fz <- fit$fz/scaleVars$Stats$My

  L <- coxLoss(theta, z, y, d, X[, ok, drop = FALSE], w, knots, tau, type, fit$Fy, fit$Fz)
  attr(L, "npar") <- npar
  if(type != "u"){attr(L, "message") <- "Not the function being minimized"}

  # Output #########################################################

  npar0 <- ncol(X)
  coefficients1 <- rep(NA, npar0); var1 <- matrix(0, npar0, npar0)
  coefficients1[ok] <- theta[1:rank]; coefficients2 <- theta[-c(1:rank)]
  var <- covar; var1[ok,ok] <- var[1:rank, 1:rank,drop = FALSE]; var2 <- var[-c(1:rank),-c(1:rank), drop=FALSE]
  means <- colMeans(X)
  linear.predictors <- drop(X[,ok,drop = FALSE] %*% cbind(coefficients1[ok])) - sum(coefficients1[ok] * means[ok])
  score0 <- exp(linear.predictors)
  by <- plfcox(y, knots, deriv = 0); h0 <- c(by %*% exp(coefficients2)) * exp(sum(coefficients1[ok] * means[ok]));
  resid <- 1 * (d > 0) - h0 * score0
  concordance <- concordancefit(Y, linear.predictors, NULL, NULL, reverse = TRUE, timefix = FALSE)
  concordance <- c(concordance$count, concordance = concordance$concordance, std = sqrt(concordance$var))

  names(coefficients1) <- colnames(var1) <- rownames(var1) <- nms
  names(linear.predictors) <- names(resid) <- 1:n

  internal <- list(z = z, y = y, d = d, w = w, data = NULL, type = type, tau = tau, knots = knots, gamma = coefficients2,
                   var_gamma = var2, ee = fit$ee, jacobian = fit$jacobian, CDF = fit$Fy, PDF = fit$fy, converged = fit$converged)
  if(!missing(data)) internal$data <- data
  out <- list(coefficients = coefficients1, var = var1, obj.function = L,
              iter = fit$n.it, linear.predictors = linear.predictors,
              residuals = resid, means = means, n = n, nevent = sum(d),
              terms = Terms, assign = attrassign(X, Terms), concordance = concordance, model = mf, X = X,
              y = Y, formula = formula, call = match.call(), internal = internal)

  class(out) <- c("Qcoxph", "coxph", "Qest")
  return(out)
}


