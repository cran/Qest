
Qest <- function(formula, Q, weights, start, data, ntau = 199, wtau = NULL,
                 control = Qest.control(), ...){

  cl <- match.call()
  okfamily <- c("Qnorm", "Qgamma", "Qpois", "Qunif")
  if (is.character(Q))
    if(all(!(Q %in% okfamily))){print(Q); stop("'Qfamily' not recognized")}
  if (is.character(Q)) Q <- get(Q, mode = "function", envir = parent.frame())
  if (is.function(Q)) Q <- if(inherits(try(Q(), silent=TRUE), "Qfamily")) Q() else Q

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  if(any((w <- model.weights(mf)) < 0)){stop("negative 'weights'")}
  if (!is.null(w) && !is.numeric(w)) stop("'weights' must be a numeric vector")
  if(is.null(w)){w <- rep.int(1, nrow(mf)); alarm <- FALSE}
  else{
    alarm <- (w == 0)
    sel <- which(!alarm)
    mf <- mf[sel,]
    w <- w[sel]
    w <- w / mean(w)
  }
  if(any(alarm)){warning("observations with null weight will be dropped")}
  if((nobs <- nrow(mf)) == 0){stop("zero non-NA cases")}
  if (missing(data)) data <- mf

  y <- model.response(mf)
  if (length(dim(y)) == 1L) {
    nm <- rownames(y)
    dim(y) <- NULL
    if (!is.null(nm)) names(y) <- nm
  }
  x <- if (!is.empty.model(mt)) model.matrix(mt, mf, NULL) else matrix(, nobs, 0L)
  nobs <- NROW(y)

  # Response variable
  type <- attributes(y)$type
  zyd <- cbind(y)
  if(is.null(type)){
    y <- zyd[, 1]
    z <- rep.int(-Inf, nobs)
    d <- rep.int(1, nobs)
    type <- "u"
  }
  else if (type == "right"){
    y <- zyd[, 1]
    z <- rep.int(-Inf, nobs)
    d <- zyd[, 2]
    type <- (if (any(d == 0)) "c" else "u")
  }
  else if(type == "counting"){
    z <- zyd[, 1]
    y <- zyd[, 2]
    d <- zyd[, 3]
    z <- pmax(z, min(y) - 1) # just avoid z = -Inf, that breaks findp
    type <- "ct"
    if(all(z < min(y))){type <- (if(any(d == 0)) "c" else "u")}
  }
  if(!(any(d == 1))) {stop("all data are censored")}
  attributes(type)$model <- "Qest"

  # Optional arguments for Q or wtau
  opt <- list(...)
  argwtau <- (if(!is.null(wtau)) formalArgs(wtau) else NULL)
  wtauoptions <- opt[names(opt) %in% argwtau]

  # tau and wtau. To see the reason behind omicron*1.1, see the comment to findAp
  tau <- buildTau(ntau, wtau, nobs, wtauoptions)

  if(!inherits(Q, "Qfamily")) {
    temp <- start.Qest(z, y, d, x, w, tau, Q, opt, start, data, type, control)
  }
  else {
    temp <- start.Qest.family(z, y, d, x, w, tau, wtau, wtauoptions, Q, opt, start, data, type, control)
  }
  z <- temp$z; y <- temp$y; x <- temp$x;
  tau <- temp$tau; Q <- temp$Q;
  theta <- temp$theta; eps <- temp$eps; alpha <- temp$alpha

  npar <- length(theta)

  # Fit the model #########################################################################################

  # Control arguments
  tol <- control$tol
  safeit <- control$safeit; if(is.null(safeit)){safeit <- 10 + 5*npar}
  maxit <- control$maxit; if(is.null(maxit)){maxit <- 50 + 25*npar}
  # safeit <- 25 + 5*npar
  display <- control$display
  alpha0 <- control$alpha0; if(is.null(alpha0)){alpha0 <- alpha} # For gradient step
  fit.ok <- count <- FALSE

  while(!fit.ok){
    if(count > 0 & display){cat("\n\n Non-invertible Jacobian matrix. Restarting... \n\n")}
    fit <- try(Qest.newton(theta = theta, type = type, tol = tol, maxit = maxit, safeit = safeit,
      alpha0 = alpha0, display = display, eps = eps, z = z, y = y, d = d, w = w,
      Q = Q, data = data, tau = tau), silent = FALSE)

    fit.ok <- (!inherits(fit, "try-error"))
    if(fit.ok){
      covar <- try(Qest.covar(fit, fit$eps, w), silent = FALSE)
      fit.ok <- !(inherits(covar, "try-error"))
    }

    count <- count + 1
    if(count == 5){stop("Could not fit the model (non-invertible Jacobian matrix).
        Is the model identifiable? Try with better starting points.
        Please, do not use Qest to fit discrete distributions.")}
    alpha0 <- alpha0/5
    maxit <- round(maxit * 1.25)
    safeit <- round(safeit * 1.5)
    eps <- fit$eps*1.1
    theta <- theta - fit$eps/5
  }

  # Finish ##################################################################################

  theta <- drop(fit$theta)
  L <- Loss(w, d, tau, type, fit$Fy, fit$Fz)
  attr(L, "npar") <- npar
  if(type != "u"){attr(L, "message") <- "Not the function being minimized"}

  se <- sqrt(diag(covar))

  nms <- names(theta)
  if(is.null(nms)) {
    nms <- paste0("Theta", seq_len(npar))
    names(theta) <- names(se) <- nms
  }
  rownames(covar) <- colnames(covar) <- rownames(fit$jacobian) <- colnames(fit$jacobian) <- names(fit$ee) <- nms

  if(inherits(Q, "Qfamily")){
    ok <- attr(Q, "ok")
    if(attr(Q, "glmFamily")$family %in% c("gaussian", "Gamma")) {
      ok <- c(TRUE, ok)
      names(ok)[1] <- names(theta)[1]
    }
    if(attr(Q, "glmFamily")$family %in% "uniform") {
      nms2 <- names(ok)
      ok <- rep(ok, 1 + attr(attributes(Q)$X, "min"))
      names(ok) <- if(attr(attributes(Q)$X, "min")) c(paste0("min: ", nms2), paste0("max: ", nms2)) else paste0("max: ", nms2)
    }
    se.new <- theta.new <- ok
    theta.new[ok] <- theta; theta.new[!ok] <- NA
    se.new[ok] <- se; se.new[!ok] <- NA
    theta <- theta.new
    se <- se.new
  }

  internal <- list(z = z, y = y, d = d, w = w, Q = Q, type = type, tau = tau)
  CDF <- fit$Fy$Fy; PDF <- fit$fy; names(CDF) <- names(PDF) <- rownames(mf)
  out <- list(coefficients = drop(theta), std.errs = se, covar = covar,
              obj.function = L, ee = fit$ee, jacobian = fit$jacobian, CDF = CDF, PDF = PDF,
              converged = fit$converged, n.it = fit$n.it, internal = internal,
              formula = formula, terms = mt, data = data, control = control,
              na.action = attr(mf, "na.action"), xlevels = .getXlevels(mt, mf),
              contrasts = attr(x, "contrasts"), call = cl)
  class(out) <- "Qest"
  return(out)
}







