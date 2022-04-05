
#### S3 methods for Qcoxph

print.Qcoxph <- function (x, digits = max(1L, getOption("digits") - 3L), signif.stars = FALSE, ...){
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
    cat("\n")
  }
  if (!is.null(x$fail)) {
    cat(" Coxph failed.", x$fail, "\n")
    return()
  }
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  coef <- x$coefficients
  se <- sqrt(diag(x$var))
  if (is.null(coef) | is.null(se)) stop("Input is not valid")
  tmp <- cbind(coef, exp(coef), se, coef/se, pchisq((coef/se)^2, 1, lower.tail = FALSE))
  dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", "se(coef)", "z", "p"))
  printCoefmat(tmp, digits = digits, P.values = TRUE, has.Pvalue = TRUE,
               signif.stars = signif.stars, ...)
  if (is.null(x$df)) df <- sum(!is.na(coef)) else df <- round(sum(x$df), 2)
  cat("\n")
  cat("Objective function = ", format(round(x$obj.function, 2)), " on ", df, " df\n", sep = "")
  omit <- x$na.action
  cat("n =", x$n)
  if (!is.null(x$nevent)) cat(", number of events =", x$nevent, "\n") else cat("\n")
  if (length(omit)) cat("   (", naprint(omit), ")\n", sep = "")
  invisible(x)
}

summary.Qcoxph <- function (object, conf.int = 0.95, scale = 1, ...){
  cox <- object
  beta <- cox$coefficients * scale
  if (is.null(cox$coefficients)) return(object)
  nabeta <- !(is.na(beta))
  beta2 <- beta[nabeta]
  if (is.null(beta) | is.null(cox$var)) stop("Input is not valid")
  se <- sqrt(diag(cox$var)) * scale
  rval <- list(call = cox$call, fail = cox$fail, na.action = cox$na.action,
               n = cox$n, obj.function = cox$obj.function)
  if (!is.null(cox$nevent)) rval$nevent <- cox$nevent
  tmp <- cbind(beta, exp(beta), se, beta/se, pchisq((beta/se)^2, 1, lower.tail = FALSE))
  dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)"))
  rval$coefficients <- tmp
  if (conf.int) {
    z <- qnorm((1 + conf.int)/2, 0, 1)
    tmp <- cbind(exp(beta), exp(-beta), exp(beta - z * se), exp(beta + z * se))
    dimnames(tmp) <- list(names(beta), c("exp(coef)", "exp(-coef)",
                                         paste("lower .", round(100 * conf.int, 2), sep = ""),
                                         paste("upper .", round(100 * conf.int, 2), sep = "")))
    rval$conf.int <- tmp
  }
  df <- length(beta2)
  if (!is.null(cox$concordance)) {
    rval$concordance <- cox$concordance[6:7]
    names(rval$concordance) <- c("C", "se(C)")
  }
  rval$type <- cox$internal$type
  class(rval) <- "summary.Qcoxph"
  rval
}

print.summary.Qcoxph <- function (x, digits = max(getOption("digits") - 3, 3), signif.stars = getOption("show.signif.stars"), ...){
  if (!is.null(x$call)) {
    cat("Call:\n")
    dput(x$call)
    cat("\n")
  }
  if (!is.null(x$fail)) {
    cat(" Coxreg failed.", x$fail, "\n")
    return()
  }
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  omit <- x$na.action
  cat("  n =", x$n)
  if (!is.null(x$nevent)) cat(", number of events =", x$nevent, "\n") else cat("\n")
  if (length(omit)) cat("   (", naprint(omit), ")\n", sep = "")
  if (nrow(x$coef) == 0) {
    cat("   Null model\n")
    return()
  }
  if (!is.null(x$coefficients)) {
    cat("\n")
    printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, ...)
  }
  if (!is.null(x$conf.int)) {
    cat("\n")
    print(x$conf.int)
  }
  cat("\n")
  if (!is.null(x$concordance)) cat("Concordance = ", format(round(x$concordance[1], 3)), " (se = ", format(round(x$concordance[2], 3)), ")", sep = "")
  cat("\nObjective function = ", format(round(x$obj.function, 2), digits = digits), sep = "")
  if(x$type != "u"){cat(" (not the function being minimized)\n")}
  cat("\n")
  invisible()
}

basehaz.Qcoxph <- function(fit, centered = TRUE, se.fit = FALSE){
  deltamethod <- function(vgamma, gamma){
    n <- length(gamma)
    out <- matrix(NA, n, n)
    for(i in 1:n)
      for(j in i:n)
        out[i, j] <- out[j, i] <- exp(gamma)[i] %*% vgamma[i,j] %*% exp(gamma)[j]
    out
  }
  coxph.getdata <- getFromNamespace("coxph.getdata", "survival")
  mf <- coxph.getdata(fit)
  x <- mf$x; y <- mf$y
  # tab <- data.frame(table(y[y[, 2] == 1, 1]))
  # t <- as.numeric(levels(tab[, 1]))[tab[, 1]]
  t <- sort(unique(y[,1]))
  knots <- fit$internal$knots
  gamma <- fit$internal$gamma
  by <- plfcox(t, knots, deriv = 0)
  h0 <- c(by %*% exp(gamma))
  if (centered) {
    zcoef <- ifelse(is.na(coef(fit)), 0, coef(fit))
    offset <- sum(fit$means * zcoef)
    h0 <- h0 * exp(offset)
  }
  out <- data.frame(hazard = h0, time = t)
  if(se.fit) {
    vgamma <- deltamethod(fit$internal$var_gamma, gamma)
    V <- by %*% vgamma %*% t(by)
    varhaz <- diag(V)
    out$varhaz <- varhaz
  }
  out
}

agsurv.Qcoxph <- function (y, x, wt, risk, fit) {
  coxph.getdata <- getFromNamespace("coxph.getdata", "survival")
  mf <- coxph.getdata(fit)
  x <- mf$x; y <- mf$y; wt <- mf$weights
  risk <- exp(fit$linear.predictors)

  nvar <- ncol(as.matrix(x))
  status <- y[, ncol(y)]
  dtime <- y[, ncol(y) - 1]
  death <- (status == 1)
  time <- sort(unique(dtime))
  nevent <- as.vector(rowsum(wt * death, dtime))
  ncens <- as.vector(rowsum(wt * (!death), dtime))
  wrisk <- wt * risk
  rcumsum <- function(x) rev(cumsum(rev(x)))
  nrisk <- rcumsum(rowsum(wrisk, dtime))
  irisk <- rcumsum(rowsum(wt, dtime))
  if (ncol(y) == 2) {
    temp2 <- rowsum(wrisk * x, dtime)
    xsum <- apply(temp2, 2, rcumsum)
  }
  else {
    delta <- min(diff(time))/2
    etime <- c(sort(unique(y[, 1])), max(y[, 1]) + delta)
    indx <- approx(etime, 1:length(etime), time, method = "constant",
                   rule = 2, f = 1)$y
    esum <- rcumsum(rowsum(wrisk, y[, 1]))
    nrisk <- nrisk - c(esum, 0)[indx]
    irisk <- irisk - c(rcumsum(rowsum(wt, y[, 1])), 0)[indx]
    xout <- apply(rowsum(wrisk * x, y[, 1]), 2, rcumsum)
    xin <- apply(rowsum(wrisk * x, dtime), 2, rcumsum)
    xsum <- xin - (rbind(xout, 0))[indx, , drop = F]
  }
  ndeath <- rowsum(status, dtime)
  ntime <- length(time)

  # Cagsurv5 <- getFromNamespace("Cagsurv5", "survival")
  xsum2 <- rowsum((wrisk * death) * x, dtime)
  erisk <- rowsum(wrisk * death, dtime)
  tsum <- .C("Cagsurv5", as.integer(length(nevent)), as.integer(nvar),
             as.integer(ndeath), as.double(nrisk), as.double(erisk),
             as.double(xsum), as.double(xsum2), sum1 = double(length(nevent)),
             sum2 = double(length(nevent)), xbar = matrix(0, length(nevent), nvar))

  tsum2 <- basehaz.Qcoxph(fit, se.fit=TRUE)
  haz <- tsum2$hazard
  varhaz <- tsum2$varhaz
  xbar <- nevent * tsum$xbar
  result <- list(n = nrow(y), time = time, n.event = nevent,
                 n.risk = irisk, n.censor = ncens, hazard = diff(c(0, haz)), cumhaz = haz,
                 varhaz = varhaz, ndeath = ndeath, xbar = apply(matrix(xbar, ncol = nvar), 2, cumsum))
  result
}

coxsurv.fit.Qcoxph <- function (ctype, stype, se.fit, varmat, cluster, y, x, wt, risk,
                                position, strata, oldid, y2, x2, risk2, strata2, id2,
                                unlist = TRUE, fit) {
  if (missing(strata) || length(strata) == 0) strata <- rep(0L, nrow(y))
  if (is.factor(strata))
    ustrata <- levels(strata)
  else ustrata <- sort(unique(strata))
  nstrata <- length(ustrata)
  survlist <- vector("list", nstrata)
  names(survlist) <- ustrata
  survtype <- if (stype == 1)
    1
  else ctype + 1
  vartype <- survtype
  for (i in 1:nstrata) {
    indx <- which(strata == ustrata[i])
    survlist[[i]] <- agsurv.Qcoxph(y[indx, , drop = F], x[indx, , drop = F],
                                   wt[indx], risk[indx], fit)
  }
  expand <- function(fit, x2, varmat, se.fit) {
    # if (survtype == 1)
    #   surv <- cumprod(fit$surv)
    # else surv <- exp(-fit$cumhaz)
    surv <- exp(-fit$cumhaz)
    if (is.matrix(x2) && nrow(x2) > 1) {
      fit$surv <- outer(surv, risk2, "^")
      dimnames(fit$surv) <- list(NULL, row.names(x2))
      if (se.fit) {
        varh <- matrix(0, nrow = length(fit$varhaz), ncol = nrow(x2))
        for (i in 1:nrow(x2)) {
          dt <- outer(fit$cumhaz, x2[i, ], "*") - fit$xbar
          varh[, i] <- (cumsum(fit$varhaz) + rowSums((dt %*% varmat) * dt)) * risk2[i]^2
        }
        fit$std.err <- sqrt(varh)
      }
      fit$cumhaz <- outer(fit$cumhaz, risk2, "*")
    }
    else {
      fit$surv <- surv^risk2
      if (se.fit) {
        dt <- outer(fit$cumhaz, c(x2)) - fit$xbar
        varh <- (cumsum(fit$varhaz) + rowSums((dt %*% varmat) * dt)) * risk2^2
        fit$std.err <- sqrt(varh)
      }
      fit$cumhaz <- fit$cumhaz * risk2
    }
    fit
  }
  if (missing(id2) || is.null(id2))
    result <- lapply(survlist, expand, x2, varmat, se.fit)
  else {
    onecurve <- function(slist, x2, y2, strata2, risk2, se.fit) {
      ntarget <- nrow(x2)
      surv <- vector("list", ntarget)
      n.event <- n.risk <- n.censor <- varh1 <- varh2 <- time <- surv
      hazard <- vector("list", ntarget)
      stemp <- as.integer(strata2)
      timeforward <- 0
      for (i in 1:ntarget) {
        slist <- survlist[[stemp[i]]]
        indx <- which(slist$time > y2[i, 1] & slist$time <= y2[i, 2])
        if (length(indx) == 0) {
          timeforward <- timeforward + y2[i, 2] - y2[i, 1]
        }
        else {
          time[[i]] <- diff(c(y2[i, 1], slist$time[indx]))
          time[[i]][1] <- time[[i]][1] + timeforward
          timeforward <- y2[i, 2] - max(slist$time[indx])
          hazard[[i]] <- slist$hazard[indx] * risk2[i]
          if (survtype == 1) surv[[i]] <- slist$surv[indx]^risk2[i]
          n.event[[i]] <- slist$n.event[indx]
          n.risk[[i]] <- slist$n.risk[indx]
          n.censor[[i]] <- slist$n.censor[indx]
          dt <- outer(slist$cumhaz[indx], x2[i, ]) - slist$xbar[indx, , drop = F]
          varh1[[i]] <- slist$varhaz[indx] * risk2[i]^2
          varh2[[i]] <- rowSums((dt %*% varmat) * dt) * risk2[i]^2
        }
      }
      cumhaz <- cumsum(unlist(hazard))
      if (survtype == 1)
        surv <- cumprod(unlist(surv))
      else surv <- exp(-cumhaz)
      if (se.fit)
        list(n = as.vector(table(strata)[stemp[1]]),
             time = cumsum(unlist(time)), n.risk = unlist(n.risk),
             n.event = unlist(n.event), n.censor = unlist(n.censor),
             surv = surv, cumhaz = cumhaz, std.err = sqrt(cumsum(unlist(varh1)) + unlist(varh2)))
      else list(n = as.vector(table(strata)[stemp[1]]),
                time = cumsum(unlist(time)), n.risk = unlist(n.risk),
                n.event = unlist(n.event), n.censor = unlist(n.censor),
                surv = surv, cumhaz = cumhaz)
    }
    if (all(id2 == id2[1])) {
      result <- list(onecurve(survlist, x2, y2, strata2, risk2, se.fit))
    }
    else {
      uid <- unique(id2)
      result <- vector("list", length = length(uid))
      for (i in 1:length(uid)) {
        indx <- which(id2 == uid[i])
        result[[i]] <- onecurve(survlist, x2[indx, ,drop = FALSE], y2[indx, , drop = FALSE],
                                strata2[indx], risk2[indx], se.fit)
      }
      names(result) <- uid
    }
  }
  if (unlist) {
    if (length(result) == 1) {
      if (se.fit)
        result[[1]][c("n", "time", "n.risk", "n.event",
                      "n.censor", "surv", "cumhaz", "std.err")]
      else result[[1]][c("n", "time", "n.risk", "n.event",
                         "n.censor", "surv", "cumhaz")]
    }
    else {
      temp <- list(n = unlist(lapply(result, function(x) x$n), use.names = FALSE),
                   time = unlist(lapply(result, function(x) x$time), use.names = FALSE),
                   n.risk = unlist(lapply(result, function(x) x$n.risk), use.names = FALSE),
                   n.event = unlist(lapply(result, function(x) x$n.event), use.names = FALSE),
                   n.censor = unlist(lapply(result, function(x) x$n.censor), use.names = FALSE),
                   strata = sapply(result, function(x) length(x$time)))
      names(temp$strata) <- names(result)
      if ((missing(id2) || is.null(id2)) && nrow(x2) > 1) {
        temp$surv <- t(matrix(unlist(lapply(result, function(x) t(x$surv)), use.names = FALSE), nrow = nrow(x2)))
        dimnames(temp$surv) <- list(NULL, row.names(x2))
        temp$cumhaz <- t(matrix(unlist(lapply(result, function(x) t(x$cumhaz)), use.names = FALSE), nrow = nrow(x2)))
        if (se.fit)
          temp$std.err <- t(matrix(unlist(lapply(result, function(x) t(x$std.err)), use.names = FALSE), nrow = nrow(x2)))
      }
      else {
        temp$surv <- unlist(lapply(result, function(x) x$surv), use.names = FALSE)
        temp$cumhaz <- unlist(lapply(result, function(x) x$cumhaz), use.names = FALSE)
        if (se.fit)
          temp$std.err <- unlist(lapply(result, function(x) x$std.err), use.names = FALSE)
      }
      temp
    }
  }
  else {
    names(result) <- ustrata
    result
  }
}

survfit.Qcoxph <- function (formula, newdata, se.fit = TRUE, conf.int = 0.95, individual = FALSE,
                            stype = 2, ctype, conf.type = c("log", "log-log", "plain","none", "logit", "arcsin"),
                            censor = TRUE, start.time, id, influence = FALSE, na.action = na.pass, type, ...) {
  Call <- match.call()
  Call[[1]] <- as.name("survfit")
  object <- formula
  Terms <- terms(object)
  robust <- !is.null(object$naive.var)
  if (!is.null(attr(object$terms, "specials")$tt))
    stop("The survfit function can not process coxph models with a tt term")
  if (!missing(type)) {
    if (!missing(stype) || !missing(ctype))
      warning("type argument ignored")
    else {
      temp1 <- c("kalbfleisch-prentice", "aalen", "efron",
                 "kaplan-meier", "breslow", "fleming-harrington",
                 "greenwood", "tsiatis", "exact")
      survtype <- match(match.arg(type, temp1), temp1)
      stype <- c(1, 2, 2, 1, 2, 2, 2, 2, 2)[survtype]
      if (stype != 1) ctype <- c(1, 1, 2, 1, 1, 2, 1, 1, 1)[survtype]
    }
  }
  if (missing(ctype)) {
    temp1 <- match(object$method, c("exact", "breslow", "efron"))
    ctype <- c(1, 1, 2)[temp1]
  }
  else if (!(ctype %in% 1:2))
    stop("ctype must be 1 or 2")
  if (!se.fit)
    conf.type <- "none"
  else conf.type <- match.arg(conf.type)
  tfac <- attr(Terms, "factors")
  temp <- attr(Terms, "specials")$strata
  has.strata <- !is.null(temp)
  if (has.strata) {
    stangle = untangle.specials(Terms, "strata")
    for (i in temp) tfac <- tfac[, tfac[i, ] == 0]
  }
  if (any(tfac > 1))
    stop("not able to create a curve for models that contain an interaction without the lower order effect")
  Terms <- object$terms
  n <- object$n[1]
  if (!has.strata)
    strata <- rep(0L, n)
  else strata <- object$strata
  missid <- missing(id)
  if (!missid & !missing(individual)) warning("the `id' option supersedes `individual'")
  if (!missid)
    individual <- TRUE
  else if (missid && individual)
    id <- rep(0, n)
  else id <- NULL
  if (individual & missing(newdata)) {
    stop("the id and/or individual options only make sense with new data")
  }
  if (has.strata) {
    temp <- attr(Terms, "specials")$strata
    factors <- attr(Terms, "factors")[temp, ]
    strata.interaction <- any(t(factors) * attr(Terms, "order") > 1)
  }
  coxms <- inherits(object, "coxphms")
  if (coxms || is.null(object$y) || is.null(object[["x"]]) ||
      !is.null(object$call$weights) || !is.null(object$call$id) ||
      (has.strata && is.null(object$strata)) || !is.null(attr(object$terms, "offset"))) {
    mf <- stats::model.frame(object)
  }
  else mf <- NULL
  position <- NULL
  Y <- object[["y"]]
  if (is.null(mf)) {
    weights <- rep(1, n)
    offset <- rep(0, n)
    X <- object[["x"]]
  }
  else {
    model.matrix.coxph <- getFromNamespace("model.matrix.coxph", "survival")
    weights <- model.weights(mf)
    if (is.null(weights)) weights <- rep(1, n)
    offset <- model.offset(mf)
    if (is.null(offset)) offset <- rep(0, n)
    X <- model.matrix.coxph(object, data = mf)
    if (is.null(Y) || coxms) {
      Y <- model.response(mf)
      if (is.null(object$timefix) || object$timefix) Y <- aeqSurv(Y)
    }
    oldid <- model.extract(mf, "id")
    survflag <- getFromNamespace("survflag", "survival")
    if (length(oldid) && ncol(Y) == 3)
      position <- survflag(Y, oldid)
    else position <- NULL
    if (!coxms && (nrow(Y) != object$n[1]))
      stop("Failed to reconstruct the original data set")
    if (has.strata) {
      if (length(strata) == 0) {
        if (length(stangle$vars) == 1)
          strata <- mf[[stangle$vars]]
        else strata <- strata(mf[, stangle$vars], shortlabel = TRUE)
      }
    }
  }
  if (FALSE) {
    if (!is.null(mf)) {
      y2 <- object[["y"]]
      if (!is.null(y2)) {
        if (ncol(y2) != ncol(Y) || length(y2) != length(Y))
          stop("Could not reconstruct the y vector")
      }
    }
  }
  type <- attr(Y, "type")
  if (!type %in% c("right", "counting", "mright", "mcounting"))
    stop("Cannot handle \"", type, "\" type survival data")
  if (type == "right" || type == "mright")
    individual <- FALSE
  if (!missing(start.time)) {
    if (!is.numeric(start.time) || length(start.time) > 1) stop("start.time must be a single numeric value")
    if (ncol(Y) == 3) {
      keep <- Y[, 2] > start.time
      Y[keep, 1] <- pmax(Y[keep, 1], start.time)
    }
    else keep <- Y[, 1] > start.time
    if (!any(Y[keep, ncol(Y)] == 1)) stop("start.time argument has removed all endpoints")
    Y <- Y[keep, , drop = FALSE]
    X <- X[keep, , drop = FALSE]
    offset <- offset[keep]
    strata <- strata[keep]
    if (length(id) > 0) id <- id[keep]
    if (length(position)) position <- position[keep]
    n <- nrow(Y)
  }
  if (length(object$means) == 0) {
    X <- matrix(0, nrow = n, ncol = 1)
    xcenter <- mean(offset)
    coef <- 0
    varmat <- matrix(0, 1, 1)
    risk <- rep(exp(offset - mean(offset)), length = n)
  }
  else {
    varmat <- object$var
    beta <- ifelse(is.na(object$coefficients), 0, object$coefficients)
    xcenter <- sum(object$means * beta) + mean(offset)
    if (!is.null(object$frail)) {
      keep <- !grepl("frailty(", dimnames(X)[[2]], fixed = TRUE)
      X <- X[, keep, drop = F]
    }
    risk <- c(exp(X %*% beta + offset - xcenter))
  }
  if (missing(newdata)) {
    if (any(attr(Terms, "order") > 1))
      warning("the model contains interactions; the default curve based on columm means of the X matrix is almost certainly not useful. Consider adding a newdata argument.")
    if (length(object$means)) {
      mf2 <- as.list(object$means)
      names(mf2) <- names(object$coefficients)
      mf2 <- as.data.frame(mf2)
      x2 <- matrix(object$means, 1)
    }
    else {
      mf2 <- data.frame(X = 0)
      x2 <- 0
    }
    offset2 <- 0
    found.strata <- FALSE
  }
  else {
    if (!is.null(object$frail)) stop("Newdata cannot be used when a model has frailty terms")
    Terms2 <- Terms
    if (!individual) Terms2 <- delete.response(Terms)
    if (is.vector(newdata, "numeric")) {
      if (individual) stop("newdata must be a data frame")
      if (is.null(names(newdata))) {
        stop("Newdata argument must be a data frame")
      }
      newdata <- data.frame(as.list(newdata), stringsAsFactors = FALSE)
    }
    if (has.strata) {
      found.strata <- TRUE
      tempenv <- new.env(, parent = emptyenv())
      assign("strata", function(..., na.group, shortlabel, sep) list(...), envir = tempenv)
      assign("list", list, envir = tempenv)
      for (svar in stangle$vars) {
        temp <- try(eval(parse(text = svar), newdata, tempenv), silent = TRUE)
        if (!is.list(temp) || any(unlist(lapply(temp, class)) == "function")) found.strata <- FALSE
      }
      if (!found.strata) {
        ss <- untangle.specials(Terms2, "strata")
        Terms2 <- Terms2[-ss$terms]
      }
    }
    tcall <- Call[c(1, match(c("id", "na.action"), names(Call), nomatch = 0))]
    tcall$data <- newdata
    tcall$formula <- Terms2
    tcall$xlev <- object$xlevels[match(attr(Terms2, "term.labels"), names(object$xlevels), nomatch = 0)]
    tcall[[1L]] <- quote(stats::model.frame)
    mf2 <- eval(tcall)
  }
  if (has.strata && found.strata) {
    temp <- untangle.specials(Terms2, "strata")
    strata2 <- strata(mf2[temp$vars], shortlabel = TRUE)
    strata2 <- factor(strata2, levels = levels(strata))
    if (any(is.na(strata2))) stop("New data set has strata levels not found in the original")
    if (length(temp$terms) > 0) Terms2 <- Terms2[-temp$terms]
  }
  else strata2 <- factor(rep(0, nrow(mf2)))
  if (!robust) cluster <- NULL
  if (individual) {
    if (missing(newdata)) stop("The newdata argument must be present when individual=TRUE")
    if (!missid) {
      id2 <- model.extract(mf2, "id")
      if (is.null(id2)) stop("id=NULL is an invalid argument")
    }
    else id2 <- rep(1, nrow(mf2))
    x2 <- model.matrix(Terms2, mf2)[, -1, drop = FALSE]
    if (length(x2) == 0) stop("Individual survival but no variables")
    offset2 <- model.offset(mf2)
    if (length(offset2) == 0) offset2 <- 0
    y2 <- model.extract(mf2, "response")
    if (attr(y2, "type") != type) stop("Survival type of newdata does not match the fitted model")
    if (attr(y2, "type") != "counting") stop("Individual=TRUE is only valid for counting process data")
    y2 <- y2[, 1:2, drop = F]
  }
  else if (missing(newdata)) {
    if (has.strata && strata.interaction)
      stop("Models with strata by covariate interaction terms require newdata")
    offset2 <- 0
    if (length(object$means)) {
      x2 <- matrix(object$means, nrow = 1, ncol = ncol(X))
    }
    else {
      x2 <- matrix(0, nrow = 1, ncol = 1)
    }
  }
  else {
    offset2 <- model.offset(mf2)
    if (length(offset2) > 0)
      offset2 <- offset2
    else offset2 <- 0
    x2 <- model.matrix(Terms2, mf2)[, -1, drop = FALSE]
  }
  if (missing(newdata))
    risk2 <- 1
  else {
    if (length(object$means))
      risk2 <- exp(c(x2 %*% beta) + offset2 - xcenter)
    else risk2 <- exp(offset2 - xcenter)
  }
  if (individual) {
    result <- coxsurv.fit.Qcoxph(ctype, stype, se.fit, varmat, cluster,
                                 Y, X, weights, risk, position, strata, oldid, y2,
                                 x2, risk2, strata2, id2, fit=object)
  }
  else {
    result <- coxsurv.fit.Qcoxph(ctype, stype, se.fit, varmat, cluster,
                                 Y, X, weights, risk, position, strata, oldid, y2,
                                 x2, risk2, fit=object)
    if (has.strata && found.strata) {
      if (is.matrix(result$surv)) {
        nr <- nrow(result$surv)
        indx1 <- split(1:nr, rep(1:length(result$strata), result$strata))
        rows <- indx1[as.numeric(strata2)]
        indx2 <- unlist(rows)
        indx3 <- as.integer(strata2)
        for (i in 2:length(rows)) rows[[i]] <- rows[[i]] + (i - 1) * nr
        indx4 <- unlist(rows)
        temp <- result$strata[indx3]
        names(temp) <- row.names(mf2)
        new <- list(n = result$n[indx3], time = result$time[indx2],
                    n.risk = result$n.risk[indx2], n.event = result$n.event[indx2],
                    n.censor = result$n.censor[indx2], strata = temp,
                    surv = result$surv[indx4], cumhaz = result$cumhaz[indx4])
        if (se.fit) new$std.err <- result$std.err[indx4]
        result <- new
      }
    }
  }
  if (!censor) {
    kfun <- function(x, keep) {
      if (is.matrix(x))
        x[keep, , drop = F]
      else if (length(x) == length(keep))
        x[keep]
      else x
    }
    keep <- (result$n.event > 0)
    if (!is.null(result$strata)) {
      temp <- factor(rep(names(result$strata), result$strata), levels = names(result$strata))
      result$strata <- c(table(temp[keep]))
    }
    result <- lapply(result, kfun, keep)
  }
  result$logse = TRUE
  if (se.fit && conf.type != "none") {
    survfit_confint <- getFromNamespace("survfit_confint", "survival")
    ci <- survfit_confint(result$surv, result$std.err, logse = result$logse, conf.type, conf.int)
    result <- c(result, list(lower = ci$lower, upper = ci$upper, conf.type = conf.type, conf.int = conf.int))
  }
  if (!missing(start.time)) result$start.time <- start.time
  result$call <- Call
  class(result) <- c("survfitQcox", "survfitcox", "survfit")
  result
}

residuals.Qcoxph <- function (object, type = c("martingale", "deviance", "score",
                                               "schoenfeld", "dfbeta", "dfbetas", "scaledsch", "partial"),
                              collapse = FALSE, weighted = FALSE, ...) {
  type <- match.arg(type)
  otype <- type
  if(type %in% c("score", "schoenfeld", "dfbeta", "dfbetas", "scaledsch"))
    stop("type not yet implemented")
  # if (type == "dfbeta" || type == "dfbetas") {
  #   otype <- type
  #   type <- "score"
  #   if (missing(weighted)) weighted <- TRUE
  # }
  # if (type == "scaledsch") type <- "schoenfeld"
  n <- length(object$residuals)
  rr <- object$residuals
  y <- object$y
  x <- object[["x"]]
  vv <- drop(object$naive.var)
  if (is.null(vv)) vv <- drop(object$var)
  weights <- object$weights
  if (is.null(weights)) weights <- rep(1, n)
  strat <- object$strata
  # method <- object$method
  # if (method == "exact" && (type == "score" || type == "schoenfeld"))
  #   stop(paste(otype, "residuals are not available for the exact method"))
  if (type == "martingale" || type == "partial")
    rr <- object$residuals
  else {
    Terms <- object$terms
    if (!inherits(Terms, "terms")) stop("invalid terms component of object")
    strats <- attr(Terms, "specials")$strata
    if (is.null(y) || (is.null(x) && type != "deviance")) {
      coxph.getdata <- getFromNamespace("coxph.getdata", "survival")
      temp <- coxph.getdata(object, y = TRUE, x = TRUE, stratax = TRUE)
      y <- temp$y
      x <- temp$x
      if (length(strats)) strat <- temp$strata
    }
    ny <- ncol(y)
    status <- y[, ny, drop = TRUE]
    if (type != "deviance") {
      nstrat <- as.numeric(strat)
      nvar <- ncol(x)
      if (is.null(strat)) {
        ord <- order(y[, ny - 1], -status)
        newstrat <- rep(0, n)
      }
      else {
        ord <- order(nstrat, y[, ny - 1], -status)
        newstrat <- c(diff(as.numeric(nstrat[ord])) != 0, 1)
      }
      newstrat[n] <- 1
      x <- x[ord, ]
      y <- y[ord, ]
      score <- exp(object$linear.predictors)[ord]
    }
  }
  # if (type == "schoenfeld") {
  #   if (ny == 2) {
  #     mintime <- min(y[, 1])
  #     if (mintime < 0)
  #       y <- cbind(2 * mintime - 1, y)
  #     else y <- cbind(-1, y)
  #   }
  #   temp <- .C(Ccoxscho, n = as.integer(n), as.integer(nvar), as.double(y), resid = as.double(x),
  #              as.double(score * weights[ord]), as.integer(newstrat), as.integer(method == "efron"),
  #              double(3 * nvar))
  #   deaths <- y[, 3] == 1
  #   if (nvar == 1)
  #     rr <- temp$resid[deaths]
  #   else rr <- matrix(temp$resid[deaths], ncol = nvar)
  #   if (weighted) rr <- rr * weights[deaths]
  #   if (length(strats)) attr(rr, "strata") <- table((strat[ord])[deaths])
  #   time <- c(y[deaths, 2])
  #   if (is.matrix(rr))
  #     dimnames(rr) <- list(time, names(object$coefficients))
  #   else names(rr) <- time
  #   if (otype == "scaledsch") {
  #     ndead <- sum(deaths)
  #     coef <- ifelse(is.na(object$coefficients), 0, object$coefficients)
  #     if (nvar == 1)
  #       rr <- rr * vv * ndead + coef
  #     else rr <- drop(rr %*% vv * ndead + rep(coef, each = nrow(rr)))
  #   }
  #   return(rr)
  # }
  # if (type == "score") {
  #   storage.mode(y) <- storage.mode(x) <- "double"
  #   storage.mode(newstrat) <- "integer"
  #   storage.mode(score) <- storage.mode(weights) <- "double"
  #   if (ny == 2) {
  #     resid <- .Call(Ccoxscore2, y, x, newstrat, score, weights[ord], as.integer(method == "efron"))
  #   }
  #   else {
  #     resid <- .Call(Cagscore2, y, x, newstrat, score, weights[ord], as.integer(method == "efron"))
  #   }
  #   if (nvar > 1) {
  #     rr <- matrix(0, n, nvar)
  #     rr[ord, ] <- resid
  #     dimnames(rr) <- list(names(object$residuals), names(object$coefficients))
  #   }
  #   else rr[ord] <- resid
  #   if (otype == "dfbeta") {
  #     if (is.matrix(rr))
  #       rr <- rr %*% vv
  #     else rr <- rr * vv
  #   }
  #   else if (otype == "dfbetas") {
  #     if (is.matrix(rr))
  #       rr <- (rr %*% vv) %*% diag(sqrt(1/diag(vv)))
  #     else rr <- rr * sqrt(vv)
  #   }
  # }
  if (weighted)
    rr <- rr * weights
  if (!is.null(object$na.action)) {
    rr <- naresid(object$na.action, rr)
    if (is.matrix(rr))
      n <- nrow(rr)
    else n <- length(rr)
    if (type == "deviance")
      status <- naresid(object$na.action, status)
  }
  if (type == "partial") {
    rr <- rr + predict(object, type = "terms")
  }
  if (!missing(collapse)) {
    if (length(collapse) != n) stop("Wrong length for 'collapse'")
    rr <- drop(rowsum(rr, collapse))
    if (type == "deviance") status <- drop(rowsum(status, collapse))
  }
  if (type == "deviance")
    sign(rr) * sqrt(-2 * (rr + ifelse(status == 0, 0, status * log(status - rr))))
  else rr
}

predict.Qcoxph <- function (object, newdata, type = c("lp", "risk", "expected", "terms", "survival"),
                            se.fit = FALSE, na.action = na.pass,
                            terms = names(object$assign), collapse, reference = c("strata", "sample"), ...) {
  if (!inherits(object, "Qcoxph")) stop("Primary argument much be a Qcoxph object")
  Call <- match.call()
  type <- match.arg(type)
  if (type == "survival") {
    survival <- TRUE
    type <- "expected"
  }
  else survival <- FALSE
  n <- object$n
  Terms <- object$terms
  if (!missing(terms)) {
    if (is.numeric(terms)) {
      if (any(terms != floor(terms) | terms > length(object$assign) | terms < 1))
        stop("Invalid terms argument")
    }
    else if (any(is.na(match(terms, names(object$assign)))))
      stop("a name given in the terms argument not found in the model")
  }
  if (length(attr(Terms, "specials")$cluster)) {
    temp <- untangle.specials(Terms, "cluster", 1)
    Terms <- object$terms[-temp$terms]
  }
  else Terms <- object$terms
  if (type != "expected")
    Terms2 <- delete.response(Terms)
  else Terms2 <- Terms
  has.strata <- !is.null(attr(Terms, "specials")$strata)
  has.offset <- !is.null(attr(Terms, "offset"))
  has.weights <- any(names(object$call) == "weights")
  na.action.used <- object$na.action
  n <- length(object$residuals)
  if (missing(reference) && type == "terms")
    reference <- "sample"
  else reference <- match.arg(reference)
  have.mf <- FALSE
  if (type == "expected") {
    y <- object[["y"]]
    if (is.null(y)) {
      mf <- stats::model.frame(object)
      y <- model.extract(mf, "response")
      have.mf <- TRUE
    }
  }
  strat.term <- untangle.specials(Terms, "strata")
  if (se.fit || type == "terms" || (!missing(newdata) && type == "expected") ||
      (has.strata && (reference == "strata") || type == "expected")) {
    use.x <- TRUE
    if (is.null(object[["x"]]) || has.weights || has.offset || (has.strata && is.null(object$strata))) {
      if (!have.mf) mf <- stats::model.frame(object)
      if (nrow(mf) != n) stop("Data is not the same size as it was in the original fit")
      x <- model.matrix(object, data = mf)
      if (has.strata) {
        if (!is.null(object$strata))
          oldstrat <- object$strata
        else {
          if (length(strat.term$vars) == 1)
            oldstrat <- mf[[strat.term$vars]]
          else oldstrat <- strata(mf[, strat.term$vars], shortlabel = TRUE)
        }
      }
      else oldstrat <- rep(0L, n)
      weights <- model.weights(mf)
      if (is.null(weights)) weights <- rep(1, n)
      offset <- model.offset(mf)
      if (is.null(offset)) offset <- 0
    }
    else {
      x <- object[["x"]]
      if (has.strata)
        oldstrat <- object$strata
      else oldstrat <- rep(0L, n)
      weights <- rep(1, n)
      offset <- 0
    }
  }
  else {
    if (has.strata) {
      stemp <- untangle.specials(Terms, "strata", 1)
      Terms2 <- Terms2[-stemp$terms]
      has.strata <- FALSE
    }
    oldstrat <- rep(0L, n)
    offset <- 0
    use.x <- FALSE
  }
  if (!missing(newdata)) {
    use.x <- TRUE
    tcall <- Call[c(1, match(c("newdata", "collapse"), names(Call), nomatch = 0))]
    names(tcall)[2] <- "data"
    tcall$formula <- Terms2
    tcall$na.action <- na.action
    tcall[[1L]] <- quote(stats::model.frame)
    if (!is.null(attr(Terms, "specials")$strata) && !has.strata) {
      temp.lev <- object$xlevels
      temp.lev[[strat.term$vars]] <- NULL
      tcall$xlev <- temp.lev
    }
    else tcall$xlev <- object$xlevels
    mf2 <- eval(tcall, parent.frame())
    collapse <- model.extract(mf2, "collapse")
    n2 <- nrow(mf2)
    if (has.strata) {
      if (length(strat.term$vars) == 1)
        newstrat <- mf2[[strat.term$vars]]
      else newstrat <- strata(mf2[, strat.term$vars], shortlabel = TRUE)
      if (any(is.na(match(newstrat, oldstrat))))
        stop("New data has a strata not found in the original model")
      else newstrat <- factor(newstrat, levels = levels(oldstrat))
      if (length(strat.term$terms))
        newx <- model.matrix(Terms2[-strat.term$terms], mf2, contr = object$contrasts)[, -1, drop = FALSE]
      else newx <- model.matrix(Terms2, mf2, contr = object$contrasts)[, -1, drop = FALSE]
    }
    else {
      newx <- model.matrix(Terms2, mf2, contr = object$contrasts)[, -1, drop = FALSE]
      newstrat <- rep(0L, nrow(mf2))
    }
    newoffset <- model.offset(mf2)
    if (is.null(newoffset)) newoffset <- 0
    if (type == "expected") {
      newy <- model.response(mf2)
      if (attr(newy, "type") != attr(y, "type")) stop("New data has a different survival type than the model")
    }
    na.action.used <- attr(mf2, "na.action")
  }
  else n2 <- n
  if (type == "expected") {
    if (missing(newdata)) pred <- y[, ncol(y)] - object$residuals
    if (!missing(newdata) || se.fit) {
      ustrata <- unique(oldstrat)
      risk <- exp(object$linear.predictors)
      x <- x - rep(object$means, each = nrow(x))
      if (missing(newdata))
        se <- double(n)
      else {
        pred <- se <- double(nrow(mf2))
        newx <- newx - rep(object$means, each = nrow(newx))
        newrisk <- c(exp(newx %*% object$coef) + newoffset)
      }
      survtype <- ifelse(object$method == "efron", 3, 2)
      for (i in ustrata) {
        indx <- which(oldstrat == i)
        afit <- agsurv.Qcoxph(y[indx, , drop = F], x[indx, , drop = F], weights[indx], risk[indx], object)
        afit.n <- length(afit$time)
        if (missing(newdata)) {
          j1 <- approx(afit$time, 1:afit.n, y[indx, 1], method = "constant", f = 0, yleft = 0, yright = afit.n)$y
          chaz <- c(0, afit$cumhaz)[j1 + 1]
          varh <- c(0, cumsum(afit$varhaz))[j1 + 1]
          xbar <- rbind(0, afit$xbar)[j1 + 1, , drop = F]
          if (ncol(y) == 2) {
            dt <- (chaz * x[indx, ]) - xbar
            se[indx] <- sqrt(varh + rowSums((dt %*% object$var) * dt)) * risk[indx]
          }
          else {
            j2 <- approx(afit$time, 1:afit.n, y[indx, 2], method = "constant", f = 0, yleft = 0, yright = afit.n)$y
            chaz2 <- c(0, afit$cumhaz)[j2 + 1]
            varh2 <- c(0, cumsum(afit$varhaz))[j2 + 1]
            xbar2 <- rbind(0, afit$xbar)[j2 + 1, , drop = F]
            dt <- (chaz * x[indx, ]) - xbar
            v1 <- varh + rowSums((dt %*% object$var) * dt)
            dt2 <- (chaz2 * x[indx, ]) - xbar2
            v2 <- varh2 + rowSums((dt2 %*% object$var) * dt2)
            se[indx] <- sqrt(v2 - v1) * risk[indx]
          }
        }
        else {
          use.x <- TRUE
          indx2 <- which(newstrat == i)
          j1 <- approx(afit$time, 1:afit.n, newy[indx2, 1], method = "constant", f = 0, yleft = 0, yright = afit.n)$y
          chaz <- c(0, afit$cumhaz)[j1 + 1]
          pred[indx2] <- chaz * newrisk[indx2]
          if (se.fit) {
            varh <- c(0, cumsum(afit$varhaz))[j1 + 1]
            xbar <- rbind(0, afit$xbar)[j1 + 1, , drop = F]
          }
          if (ncol(y) == 2) {
            if (se.fit) {
              dt <- (chaz * newx[indx2, ]) - xbar
              se[indx2] <- sqrt(varh + rowSums((dt %*% object$var) * dt)) * newrisk[indx2]
            }
          }
          else {
            j2 <- approx(afit$time, 1:afit.n, newy[indx2, 2], method = "constant", f = 0, yleft = 0, yright = afit.n)$y
            chaz2 <- approx(-afit$time, afit$cumhaz, -newy[indx2, 2], method = "constant", rule = 2, f = 0)$y
            chaz2 <- c(0, afit$cumhaz)[j2 + 1]
            pred[indx2] <- (chaz2 - chaz) * newrisk[indx2]
            if (se.fit) {
              varh2 <- c(0, cumsum(afit$varhaz))[j1 + 1]
              xbar2 <- rbind(0, afit$xbar)[j1 + 1, , drop = F]
              dt <- (chaz * newx[indx2, ]) - xbar
              dt2 <- (chaz2 * newx[indx2, ]) - xbar2
              v2 <- varh2 + rowSums((dt2 %*% object$var) * dt2)
              v1 <- varh + rowSums((dt %*% object$var) * dt)
              se[indx2] <- sqrt(v2 - v1) * risk[indx2]
            }
          }
        }
      }
    }
    if (survival) {
      if (se.fit) se <- se * exp(-pred)
      pred <- exp(-pred)
    }
  }
  else {
    if (is.null(object$coefficients))
      coef <- numeric(0)
    else {
      coef <- ifelse(is.na(object$coefficients), 0, object$coefficients)
    }
    if (missing(newdata)) {
      offset <- offset - mean(offset)
      if (has.strata && reference == "strata") {
        xmeans <- rowsum(x * weights, oldstrat)/c(rowsum(weights, oldstrat))
        newx <- x - xmeans[match(oldstrat, row.names(xmeans)),]
      }
      else if (use.x)
        newx <- x - rep(object$means, each = nrow(x))
    }
    else {
      offset <- newoffset - mean(offset)
      if (has.strata && reference == "strata") {
        xmeans <- rowsum(x * weights, oldstrat)/c(rowsum(weights, oldstrat))
        newx <- newx - xmeans[match(newstrat, row.names(xmeans)), ]
      }
      else newx <- newx - rep(object$means, each = nrow(newx))
    }
    if (type == "lp" || type == "risk") {
      if (use.x)
        pred <- drop(newx %*% coef) + offset
      else pred <- object$linear.predictors
      if (se.fit) se <- sqrt(rowSums((newx %*% object$var) * newx))
      if (type == "risk") {
        pred <- exp(pred)
        if (se.fit) se <- se * sqrt(pred)
      }
    }
    else if (type == "terms") {
      asgn <- object$assign
      nterms <- length(asgn)
      pred <- matrix(ncol = nterms, nrow = NROW(newx))
      dimnames(pred) <- list(rownames(newx), names(asgn))
      if (se.fit) se <- pred
      for (i in 1:nterms) {
        tt <- asgn[[i]]
        tt <- tt[!is.na(object$coefficients[tt])]
        xtt <- newx[, tt, drop = F]
        pred[, i] <- xtt %*% object$coefficient[tt]
        if (se.fit) se[, i] <- sqrt(rowSums((xtt %*% object$var[tt, tt]) * xtt))
      }
      pred <- pred[, terms, drop = F]
      if (se.fit) se <- se[, terms, drop = F]
      attr(pred, "constant") <- sum(object$coefficients * object$means, na.rm = T)
    }
  }
  if (type != "terms") {
    pred <- drop(pred)
    if (se.fit) se <- drop(se)
  }
  if (!is.null(na.action.used)) {
    pred <- napredict(na.action.used, pred)
    if (is.matrix(pred))
      n <- nrow(pred)
    else n <- length(pred)
    if (se.fit) se <- napredict(na.action.used, se)
  }
  if (!missing(collapse) && !is.null(collapse)) {
    if (length(collapse) != n2) stop("Collapse vector is the wrong length")
    pred <- rowsum(pred, collapse)
    if (se.fit) se <- sqrt(rowsum(se^2, collapse))
    if (type != "terms") {
      pred <- drop(pred)
      if (se.fit) se <- drop(se)
    }
  }
  if (se.fit)
    list(fit = pred, se.fit = se)
  else pred
}

