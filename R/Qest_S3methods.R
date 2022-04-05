

#### S3 methods for Qest
print.Qest <- function (x, digits = max(3L, getOption("digits") - 3L), ...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(x)
}

summary.Qest <- function (object, covar = FALSE, ...){
  theta <- object$coefficients
  se <- object$std.errs
  qq <- length(theta)
  coe <- matrix(0, nrow = qq, ncol = 4, dimnames=list(names(theta), c("Estimate", "std.err", "z value", "p(>|z|)")))
  coe[, 1] <- theta
  coe[, 2] <- se
  coe[, 3] <- coe[, 1] / coe[, 2]
  coe[, 4] <- 2 * (1 - pnorm(abs(coe[, 3])))

  out <- list()
  out$coefficients <- coe
  out$obj.function <- object$obj.function
  out$n <- length(object$internal$y)
  out$npar <- length(theta[!is.na(theta)])
  out$iter <- object$n.it
  if(covar) out$covar <- object$covar
  out$call <- object$call
  out$type <- object$internal$type
  class(out) <- "summary.Qest"
  out
}

print.summary.Qest <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("n. of iterations:", x$iter, "\n")
  cat("n. of observations:", x$n, "\n")
  cat("n. of parameters:", x$npar, "\n\n")

  cat("######################", "\n")
  cat("######################", "\n\n")

  cat("Coefficients:\n")
  printCoefmat(x$coe, digits = max(3L, getOption("digits") - 3L), signif.stars = TRUE, cs.ind = 1:2, tst.ind = 3, P.values = TRUE, has.Pvalue = TRUE)
  cat("\n")

  cat("######################", "\n")
  cat("######################", "\n")

  cat("\n")
  if(x$type == "u"){m <- "Minimized loss function"}
  else{m <- "Loss (not the function being minimized)"}
  cat(m, x$obj.function)

  if (!is.null(x$covar)) {
    cat("\n\n")
    cat("######################", "\n")
    cat("######################", "\n")

    cat("\n")
    cat("Covar:\n")
    print.default(format(x$covar, digits = digits), print.gap = 2L, quote = FALSE)
  }

  cat("\n\n")

  invisible(x)
}

confint.Qest <- function(object, parm, level = 0.95, ...) {
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm))
    parm <- pnames
  else if (is.numeric(parm))
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- format.perc(a, 3)
  fac <- qnorm(a)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  ses <- object$std.errs[parm]
  ci[] <- cf[parm] + ses %o% fac
  ci
}

vcov.Qest <- function(object, ...) object$covar







