\name{internals}
\alias{internals}
\alias{print.Qest}
\alias{print.summary.Qest}
\alias{confint.Qest}
\alias{vcov.Qest}

\alias{summary.Qlm}
\alias{print.summary.Qlm}
\alias{vcov.Qlm}

\alias{print.Qcoxph}
\alias{summary.Qcoxph}
\alias{print.summary.Qcoxph}
\alias{basehaz.Qcoxph}
\alias{survfit.Qcoxph}
\alias{residuals.Qcoxph}
\alias{predict.Qcoxph}

\title{Internal Functions}
\description{
Functions for internal use only, or not yet documented.
}
\usage{
expit(x)
logit(x)
pmax0(x)
num.fun(dx,fx)
format.perc(probs, digits)
Ltau(opt, tau)
minabs(x1,x2)
invJ(J, type)
tensorX(X)
buildTau(ntau, wtau = NULL, nobs, wtauoptions = NULL)
callwtau(wtau, tau, opt)
callQ(Q, theta, tau, data)
start.Qest(z, y, d, x, w, tau, Q, opt, start, data, type, control)
start.Qest.family(z, y, d, x, w, tau, wtau, wtauoptions,
  Q, opt, start, data, type, control)
rq.fit.br2(x, y, tau = 0.5)

% Optimizers for Qest and Qcoxph
Qest.sgs.internal(theta, type, tol, maxit, alpha0, ee, display, eps, n.it, \ldots)
Qest.gs.internal(theta, type, tol, maxit, alpha0, ee, display, eps, n.it, \ldots)
Qest.gs(theta, type, tol, maxit, alpha0, ee, display, eps, \ldots)
Qest.newton(theta, type, tol, maxit, safeit, alpha0, display, eps, \ldots)

% Qlm utilities
Qlm.bfun(wtau, \ldots)
start.Qlm(x, y, w, start, ok, Stats)
scalevars.Qlm(X,y)
descalecoef.Qlm(theta, Stats)
% Optimizers for Qlm
Qlm.sgs.internal(theta, type, tol, maxit, alpha0, ee, display, n.it, y, X, w, bfun)
Qlm.gs.internal(theta, type, tol, maxit, alpha0, ee, display, n.it, y, X, w, bfun)
Qlm.gs(theta, type, tol, maxit, alpha0, ee, display, y, X, w, bfun)
Qlm.newton(theta, type = "u", tol, maxit, safeit, alpha0, display, y, X, w, bfun)

% Qcoxph utilities
plfcox(y, knots, deriv = 0)
scalevars.Qcoxph(X,z,y,knots)
descalecoef.Qcoxph(theta, Stats)
check.singularities(X, scaleVars)
starting.points.Qcox(X, Y, n, w, mf, knots)
adjust.coef(theta)
agsurv.Qcoxph(y, x, wt, risk, fit)
coxsurv.fit.Qcoxph(ctype, stype, se.fit, varmat, cluster,
  y, x, wt, risk, position, strata, oldid, y2, x2, risk2,
  strata2, id2, unlist = TRUE, fit)
seg.lm.fit1(y,XREG,Z,PSI,return.all.sol=FALSE)

% Qest utilities
gs(theta0, f, \ldots, tol = 1e-4, maxit = 100)
myg(theta, f, f0, eps, \ldots)
dlist(x1,x2)
omega(d, tau, type, Fy, Fz)
choose_eps(Q, theta, y, data, eps0, obj = 0.01)
derQtheta(theta, eps, Q, Q1, data, tau, ind)
der2Qtheta(theta, eps, Q, Qtheta1, data, tau)
intA(A, tau, ifun = TRUE)
findp(y, tau, Q1)
findAp(p,tau,A)

% Qest building blocks
A_beta_fun(BB)
A_gamma_fun(BB)
A_beta_beta_fun(BB)
A_gamma_gamma_mix_fun(BB)
A_gamma_gamma_fun(BB)
A_beta_gamma_fun(BB)
coxBB(theta, y, X, knots, tau)

% Qest Gamma
derQtheta.gamma(Q)
der2Qtheta.gamma(Q, Qtheta)
findAp.gamma(atau, tau, dtau, p, int = TRUE)
QestGamma.ee.u(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE)
QestGamma.ee.c(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE)
QestGamma.ee.ct(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE)

% Qest Poisson
tau.pois(tau)
ppoisC(y, lambda)
dpoisC(y, lambda)
qpoisC.955(z, lambda)
qpoisC.me(log.lambda, A, B)
qpoisC.bisec(tau, lambda)
qpoisC(obj)

derQtheta.pois(Q)
der2Qtheta.pois(Q, Qtheta)
findp.pois(y, tau, Q1, Fy, theta)
findAp.pois(p, tau, A)
QestPois.ee.u(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE)

% Qest Uniform
QestUnif.ee.u(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE)

% Qest Normal
QestNorm.ee.u(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE)
QestNorm.ee.c(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE)
QestNorm.ee.ct(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE)

% Estimating equations
Qest.ee.u(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE)
Qest.ee.c(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE)
Qest.ee.ct(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE)
QCox.ee.c(theta, eps, z, y, d, X, w, knots, tau, J = FALSE, EE)
QCox.ee.ct(theta, eps, z, y, d, X, w, knots, tau, J = FALSE, EE)
Qlm.ee.u(theta, X, w, bfun, EE, J = FALSE)

% Covariance matrix
Qest.covar(fit, eps, w)
Qcox.covar(theta, z, y, d, X, w, knots, tau, type)
Qlm.covar(g.i, w, H)

% Loss function
Loss(w, d, tau, type, Fy, Fz)
coxLoss(theta, z, y, d, X, w, knots, tau, type, Fy, Fz)
qlmLoss(theta, y, X, w, bfun)

\method{print}{Qest}(x, digits = max(3L, getOption("digits") - 3L), \ldots)
\method{print}{summary.Qest}(x, digits = max(3L, getOption("digits") - 3L), \ldots)
\method{confint}{Qest}(object, parm, level = 0.95, \ldots)
\method{vcov}{Qest}(object, \ldots)

\method{summary}{Qlm}(object, correlation = FALSE,
  symbolic.cor = FALSE, \ldots)
\method{print}{summary.Qlm}(x, digits = max(3L, getOption("digits") - 3L),
  symbolic.cor = x$symbolic.cor, signif.stars = getOption("show.signif.stars"),
  \ldots)
\method{vcov}{Qlm}(object, \ldots)

\method{print}{Qcoxph}(x, digits = max(1L, getOption("digits") - 3L),
  signif.stars = FALSE, \ldots)
\method{summary}{Qcoxph}(object, conf.int = 0.95, scale = 1, \ldots)
\method{print}{summary.Qcoxph}(x, digits = max(getOption("digits") - 3, 3),
  signif.stars = getOption  ("show.signif.stars"), \ldots)
\method{basehaz}{Qcoxph}(fit, centered = TRUE, se.fit = FALSE)
\method{survfit}{Qcoxph}(formula, newdata, se.fit = TRUE, conf.int = 0.95,
  individual = FALSE, stype = 2, ctype, conf.type = c("log", "log-log",
  "plain","none", "logit", "arcsin"), censor = TRUE, start.time, id,
  influence = FALSE, na.action = na.pass, type, \ldots)
\method{residuals}{Qcoxph}(object, type = c("martingale", "deviance", "score",
  "schoenfeld", "dfbeta", "dfbetas", "scaledsch", "partial"),
  collapse = FALSE, weighted = FALSE, \ldots)
\method{predict}{Qcoxph}(object, newdata, type = c("lp", "risk", "expected",
  "terms", "survival"), se.fit = FALSE, na.action = na.pass,
  terms = names(object$assign), collapse, reference = c("strata", "sample"),
  \ldots)

}
\keyword{internal}
\value{
No return value, internal functions.
}
