## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4)

## ----setup--------------------------------------------------------------------
library(Qest)
library(survival)

## ----user---------------------------------------------------------------------
Qlinmod <- function(theta, tau, data) {
  X <- model.matrix(~ x1 + x2, data = data)
  qnorm(tau, c(X %*% theta[-1]), exp(theta[1]))
}
set.seed(1)
n <- 200
dat <- data.frame(x1 = rnorm(n), x2 = runif(n, 0, 3))
dat$y <- 4 + dat$x1 + 2 * dat$x2 + rnorm(n)

m1 <- Qest(y ~ x1 + x2, Q = Qlinmod, start = rep(NA, 4), data = dat)
summary(m1)

## ----compare------------------------------------------------------------------
m2 <- Qest(y ~ x1 + x2, Q = Qnorm, data = dat)
m3 <- Qlm(y ~ x1 + x2, data = dat)
cbind(user = coef(m1), Qnorm = coef(m2), Qlm = c(log(m3$sigma), coef(m3)))

## ----weibull------------------------------------------------------------------
QWeibPH <- function(theta, tau, data) {
  shape <- exp(theta[1])
  X <- model.matrix(~ x, data = data)
  qweibull(tau, shape = shape, scale = exp(-c(X %*% theta[-1]) / shape))
}
dat <- data.frame(x = runif(n, 0, 3))
t <- QWeibPH(c(log(2), -0.5, 1), runif(n), dat)
cens <- runif(n, 0, 2)
dat$y <- pmin(t, cens)
dat$d <- as.numeric(t <= cens)
m4 <- Qest(Surv(y, d) ~ x, Q = QWeibPH, start = rep(NA, 3), data = dat)
summary(m4)

## ----plot---------------------------------------------------------------------
plot(m4)

## ----families-----------------------------------------------------------------
x <- runif(n)
yg <- rgamma(n, shape = 2, scale = exp(1 + x))
coef(Qest(yg ~ x, Q = Qgamma))
yp <- rpois(n, exp(1 + x)) + runif(n)
coef(Qest(yp ~ x, Q = Qpois))
yu <- runif(n, 0, 1 + 2 * x)
coef(Qest(yu ~ x, Q = Qunif(min = FALSE)))

## ----predict------------------------------------------------------------------
fg <- Qest(yg ~ x, Q = Qgamma)
predict(fg, newdata = data.frame(x = c(0, 1)), p = c(0.1, 0.5, 0.9))

## ----robust-------------------------------------------------------------------
set.seed(2)
x <- runif(n)
y <- 1 + 2 * x + rnorm(n)
out <- sample(n, 10)
y[out] <- y[out] + 15
rbind(lm = coef(lm(y ~ x)),
      Qlm = coef(Qlm(y ~ x)),
      Qlm.trunc = coef(Qlm(y ~ x, wtau = wtrunc, delta.left = 0.1, delta.right = 0.1)))

## ----cox----------------------------------------------------------------------
set.seed(3)
dat <- data.frame(x = runif(n, 0, 3))
t <- rweibull(n, shape = 2, scale = (1 / exp(2 + 2 * dat$x))^(1 / 2))
cens <- runif(n, 0, 1)
dat$y <- pmin(t, cens); dat$d <- as.numeric(t <= cens)
fc <- Qcoxph(Surv(y, d) ~ x, data = dat)
summary(fc)
cbind(coxph = coef(coxph(Surv(y, d) ~ x, data = dat)), Qcoxph = coef(fc))
sf <- survfit(fc, newdata = data.frame(x = c(0.5, 1.5)))
plot(sf, col = 1:2, conf.int = TRUE, xlab = "time", ylab = "survival")

