\name{Qest}
\alias{Qest}
\title{
Q-Estimation
}
\description{
An implementation of the quantile-based estimators described in Sottile and Frumento (2022).
}
\usage{
Qest(formula, Q, weights, start, data, ntau = 199, wtau = NULL,
  control = Qest.control(), \ldots)
}
\arguments{
  \item{formula}{
a two-sided formula of the form \code{y ~ x}. Note that the parametric model is identified through \code{Q}, and not through \code{formula}, that only identifies the response and the predictors. Use \code{Surv(time, event)}, if the data are right-censored, and \code{Surv(start, stop, event)}, if the data are right-censored and left-truncated (\code{start < stop}, \code{start} can be \code{-Inf}).
}
  \item{Q}{
a parametric quantile function of the form \code{Q(theta, tau, data)}. Alternatively, a
    character string naming a \code{Qfamily} function, a \code{Qfamily} function itself, or the
    result of a call to a \code{Qfamily} function. See \code{\link{Qfamily}} for details.
}
  \item{weights}{
 an optional vector of weights to be used in the fitting process. The weights will always be normalized to sum to the sample size. This implies that, for example, using double weights will \emph{not} halve the standard errors.
}
  \item{start}{
a vector of starting values. \code{NA}s are allowed, but will be internally replaced by zeroes. Make sure that the quantile function is well-defined at \code{theta = start}. The size of \code{start} is also used to identify the number of parameters in the model. You \emph{must} supply starting points, unless you are fitting a model defined by a \code{Qfamily}.
}
  \item{data}{
 an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{Qest} is called.
}
  \item{ntau}{
the number of points for numerical integration (see \dQuote{Details}). Default \code{ntau = 199}.
}
  \item{wtau}{
an optional function that assigns a different weight to each quantile. By default, all quantiles in (0,1) have the same weight. Please check the documentation of \code{\link{wtrunc}} for built-in weighting functions.
}
  \item{control}{
a list of operational parameters. This is usually passed through \code{\link{Qest.control}}.
}
  \item{\ldots}{
additional arguments for \code{wtau} and \code{Q}.
}
}
\details{
A parametric model, \eqn{Q(\tau | \theta, x)}, is used to describe the conditional quantile function of an outcome \eqn{Y}, given a vector \eqn{x} of covariates. The model parameters, \eqn{\theta}, are estimated by minimizing the (weighted) integral, with respect to \eqn{\tau}, of the loss function of standard quantile regression. If the data are censored or truncated, \eqn{\theta} is estimated by solving a set of estimating equations. In either case, numerical integration is required to calculate the objective function: a grid of \code{ntau} points in \code{(0,1)} is used.  The estimation algorithm is briefly described in the documentation of \code{\link{Qest.control}}.

The optional argument \code{wtau} can be used to attribute a different weight to each quantile. Although it is possible to choose \code{wtau} to be a discontinuous function (e.g., \code{wtau = function(tau){tau < 0.95}}), this may occasionally result in poorly estimated standard errors.

The quantile function \code{Q} must have at least the following three arguments: \code{theta, tau, data}, in this order. The first argument, \code{theta}, is a vector (not a matrix) of parameters' values. The second argument, \code{tau}, is the order of the quantile. When \code{Q} receives a \code{n*ntau} matrix of \code{tau} values, it must return a \code{n*ntau} matrix of quantiles. The third argument, \code{data}, is a data frame that includes the predictors used by \code{Q}.

If \code{Q} is identified by one \code{\link{Qfamily}}, everything becomes much simpler. It is not necessary to implement your own quantile function, and the starting points are not required. Note that \code{ntau} is ignored if \code{Q = Qnorm} or \code{Q = Qunif}.

Please check the documentation of \code{\link{Qfamily}} to see the available built-in distributions. A convenient Q-based implementation of the standard linear regression model is offered by \code{\link{Qlm}}. Proportional hazards models are implemented in \code{\link{Qcoxph}}.
}
\value{
a list with the following elements:
\item{coefficients}{a named vector of coefficients.}
\item{std.errs}{a named vector of estimated standard errors.}
\item{covar}{the estimated covariance matrix of the estimators.}
\item{obj.function}{the value of the minimized loss function. If the data are censored or truncated, a meaningful loss function which, however, is not the function being minimized (see \dQuote{Note}).}
\item{ee}{the values of the estimating equations at the solution. If the data are neither censored nor truncated, the partial derivatives of the loss function.}
\item{jacobian}{the jacobian at the solution. If the data are neither censored nor truncated, the matrix of second derivatives of the loss function.}
\item{CDF, PDF}{the fitted values of the cumulative distribution function (CDF) and the probability density function (PDF).}
\item{converged}{logical. The convergence status.}
\item{n.it}{the number of iterations.}
\item{internal}{internal elements.}
\item{call}{the matched call.}
}
\references{
Sottile G, and Frumento P (2022). \emph{Robust estimation and regression with parametric quantile functions.} Computational Statistics and Data Analysis. <doi:10.1016/j.csda.2022.107471>
}
\author{
Paolo Frumento <paolo.frumento@unipi.it>, Gianluca Sottile <gianluca.sottile@unipa.it>
}
\note{
NOTE 1. If the data are censored or truncated, estimation is carried out by solving estimating equations, and no associated loss is present. In this case, a meaningful value of \code{obj.function} is the integrated loss [equation 1 of Sottile and Frumento (2022)] in which the indicator function \eqn{I(y \le Q(\tau | \theta, x))} has been replaced with one of the expressions presented in equations 6 and 7 of the paper. The resulting loss, however, is not the function being minimized.


NOTE 2. To prevent computational problems, avoid situations in which some of the estimated parameters are expected to be very small or very large. For example, standardize the predictors, and normalize the response. Avoid as much as possible parameters with bounded support. For example, model a variance/rate/shape parameter on the log scale, e.g., \eqn{\sigma = exp(\theta)}. Carefully select the starting points, and make sure that \code{Q(start, ...)} is well-defined. If \code{Q} is identified by one \code{\link{Qfamily}}, all these recommendations can be ignored.


NOTE 3. You should \emph{not} use \code{Qest} to fit parametric models describing discrete distributions, where the quantile function is piecewise constant. You can try, but the optimization algorithm will most likely fail. The predefined family \code{Qpois} allows to fit a Poisson distribution by using a continuous version of its quantile function (see \code{\link{Qfamily}}).
}
\seealso{
\code{\link{Qest.control}}, for operational parameters, and \code{\link{summary.Qest}}, for model summary. \code{\link{Qfamily}}, for the available built-in distributions. \code{\link{wtrunc}} for built-in weighting functions (\code{wtau} argument). \code{\link{Qlm}}, for Q-estimation of the standard normal (linear) regression model; \code{\link{Qcoxph}}, for proportional hazards models.
}
\examples{
# Ex1. Normal model

# Quantile function of a linear model
Qlinmod <- function(theta, tau, data){
  sigma <- exp(theta[1])
  beta <- theta[-1]
  X <- model.matrix( ~ x1 + x2, data = data)
  qnorm(tau, X \%*\% beta, sigma)
}

n <- 100
x1 <- rnorm(n)
x2 <- runif(n,0,3)
theta <- c(1,4,1,2)
y <- Qlinmod(theta, runif(n), data.frame(x1,x2)) # generate the data

m1 <- Qest(y ~ x1 + x2, Q = Qlinmod, start = c(NA,NA,NA,NA)) # User-defined quantile function
summary(m1)

m2 <- Qest(y ~ x1 + x2, Q = Qnorm) # Qfamily
summary(m2)

m3 <- Qlm(y ~ x1 + x2)
summary(m3) # using 'Qlm' is much simpler and faster, with identical results




\donttest{
# Ex2. Weibull model with proportional hazards

# Quantile function
QWeibPH <- function(theta, tau, data){
  shape <- exp(theta[1])
  beta <- theta[-1]
  X <- model.matrix(~ x1 + x2, data = data)
  qweibull(tau, shape = shape, scale = (1/exp(X \%*\% beta))^(1/shape))
}

n <- 100
x1 <- rbinom(n,1,0.5)
x2 <- runif(n,0,3)
theta <- c(2,-0.5,1,1)

t <- QWeibPH(theta, runif(n), data.frame(x1,x2)) # time-to-event
c <- runif(n,0.5,1.5) # censoring variable
y <- pmin(t,c) # observed response
d <- (t <= c) # event indicator

m1 <- Qest(Surv(y,d) ~ x1 + x2, Q = QWeibPH, start = c(NA,NA,NA,NA))
summary(m1)

m2 <- Qcoxph(Surv(y,d) ~ x1 + x2)
summary(m2) # using 'Qcoxph' is much simpler and faster (but not identical)




# Ex3. A Gamma model

# Quantile function
Qgm <- function(theta, tau, data){
  a <- exp(theta[1])
  b <- exp(theta[2])
  qgamma(tau, shape = a, scale = b)
}
n <- 100
theta <- c(2,-1)
y <- rgamma(n, shape = exp(theta[1]), scale = exp(theta[2]))

m1 <- Qest(y ~ 1, Q = Qgm, start = c(NA, NA)) # User-defined quantile function
m2 <- Qest(y ~ 1, Q = Qgamma) # Qfamily
m3 <- Qest(y ~ 1, Q = Qgamma, wtau = function(tau, h) dnorm((tau - 0.5)/h), h = 0.2)
# In m3, more weight is assigned to quantiles around the median




# Ex4. A Poisson model

# Quantile function
n <- 100
x1 <- runif(n)
x2 <- rbinom(n,1,0.5)
y <- rpois(n, exp(1.5 -0.5*x1 + x2))
m1 <- Qest(y ~ x1 + x2, Q = Qpois) # Use a Qfamily! See "Note"
m2 <- Qest(y + runif(n) ~ x1 + x2, Q = Qpois) # Use jittering! See the documentation of "Qfamily"
}
}

\keyword{models}
\keyword{regression}
