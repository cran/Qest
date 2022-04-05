\name{Qfamily}
\alias{Qfamily}
\alias{Qnorm}
\alias{Qgamma}
\alias{Qpois}
\alias{Qunif}
%\alias{print.family}

\title{Family Objects for Qest}
\usage{
Qnorm()
Qgamma()
Qpois(offset = NULL)
Qunif(min = TRUE)
}
\arguments{
  \item{offset}{an optional vector of offsets for a Poisson model.}
  \item{min}{logical. If \kbd{TRUE}, fit a \eqn{U(a, b)} distribution. If \kbd{FALSE}, fit a \eqn{U(0, b)} distribution.}
}
\description{
  Family objects are used to specify the model to be fitted by \code{\link{Qest}}.
}
\details{
  A \code{Qfamily} object can be used to identify a certain type of distribution within a call to \code{Qest}. You can supply either the name of the family, or the function itself, or a call to it. For example, the following are equivalent: \code{Qest(formula, "Qpois")}, \code{Qest(formula, Qpois)}, and \code{Qest(formula, Qpois())}. The latter syntax can be used to pass additional arguments, if any.

  The \code{Qnorm} family fits a normal homoskedastic model in which the mean is described by a linear predictor. The parameters are: \code{log(sigma), beta}. \code{Qest(formula, Qnorm)} is equivalent to \code{Qlm(formula)}, but returns a very basic output. However, \code{Qest} allows for censored and truncated data, while \code{Qlm} does not.

  The \code{Qgamma} family fits a Gamma distribution in which the log-scale is modeled by a linear predictor. The model parameters are: \code{log(shape), beta}.

  The \code{Qpois} family fits a Poisson distribution in which the log-rate is modeled by a linear predictor. In reality, to obtain a continuous quantile function, \code{qpois} is replaced by the inverse, \emph{with respect to} \eqn{y}, of the upper regularized gamma function, \eqn{Q(y,\lambda)}. It is recommended to apply \code{Qpois} to a jittered response (i.e., \code{y + runif(n)}).

  The \code{Qunif} family fits a Uniform distribution \eqn{U(a,b)} in which both \eqn{a} and \eqn{b} are modeled by linear predictors. The design matrix, however, is the same for \eqn{a} and \eqn{b}. Use \code{Qunif(min = FALSE)} to fit a \eqn{U(0,b)} model. The parameters are: \code{beta_a, beta_b}, or only \code{beta_b} if \code{min = FALSE}.

  The families \code{Qnorm} and \code{Qgamma} can be used when the data are censored or truncated, while \code{Qpois} and \code{Qunif} cannot. All families can be estimated without covariates, using \code{formula = ~ 1}.
}

\value{
  An object of class \code{"Qfamily"} that contains all the necessary information to be passed to \code{Qest}.
}
\author{
  Gianluca Sottile <gianluca.sottile@unipa.it>, Paolo Frumento <paolo.frumento@unipi.it>
}
\seealso{
  \code{\link{Qest}}.
}
\examples{

n <- 250
x <- runif(n)
eta <- 1 + 2*x # linear predictor

# Normal model
y <- rnorm(n, eta, exp(1))
m1 <- Qest(y ~ x, Qnorm)
# Use Qlm(y ~ x) instead!

# Gamma model
y <- rgamma(n, shape = exp(1), scale = exp(eta))
m2 <- Qest(y ~ x, Qgamma)

# Poisson model
y <- rpois(n, exp(eta))
m3 <- Qest(y ~ x, Qpois)
m4 <- Qest(y + runif(n) ~ x, Qpois) # Jittering is recommended

# Uniform model
y <- runif(n, 0, eta)
m5 <- Qest(y ~ x, Qunif(min = TRUE))  # U(a,b)
m6 <- Qest(y ~ x, Qunif(min = FALSE)) # U(0,b)
}
\keyword{models}








