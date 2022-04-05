\name{Qlm}
\alias{Qlm}
\title{
Q-Estimation of Linear Regression Models
}
\description{
Use Q-estimation to fit a Normal model in which the mean is a linear function of the predictors, and the variance is constant.
}
\usage{
Qlm(formula, data, subset, weights, na.action, start = NULL, contrasts = NULL,
  wtau = NULL, control = Qest.control(), \ldots)
}
\arguments{
  \item{formula}{
an object of class \dQuote{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
}
  \item{data}{
an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{Qlm} is called.
}
  \item{subset}{
an optional vector specifying a subset of observations to be used in the fitting process.
}
  \item{weights}{
 an optional vector of weights to be used in the fitting process. The weights will always be normalized to sum to the sample size. This implies that, for example, using double weights will \emph{not} halve the standard errors.
}
  \item{na.action}{
a function which indicates what should happen when the data contain \code{NA}s. See \code{\link{lm}}.
}
  \item{start}{
optional starting values for the regression coefficients.
}
  \item{contrasts}{
an optional list. See the \code{contrasts.arg} of \code{model.matrix.default}.
}
  \item{wtau}{
an optional function that assigns a different weight to each quantile. By default, all quantiles in (0,1) have the same weight. Please check the documentation of \code{\link{wtrunc}} for built-in weighting functions.
}
  \item{control}{
a list of operational parameters. See \code{Qest.control} for details.
}
  \item{\ldots}{
additional arguments for \code{wtau}.
}
}
\details{
This function is used exactly as \code{\link{lm}}, but estimates the model parameters as in \code{\link{Qest}}. Using Q-estimation allows to obtain outlier-robust estimators of the regression coefficients. The optional argument \code{wtau} permits assigning a different weight to each quantile in (0,1). It is possible to choose \code{wtau} to be a discontinuous function (e.g., \code{wtau = function(tau){tau < 0.95}}). However, this may occasionally result in poorly estimated of the standard errors.

Note that \code{Qlm}, like \code{lm}, does not support censored or truncated data.
}
\value{
Qlm returns an object of classes \dQuote{Qlm}, \dQuote{lm}, and \dQuote{Qest}.
The generic accessor functions \code{summary}, \code{coefficients}, \code{fitted.values}, and \code{residuals} can be used to extract infromation from a \dQuote{Qlm} object.

An object of class \dQuote{Qlm} is a list containing at least the following elements:

\item{coefficients}{a named vector of coefficients.}
\item{std.errs}{a named vector of standard errors.}
\item{covar}{the estimated covariance matrix of the estimators.}
\item{dispersion}{the estimated dispersion parameter (residual variance).}
\item{residuals}{the working residuals.}
\item{rank}{the estimated degrees of freedom.}
\item{fitted.values}{the fitted values.}
\item{df.residual}{the residual degrees of freedom.}
\item{obj.function}{the value of the minimized loss function.}
\item{gradient}{the first derivatives of the minimized loss function.}
\item{hessian}{the matrix of second derivatives of the minimized loss function.}
\item{convergence}{logical. The convergence status.}
\item{n.it}{the number of iterations.}
\item{control}{control parameters.}
\item{xlevels}{(only where relevant) a record of the levels of the factors used in fitting.}
\item{call}{the matched call.}
\item{terms}{the \dQuote{terms} object used.}
\item{model}{if requested (the default), the model frame used.}
}
\references{
Sottile G, and Frumento P (2022). \emph{Robust estimation and regression with parametric quantile functions.} Computational Statistics and Data Analysis. <doi:10.1016/j.csda.2022.107471>
}
\author{
Gianluca Sottile <gianluca.sottile@unipa.it>, Paolo Frumento <paolo.frumento@unipi.it>
}
\seealso{
\code{\link{Qest}}, for general Q-estimation.
}
\examples{
set.seed(1234)
n <- 100
x1 <- rnorm(n)
x2 <- runif(n,0,3)
theta <- c(1,4,1,2)
y <- rnorm(n, 4 + x1 + 2*x2, 1)

m1 <- Qlm(y ~ x1 + x2)
summary(m1)
}
\keyword{models}
\keyword{regression}
