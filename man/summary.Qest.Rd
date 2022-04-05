\name{summary.Qest}
\alias{summary.Qest}
\title{
Summarizing Q-estimators
}
\description{
Summary method for class \dQuote{Qest}.
}
\usage{
\method{summary}{Qest}(object, covar = FALSE, \ldots)
}
\arguments{
  \item{object}{
an object of class \dQuote{Qest}.
}
  \item{covar}{
logical; if TRUE, the variance covariance matrix of the estimated parameters is returned.
}
  \item{\dots}{
for future methods.
}
}
\details{
This function returns a summary of the most relevant information on model parameters, standard errors, and convergence status.
}
\value{
The function \code{summary.Qest} computes and returns a list of summary statistics of the fitted model given in object, using the "call" and "terms" from its argument, plus

\item{coefficients}{a matrix with 4 columns reporting the estimated coefficients, the estimated standard errors, the corresponding z-values (coef/se), and the two-sided p-values.}
\item{obj.function}{the value of the minimized loss function (see \code{\link{Qest}} for details).}
\item{n}{the number of observations.}
\item{npar}{the number of free parameters.}
\item{iter}{the number of iterations.}
\item{covar}{only if \kbd{covar = TRUE}, the estimated covariance matrix.}
\item{call}{the matched call.}
\item{type}{a character string defined as follows: \code{"c"} for right-censored data; \code{"ct"} for left-truncated, right-censored data; and \code{"u"} otherwise.}
}
\references{
Sottile G, and Frumento P (2022). \emph{Robust estimation and regression with parametric quantile functions.} Computational Statistics and Data Analysis. <doi:10.1016/j.csda.2022.107471>
}
\author{
Gianluca Sottile <gianluca.sottile@unipa.it>
}
\seealso{
\code{\link{Qest}}, for model fitting.
}
\examples{

# Quantile function of an Exponential model
Qexp <- function(theta, tau, data){
  qexp(tau, exp(theta))
}

y <- rexp(100, exp(1))
m1 <- Qest(y ~ 1, Q = Qexp, start = NA)
summary(m1)
summary(m1, covar = TRUE)
}
\keyword{methods}
