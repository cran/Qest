\name{wtrunc}
\alias{wtrunc}

\title{Weighting Function for \code{Qest}, \code{Qlm}, and \code{Qcoxph}.}
\description{
This function can be used within a call to \code{Qest}, \code{Qlm}, or \code{Qcoxph} to assign a different weight to each quantile.
}
\usage{
wtrunc(tau, delta.left = 0.05, delta.right = 0.05, smooth = FALSE, sigma = 0.01)
}
\arguments{
  \item{tau}{a vector of quantiles.}
  \item{delta.left, delta.right}{proportion of quantiles to be removed from the left and righ tail. The weighting function is \code{1} in the interval \code{(delta.left, 1 - delta.right)}, and zero elsewhere. Default is \code{delta.left = 0.05} and \code{delta.right = 0.05}. When a weighting function is used to counteract the effect of extreme observations, \code{delta.left} is a guess for the proportion of outliers on the left tail; and \code{delta.right} is a guess for the proportion of outliers on the right tail.}
  \item{smooth}{if \code{smooth = TRUE} the indicator functions used to construct \code{wtrunc(tau)} are replaced by integrated Gaussian kernels. Default \code{smooth = FALSE}.}
  \item{sigma}{the bandwith of a Gaussian kernel. This parameter controls the smoothness of the weighting function, and is ignored if \code{smooth = FALSE}. Default \code{sigma = 0.01}.}
}

\details{
 Within a call to \code{\link{Qest}}, \code{\link{Qlm}}, or \code{\link{Qcoxph}}, one may want to assign a different weight to each quantile through the optional argument \code{wtau}. This can be done for reasons of efficiency, or to counteract the presence of outliers. While \code{wtau} can be any user-defined function, one can use \code{wtrunc} as a shortcut to construct a weighting function that truncates a certain subset of quantiles in the tails of the distribution. For instance, the estimator defined by \code{Qest(\ldots, wtau = wtrunc, delta.left = 0.05, delta.right = 0.1)} only uses quantiles in the interval (0.05, 0.90) to fit the model. In this example, \code{delta.left = 0.05} is a guess for the proportion of outliers on the left tail; and \code{delta.right} is a guess for the proportion of outliers on the right tail.
 Use \code{smooth = TRUE} to replace the indicator functions involved in \code{wtrunc} with smooth functions. Introducing a weighting function that only assigns a positive weight to the quantiles that correspond to the ``healthy'' part of the distribution allows to deal with any level of contamination by outliers.
}

\value{
  A vector of weights assigned to each quantile.
}
\author{
  Gianluca Sottile <gianluca.sottile@unipa.it>, Paolo Frumento <paolo.frumento@unipi.it>
}
\seealso{
  \code{\link{Qest}}, \code{\link{Qlm}}, \code{\link{Qcoxph}}.
}
\examples{

\dontrun{
taus <- seq(0, 1, length.out = 1000)

### zero weight to quantiles above 0.95
plot(taus, wtrunc(taus, delta.left = 0, delta.right = 0.05),
  type = "l", lwd = 1.5)
# smooth counterpart
lines(taus, wtrunc(taus, delta.left = 0, delta.right = 0.05,
  smooth = TRUE, sigma = .01), col = 2, lwd = 1.5)
lines(taus, wtrunc(taus, delta.left = 0, delta.right = 0.05,
  smooth = TRUE, sigma = .05), col = 3, lwd = 1.5)

### zero weight to quantiles below 0.05
plot(taus, wtrunc(taus, delta.left = 0.05, delta.right = 0),
  type = "l", lwd = 1.5)
# smooth counterpart
lines(taus, wtrunc(taus, delta.left = 0.05, delta.right = 0,
  smooth = TRUE, sigma = .01), col = 2, lwd = 1.5)
lines(taus, wtrunc(taus, delta.left = 0.05, delta.right = 0,
  smooth = TRUE, sigma = .05), col = 3, lwd = 1.5)


### zero weight to quantiles below 0.05 and above 0.90
plot(taus, wtrunc(taus, delta.left = 0.05, delta.right = 0.10),
  type = "l", lwd = 1.5)
# smooth counterpart
lines(taus, wtrunc(taus, delta.left = 0.05, delta.right = 0.10,
  smooth = TRUE, sigma = .01), col = 2, lwd = 1.5)
lines(taus, wtrunc(taus, delta.left = 0.05, delta.right = 0.10,
  smooth = TRUE, sigma = .05), col = 3, lwd = 1.5)


### Use wtrunc in Qest, Qlm, Qcoxph
# Qest(..., wtau = wtrunc, delta.left = 0.05, delta.right = 0.1)
}
}
\keyword{models}








