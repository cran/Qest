\name{invQ}
\alias{invQ}
\title{
Inverse of Quantile Function
}
\description{
Auxiliary function to compute cumulative distribution function (CDF) by inverting a quantile function.
}
\usage{
invQ(Q, theta, y, data, n.it = 17)
}
\arguments{
  \item{Q}{
any parametric quantile function of the form \code{Q(theta, tau, data)}.
}
  \item{theta}{
a vector of model parameters.
}
  \item{y}{
vector of observations to evaluate the CDF.
}
  \item{data}{
data frame containing the variables used in the Q() function.
}
  \item{n.it}{
the number of iterations (see \dQuote{details}).
}
}
\details{
Given a parametric quantile function \eqn{Q(\tau | \theta)}, the CDF is defined as \eqn{F(y | \theta) = Q^{-1}(y | \theta)}. Alternatively, \eqn{F(y | \theta)} corresponds to the value \eqn{\tau*} such that \eqn{Q(\tau* | \theta) = y}. Starting from \eqn{\tau = 0.5}, a bisection algorithm is used to evaluate numerically \eqn{\tau*}. The maximum error is given by \code{1/2^(n.it + 1)}.
}
\value{
a vector of CDF values.
}
\author{
Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}


\seealso{
\code{\link{Qest}}
}
\examples{
# Ex. 1 Normal model

# Quantile function of a linear model.
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

# generate the data
U <- runif(n)
y <- Qlinmod(theta, U, data.frame(x1,x2))

# Given y and theta, evaluate U = F(y)
invQ(Qlinmod, theta, y, data.frame(x1,x2))

}
\keyword{methods}
