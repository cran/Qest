\name{Qlm.fit}
\alias{Qlm.fit}
\title{
Fitter Functions for Quantile-based Linear Models
}
\description{
This is the basic computing engine called by \dQuote{Qlm} used to fit quantile-based linear models. This function should only be used directly by experienced users.
}
\usage{
Qlm.fit(y, X, w = rep(1, nobs), start = NULL, wtau = NULL,
  control = Qest.control(), \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
vector of observations of length n.
}
  \item{X}{
design matrix of dimension n * p.
}
  \item{w}{
 an optional vector of weights to be used in the fitting process.
}
  \item{start}{
starting values for the parameters in the linear predictor.
}
  \item{wtau}{
an optional function that assigns a different weight to each quantile. By default, all quantiles in (0,1) have the same weight.
}
  \item{control}{
a list of operational parameters. This is usually passed through \code{\link{Qest.control}}.
}
  \item{\ldots}{
additional arguments for \code{wtau}.
}
}

\value{
a \dQuote{list} with components

\item{coefficients}{p vector}
\item{std.errs}{p vector}
\item{covar}{p x p matrix}
\item{dispersion}{estimated dispersion parameter}
\item{residuals}{n vector}
\item{rank}{integer, giving the rank}
\item{fitted.values}{n vector}
\item{qr}{the QR decomposition, see \dQuote{qr}}
\item{df.residual}{degrees of freedom of residuals}
\item{obj.function}{the minimized loss function}
\item{gradient}{p vector}
\item{hessian}{p x p matrix}
\item{convergence}{logical. The convergence status}
\item{n.it}{the number of iterations}
\item{control}{control elements}
}
\references{
Sottile G, and Frumento P (2022). \emph{Robust estimation and regression with parametric quantile functions.} Computational Statistics and Data Analysis. <doi:10.1016/j.csda.2022.107471>
}
\author{
Gianluca Sottile <gianluca.sottile@unipa.it>, Paolo Frumento <paolo.frumento@unipi.it>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Qlm}}
}
\examples{
# Ex. 1 Normal model

set.seed(1234)
n <- 100
x1 <- rnorm(n)
x2 <- runif(n,0,3)
y <- rnorm(n, 4 + x1 + 2*x2, 1)
X <- cbind(1, x1, x2)
w <- rep.int(1, n)

m <- Qlm.fit(y = y, X = X, w = w, control = Qest.control(display = TRUE))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
