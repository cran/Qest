\name{Qest.control}
\alias{Qest.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auxiliary for Controlling Qest Fitting
}
\usage{
Qest.control(tol = 1e-8, maxit, safeit, alpha0, display = FALSE, restart = FALSE)
}
\arguments{
  \item{tol}{
tolerance for convergence of Newton-Raphson algorithm, default is 1e-8.
}
  \item{maxit}{
maximum number of iterations of Newton-Raphson algorithm. If not provided, a default is computed as \code{50 + 25*npar}, where \code{npar} is the number of parameters.
}
\item{safeit}{
maximum number of iterations of gradient-search algorithm. If not provided, a default is computed as \code{10 + 5*npar}, where \code{npar} is the number of parameters.
}
  \item{alpha0}{
step size for the preliminary gradient-based iterations. If estimation fails, you can try choosing a small value of \code{alpha0}. If \code{alpha0} is missing, an adaptive choiche will be made internally.
}
  \item{display}{
Logical. If \code{TRUE}, tracing information on the progress of the optimization is printed on screen. Default is \code{FALSE}.
}
  \item{restart}{
Logical. If \code{FALSE} (the default), step (i) is not performed, and the provided starting points are directly passed to step (ii). This may save you some time, but is not recommended unless you are confident about your choice of initial values. When \code{restart = TRUE}, the provided starting points are used to initialize step (i).
}
}
\description{
Auxiliary function for controlling \code{Qest} fitting. Estimation proceeds in three steps: (i) evaluation of starting points; (iia) stochastic gradient-based optimization (iib) standard gradient-based optimization; and (iii) Newton-Raphson. Step (i) is initialized at the provided starting values (the \code{start} argument of \code{Qest}), and utilizes a preliminary flexible model, estimated with \code{\link{pchreg}}, to generate a cheap guess of the model parameters. If you have good starting points, you can skip step (i) by setting \code{restart = FALSE}. Steps (iia) and (iib) find an approximate solution, and make sure that the Jacobian matrix is well-defined. Finally, step (iii) finds a more precise solution.
}
\note{Step (i) is not performed, and \code{restart} is ignored, if the quantile function is one of the available \code{\link{Qfamily}}.}
\details{
If called with no arguments, \code{Qest.control()} returns a list with the current settings of these parameters. Any arguments included in the call sets those parameters to the new values, and then silently returns.
}
\value{
A list with named elements as in the argument list
}
\seealso{
\code{\link{Qest}} and \code{\link{Qlm}}
}
\author{
Gianluca Sottile <gianluca.sottile@unipa.it>
Paolo Frumento <paolo.frumento@unipi.it>
}
