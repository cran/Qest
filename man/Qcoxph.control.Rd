\name{Qcoxph.control}
\alias{Qcoxph.control}
\title{
Auxiliary for Controlling Qcoxph Fitting
}
\description{
Auxiliary function for controlling \code{Qcoxph} fitting. Estimation proceeds in three steps: (i) evaluation of starting points; (iia) stochastic gradient-based optimization (iib) standard gradient-based optimization; and (iii) Newton-Raphson. Step (i) is based on a preliminary fit of a Cox model via \code{coxph}. Steps (iia) and (iib) find an approximate solution, and make sure that the Jacobian matrix is well-defined. Finally, step (iii) finds a more precise solution.
}
\usage{
Qcoxph.control(tol = 1e-8, maxit, safeit, alpha0, display = FALSE)
}
\arguments{
  \item{tol}{
tolerance for convergence of Newton-Raphson algorithm, default is 1e-8.
}
  \item{maxit}{
maximum number of iterations of Newton-Raphson algorithm. If not provided, a default is computed as \code{50 + 25*npar}, where \code{npar} is the total number of parameters.
}
\item{safeit}{
maximum number of iterations of gradient-search algorithm. If not provided, a default is computed as \code{10 + 5*npar}, where \code{npar} is the total number of parameters.
}
  \item{alpha0}{
step size for the preliminary gradient-based iterations. If estimation fails, you can try choosing a small value of \code{alpha0}. If \code{alpha0} is missing, an adaptive choiche will be made internally.
}
  \item{display}{
Logical. If \code{TRUE}, tracing information on the progress of the optimization is printed on screen. Default is \code{FALSE}.
}
}
\details{
If called with no arguments, \code{Qcoxph.control()} returns a list with the current settings of these parameters. Any arguments included in the call sets those parameters to the new values, and then silently returns.
}
\value{
A list with named elements as in the argument list
}
\seealso{
\code{\link{Qcoxph}}
}
\author{
Gianluca Sottile <gianluca.sottile@unipa.it>
Paolo Frumento <paolo.frumento@unipi.it>
}
