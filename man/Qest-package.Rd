\name{Qest-package}
\alias{Qest-package}
\docType{package}
\title{
\packageTitle{Qest}
}
\description{
\packageDescription{Qest}
}
\details{

\tabular{ll}{
Package: \tab Qest\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2022-03-04\cr
License: \tab GPL-2\cr
}

The DESCRIPTION file:
\packageDESCRIPTION{Qest}
\packageIndices{Qest}
}
\author{
\packageAuthor{Qest}

Maintainer: \packageMaintainer{Qest}
}
\references{
Sottile G, and Frumento P (2022). \emph{Robust estimation and regression with parametric quantile functions.} Computational Statistics and Data Analysis. <doi:10.1016/j.csda.2022.107471>
}

\keyword{ package }
\seealso{
\code{\link{Qest}}, \code{\link{Qlm}}, \code{\link{Qcoxph}}
}
\examples{
\dontrun{
Qest(y ~ x, Q, start) # General-purpose Q-estimator
Qlm(y ~ x) # Q-estimation of linear models
Qcoxph(Surv(time, event) ~ x) # Q-estimation of proportional hazards models
}
}
