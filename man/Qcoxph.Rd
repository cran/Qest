\name{Qcoxph}
\alias{Qcoxph}
\title{
Q-Estimation of Proportional Hazards Regression Models
}
\description{
Fit proportional hazards regression models using Q-estimation.
}
\usage{
Qcoxph(formula, weights, start, data, knots, wtau = NULL,
    control = Qcoxph.control(), \ldots)
}
\arguments{
  \item{formula}{
an object of class \dQuote{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted. Use \code{Surv(time, event) ~ x}, if the data are right-censored, and \code{Surv(time, time2, event)  ~ x}, if the data are right-censored and left-truncated (\code{time < time2}, \code{time} can be \code{-Inf}).
}
  \item{weights}{
 an optional vector of weights to be used in the fitting process. The weights will always be normalized to sum to the sample size. This implies that, for example, using double weights will \emph{not} halve the standard errors.
}
  \item{start}{
optional starting values for the coefficients of the linear predictor.
}
  \item{data}{
an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{Qcoxph} is called.
}
  \item{knots}{
knots to create the basis of a piecewise linear function. If \code{knots} is a vector of at least two elements, it is used to identify the exact position of \emph{all} knots, including boundaries. If \code{knots} is a scalar, its value is used to determine the number of internal knots (\code{knots = 0} is allowed, and fits an Exponential model). If \code{knots} is missing, by default \code{max(1, min(floor(n.events/30), 3))} internal knots are used. Unless a vector of \code{knots} is provided by the user, the \dQuote{optimal} position of the knots will be identified using the method described in Muggeo (2008). If this fails, the knots are positioned at the empirical quantiles of the observed events.
}
  \item{wtau}{
an optional function that assigns a different weight to each quantile. By default, all quantiles in (0,1) have the same weight. Please check the documentation of \code{\link{wtrunc}} for built-in weighting functions.
}
  \item{control}{
a list of operational parameters. This is usually passed through \code{\link{Qcoxph.control}}.
}
  \item{\ldots}{
additional arguments for \code{wtau}.
}
}
\details{
This function estimates a proportional hazards model, allowing for right-censored and left-truncated data. The syntax and output of \code{Qcoxph} are almost identical to those of \code{\link{coxph}}, but the parameters are estimated using Q-estimation (Sottile and Frumento, 2020). This method can be used to obtain outlier-robust estimators of the regression coefficients.

The quantile function of a proportional hazards model is given by
\deqn{Q(\tau | x) = H0^{-1}(-exp{-x'\beta} log(1 - \tau))}
where \eqn{H0} is the baseline cumulative hazard function. In \code{Qcoxph}, \eqn{H0} is parametrized by a piecewise linear function identified by the provided \code{knots}. As the number of knots increases, the baseline hazard becomes arbitrarily flexible.

Estimation is carried out by finding the zeroes of a set of integrals equation. The optional argument \code{wtau} permits assigning a different weight to each quantile in (0,1). It is possible to choose \code{wtau} to be a discontinuous function (e.g., \code{wtau = function(tau){tau < 0.95}}). However, this may occasionally result in poorly estimated of the standard errors.

The estimation algorithm is briefly described in the documentation of \code{\link{Qcoxph.control}}.
}
\value{
an object of classes \dQuote{Qcoxph}, \dQuote{coxph}, and \dQuote{Qest}. See \code{\link{coxph.object}} for details. All the S3 methods that are available for \dQuote{coxph} objects will also work with a \dQuote{Qcoxph} object.

An object of class \dQuote{Qcoxph} is a list containing at least the following components:

\item{coefficients}{a named vector of coefficients.}
\item{var}{the covariance matrix of the coefficients.}
\item{iter}{number of iterations used.}
\item{linear.predictors}{the vector of linear predictors, one per subject. Note that this vector has not been centered, see \code{\link{predict.coxph}} for details.}
\item{residuals}{the martingale residuals.}
\item{means}{vector of column means of the X matrix. Subsequent survival curves are adjusted to this value.}
\item{n}{the number of observations used in the fit.}
\item{nevent}{the number of events used in the fit.}
\item{concordance}{a vector of length 6, containing the number of pairs that are concordant, discordant, tied on x, tied on y, and tied on both, followed by the standard error of the concordance statistic.}
\item{terms, assign, formula, call, y}{other objects used for prediction.}
\item{obj.function}{the objective function of the model. Please, interpret with care: read the note in the documentation of \code{\link{Qest}}.}
\item{internal}{internal objects.}
}
\references{
Sottile G, and Frumento P (2022). \emph{Robust estimation and regression with parametric quantile functions.} Computational Statistics and Data Analysis. <doi:10.1016/j.csda.2022.107471>

Muggeo VMR (2008). \emph{Segmented: an R package to fit regression models with broken-line relationships.} R News 8/1, 20–25.
}
\author{
Paolo Frumento <paolo.frumento@unipi.it>, Gianluca Sottile <gianluca.sottile@unipa.it>
}
\seealso{
\code{\link{Qest}}, for general Q-estimation, and \code{\link{Qlm}}, for Q-estimation of linear models.
}
\examples{

# A proportional-hazards Weibull model

n <- 100
x <- runif(n,0,3)
shape <- 2
t <- rweibull(n, shape = shape, scale = (1/exp(2 + 2*x))^(1/shape)) # time-to-event
c <- runif(n,0,1) # censoring variable
y <- pmin(t,c) # observed response
d <- (t <= c) # event indicator

require(survival)
m1 <- coxph(Surv(y,d) ~ x) # standard Cox model
m2 <- Qcoxph(Surv(y,d) ~ x) # Q-estimator

}
\keyword{models}
\keyword{regression}
