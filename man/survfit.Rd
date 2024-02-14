\name{survfit}
\alias{survfit}
\title{Create survival curves}
\description{
This function creates survival curves from either a formula (e.g. the 
Kaplan-Meier), a previously fitted Cox model, or a previously fitted
accelerated failure time model.
}
\usage{
survfit(formula, ...)
}
%- coef, vcov, and confint return "not applicable" errors, on purpose
\arguments{
  \item{formula}{either a formula or a previously fitted model}
  \item{\dots}{other arguments to the specific method}
}
\details{
  A survival curve is based on a tabulation of the number at risk and
  number of events at each unique death time.  When time is a floating
  point number the definition of "unique" is subject to interpretation.
  The code uses factor() to define the set.
  For further details see the documentation for the appropriate method, i.e.,
  \code{?survfit.formula} or \code{?survfit.coxph}.  

  A survfit object may contain a single curve, a set of curves (vector), a
  matrix of curves, or even a 3 way array: \code{dim(fit)} will reveal
  the dimensions.
  Predicted curves from a \code{coxph} model have one row for each
  stratum in the Cox model fit and one column for each specified
  covariate set.
  Curves from a multi-state model have one row for each stratum and
  a column for each state, the strata correspond to predictors on the
  right hand side of the equation.  The default printing and plotting
  order for curves is by column, as with other matrices.
}
\value{
An object of class \code{survfit} containing one or more survival curves.
}
\author{Terry Therneau}
\note{Older releases of the code also allowed the specification for 
a single curve
to omit the right hand of the formula, i.e.,
\code{survfit(Surv(time, status))}, in which case the formula argument
is not actually a formula.
Handling this case required some non-standard and fairly fragile 
manipulations, and this case is no longer supported.
}
\seealso{\code{\link{survfit.formula}},
         \code{\link{survfit.coxph}},
	 \code{\link{survfit.object}}, \code{\link{print.survfit}},
	 \code{\link{plot.survfit}}, \code{\link{quantile.survfit}},
	 \code{\link{residuals.survfit}}, \code{\link{summary.survfit}}
       }

\keyword{ survival}

