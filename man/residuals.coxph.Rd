\name{residuals.coxph}
\alias{residuals.coxph.penal}
\alias{residuals.coxph.null}
\alias{residuals.coxph}
\alias{residuals.coxphms}
\title{Calculate Residuals for a `coxph' Fit
}
\description{
Calculates martingale, deviance, score or Schoenfeld residuals for a
Cox proportional hazards model.
}
\usage{
\method{residuals}{coxph}(object,
       type=c("martingale", "deviance", "score", "schoenfeld",
	      "dfbeta", "dfbetas", "scaledsch","partial"),
       collapse=FALSE, weighted= (type \%in\% c("dfbeta", "dfbetas")), 
       na.action, ...)
\method{residuals}{coxphms}(object,
       type=c("martingale", "score", "schoenfeld",
	      "dfbeta", "dfbetas", "scaledsch"),
       collapse=FALSE, weighted= FALSE, na.action, ...)
\method{residuals}{coxph.null}(object,
       type=c("martingale", "deviance","score","schoenfeld"),
       collapse=FALSE, weighted= FALSE, ...)
}
\arguments{
\item{object}{
an object inheriting from class \code{coxph}, representing a fitted Cox
regression model.
Typically this is the output from the \code{coxph} function.
}
\item{type}{
character string indicating the type of residual desired.
Possible values are \code{"martingale"}, \code{"deviance"}, \code{"score"}, \code{"schoenfeld"},
"dfbeta"', \code{"dfbetas"}, \code{"scaledsch"} and \code{"partial"}.
Only enough of the string to determine a unique match is required.
}
\item{collapse}{
  either a logical (TRUE/FALSE), or a vector indicating which rows to
  collapse (sum) over.
  In time-dependent models more than one row data can pertain
  to a single individual or group.

 If the original coxph fit contained a \code{cluster} argument this will
 be used as the collapsing index, or lacking this the \code{id} argument
 will be used. The third alternative is an explicit vector.
If there were 4 individuals represented by 3, 1, 2 and 4 rows of data
respectively, then \code{collapse=c(1,1,1, 2, 3,3, 4,4,4,4)} could be used to
obtain per subject rather than per observation residuals.
}
\item{weighted}{
if \code{TRUE} and the model was fit with case weights, then the weighted
residuals are returned.
}
\item{na.action}{either \code{"na.omit"} or \code{"na.exclude"},
  the default is the
  global option in effect when the \code{coxph} model was fit.
  The latter option causes the residuals to be the same size as the
  original data \emph{before} missings were excluded by \code{coxph}, it
  is useful when plotting residuals versus a variable in the data
  frame, for instance.  Not applicable to Schoenfeld or scaled
  Schoenfeld residuals, which have a row per event, nor when
  \code{collapse=TRUE}.
  }
  \item{...}{other unused arguments}
}

\value{a vector or matrix of residuals}

\section{Details}{

The martingale residuals are a vector for each observation, for a
multi-state model they will be a matrix with one row per observation and
one column per transition. If a particular observation was not at risk
for a given transition, that matrix value will be zero.
Deviance residuals are a transform of the martingale residuals, one
which has not proven to be useful. They have not been defined for
multistate models.

The score residuals will be a matrix with one row per observation and
one column per covariate, for a multistate model the result is an array
with dimensions of observation, covariate, and transition.
This constains each observation's contribution to the score vector
(first derivative of the log partial likelihood).
Two transformations of
this are often more useful: \code{dfbeta} is the approximate change in the
coefficient vector if that observation were dropped,
and \code{dfbetas} is the approximate change in the coefficients, scaled by
the standard error for the coefficients.
For multi-state models the result will be an array with dimensions of
(observation, covariate, transition). Some variables might not be used for all
transitions, on which case the appropriate element of the array will be zero.

For each of the above first dimension is in the same order as the original data,
either before (\code{na.action= "na.exclude"}) or after removal of
missings by the coxph function.
When \code{collapse = TRUE} the result will be in the order in which the
cluster and/or id variables are encountered in the data.

The Schoenfeld residuals have a row for each event time and a column for
each covariate. If the coxph model was stratified, the resulting matrix
will contain a \code{strata} attribute showing which stratum contributed
each of the events, and the result will be sorted by time within strata.
For multistate models the resulting matrix will have the events for
transition 1, then transition 2, etc. and will also have a
\code{transition} attribute.  Different transitions will have different
event times, so an array form of the output is not practical.
}
\references{
T. Therneau, P. Grambsch, and T. Fleming. "Martingale based residuals
for survival models", \emph{Biometrika}, March 1990.
}
\seealso{
  \code{\link{coxph}}, \code{\link{predict.coxph}}
}
\examples{

 fit <- coxph(Surv(start, stop, event) ~ (age + surgery)* transplant,
               data=heart)
 mresid <- resid(fit, collapse=heart$id)
}
\keyword{survival}
