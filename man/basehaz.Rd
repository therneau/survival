\name{basehaz}
\alias{basehaz}
\title{Alias for the survfit function}
\description{
  Compute the predicted survival curve for a Cox model.
}
\usage{
basehaz(fit, centered=TRUE)
}
\arguments{
  \item{fit}{a coxph fit}
  \item{centered}{if TRUE return data from a predicted survival curve at
    the mean values of the covariates \code{fit$mean}, if FALSE return a
    prediction for all covariates equal to zero.}
}
\details{
  This function is simply an alias for \code{survfit}, which is the
  actual function that does all the computations.  See the manual page
  for that function for the preferred use.  This function survives only
  for backwards support of prior usage.

  The function returns a data frame containing the \code{time},
  \code{cumhaz} and optionally the strata (if the fitted Cox model used
  a strata statement), which are copied the \code{survfit} result.
  If there are factor variables in the model, then the default
  predictions at the "mean" are meaningless since they do not correspond
  to any possible subject; correct results require use of the
  \code{newdata} argument of survfit.  Results for all covariates =0 are
  normally only of use as a building block for further calcuations.
}
\value{
  a data frame with variable names of \code{hazard}, \code{time} and
  optionally \code{strata}.  The first is actually the cumulative hazard.
}

\seealso{\code{\link{survfit.coxph}}}
\keyword{survival }

