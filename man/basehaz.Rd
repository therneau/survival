\name{basehaz}
\alias{basehaz}
\title{Alias for the survfit function}
\description{
  Compute the predicted survival curve for a Cox model.
}
\usage{
basehaz(fit, newdata, centered=TRUE)
}
\arguments{
  \item{fit}{a coxph fit}
  \item{newdata}{a data frame containing one row for each predicted
    survival curve, said row contains the covariate values for that curve}
  \item{centered}{ignored if the \code{newdata} argument is present.
    Otherwise, if TRUE return data from a predicted survival curve
    for the covariate values \code{fit$mean}, if FALSE return a
    prediction for all covariates equal to zero.}
}
\details{
  This function is an alias for \code{survfit.coxph}, which does the
  actual work and has a richer set of options.
  Look at that help file for more discussion and explanation.
  This alias exists primarily because some users look for predicted survival
  estimates under this name. 

  The function returns a data frame containing the \code{time},
  \code{cumhaz} and optionally the strata (if the fitted Cox model used
  a strata statement), which are copied from the \code{survfit} result.

  If H(t; z) is the predicted cumulative hazard for an observation with
  covariate vector z, then H(t;x) = H(t;z) r(x,z)
  where r(x,z)= exp(beta[1](x[1]- z[1]) + beta[2](x[2]-z[2]) + ...) =
  \code{exp(sum(coef(fit) * (x-z)))} is the
  Cox model's hazard ratio for covariate vector x vs covariate vector z.
  That is, 
  the cumulative hazard H for a single reference value z is sufficient to
  provide the hazard for any covariate values.
  The predicted survival curve is S(t; x)= exp(-H(t;x)).
  There is not a simple transformation for the variance of H, however.

  Many textbooks refer to H(t; 0) as "the" baseline hazard
  for a Cox model; this is returned by the \code{centered= FALSE}
  option.
  However, due to potential overflow or underflow in the exp() function
  this can be a very bad idea in practice. The authors do not recommend
  this option, but for users who insist: caveat emptor.
  Offset terms can pose a particular challenge for the underlying code
  and are always recentered; to override this use the newdata argument
  and include the offset as one of the variables.
}
\value{
  a data frame with variable names of \code{hazard}, \code{time} and
  optionally \code{strata}.  The first is actually the cumulative hazard.
}

\seealso{\code{\link{survfit.coxph}}}
\keyword{survival }

