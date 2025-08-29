\name{survfit_confint}
\alias{survfit_confint}
\title{
Confidence intervals for survival curves
}
\description{
Compute the confidence intervals for a survfit object, using one of
several approximations.
}
\usage{
survfit_confint(p, se, logse = TRUE, conf.type, conf.int = 0.95, selow, ulimit = TRUE)
}
\arguments{
  \item{p}{the estimate, i.e., the
    \code{surv}, \code{pstate} or \code{cumhaz} component of a survfit object
}
  \item{se}{vector or matrix of standard errors}
  \item{logse}{if TRUE, the provided \code{se} component is the std of
    log(p)}
  \item{conf.type}{
    one of "plain", "log", "log-log", "logit", or "arcsin"}
  \item{conf.int}{the confidence level, a value strictly between 0 and 1}
  \item{selow}{alternate standard error for the lower limit}
  \item{ulimit}{if TRUE, constrain the result to be between 0 and 1}
}
\details{
  This routine is used internally to create confidence intervals and would
  rarely be called by a user (but perhaps by another package).
  As such it does very little checking of the input arguments.  

  "Plain" intervals of p +- 1.96*se have been shown by many authors to
  perform poorly.  Which of the other choices is "best" depends on the
  details of the simulation study used to evaluate the method, all do
  quite well overall.

  For log intervals, the CI at p=0 is exp(log(0) + se) which becomes
  NA in R, and likewise for p= 0 or 1 when using log-log or logit.
  If the se is 0, however, the routine returns p as the upper and lower
  limit.
  (One side effect of the NA is that the plotted se bands for a
  Kaplan-Meier whose final time point is a death will stop at the death
  time and not plunge to 0; i.e., not have the final stair step, 
  a graphical aspect we prefer.)
}
\value{a list with components \code{lower} and \code{upper}.}
\seealso{\code{\link{survfit.object}}}
\keyword{survival}
