\name{anova.coxph}
\alias{anova.coxph}
\alias{anova.coxphlist}
\title{Analysis of Deviance for a Cox model.}
\usage{
\method{anova}{coxph}(object, \dots,  test = 'Chisq')
}
\description{
  Compute an analysis of deviance table for one or more Cox model fits.
}
\arguments{
  \item{object}{An object of class \code{coxph}}
  \item{\dots}{Further \code{coxph} objects}
  \item{test}{a character string. The appropriate test is a chisquare, all
other choices result in no test being done.}
}
\details{
  Specifying a single object gives a sequential analysis of deviance
  table for that fit.  That is, the reductions in the model log-likelihood
  as each term of the formula is added in turn are given in as
  the rows of a table, plus the log-likelihoods themselves.
  A robust variance estimate is normally used in situations where the
  model may be mis-specified, e.g., multiple events per subject.
  In this case a comparison of partial-likelihood values does not make
  sense, and \code{anova} will refuse to print results.
  
  If more than one object is specified, the table has a row for the
  degrees of freedom and loglikelihood for each model. For all
  but the first model, the change in degrees of freedom and loglik
  is also given. (This only make statistical sense if the models are
  nested.)  It is conventional to list the models from smallest to
  largest, but this is up to the user.
  
  The table will optionally contain test statistics (and P values)
  comparing the reduction in loglik for each row.
}
\value{
  An object of class \code{"anova"} inheriting from class \code{"data.frame"}.
}
\section{Warning}{
  The comparison between two or more models by \code{anova} or
   will only be valid if they
  are fitted to the same dataset. This may be a problem if there are
  missing values.}
\seealso{
  \code{\link{coxph}}, \code{\link{anova}}.
}
\examples{
fit <- coxph(Surv(futime, fustat) ~ resid.ds *rx + ecog.ps, data = ovarian) 
anova(fit)
fit2 <- coxph(Surv(futime, fustat) ~ resid.ds +rx + ecog.ps, data=ovarian)
anova(fit2,fit)
 }
\keyword{models}
\keyword{regression}
\keyword{survival}

