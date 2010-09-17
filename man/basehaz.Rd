\name{basehaz}
\alias{basehaz}
\title{Compute the baseline survival curve for a Cox model}
\description{
  Compute the baseline survival curve for a Cox model.
}
\usage{
basehaz(fit, centered = TRUE)
}
\arguments{
  \item{fit}{The result of a \code{coxph} fit.}
  \item{centered}{If TRUE, the resultant curve is for a hypothetical
    subject whose covariate values are the corresponding means from the
    original data, otherwise for a hypothetical subject with a mean
    vector of zero.}
}
\details{
  This function exists primarily because users will look for the phrase
  'baseline hazard' (often SAS converts looking for familiar
  keywords.) The primary function for creating a survival curve is
  \code{survfit}, which this calls.  See that manual page for more
  options, including confidence limits and the ability to use other
  covariate vectors.  The result of \code{survfit} also has print, plot
  and summary methods that make it far more useful. 
}
\value{
  a data frame with components
  \item{time}{The sorted vector of unique time points (those at which an
    event occurred}
  \item{hazard}{The baseline hazard function}
  \item{strata}{If \code{fit} was a stratified Cox model, the strata.
    There will be one survival curve per strata.}
}

\seealso{\code{\link{survfit}}}
\keyword{survival }

