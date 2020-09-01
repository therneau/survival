\name{brier}
\alias{brier}
\title{Calculate the Brier score}
\description{ The Brier score is an R^2-like statistic focused on
  the event variable, from Surv(time, event) survival data.
}
\usage{
brier(fit, times, newdata, censformula)
}
%- maybe also 'usage' for other objects documented here.
\arguments{the result of a \code{coxph} fit
  \item{fit}{the result of a \code{coxph} call}
  \item{times}{the set of time points at which to do prediction.
    The default is to use all event times in the data.}
  \item{newdata}{optional validation data set. If present, predictions
    from the model are used to calculate the Brier score on new data.}
  \item{censformula}{optional formula to model the censoring
    distribution.  The default is no predictors.}
}
\details{
  This returns the Brier score under both the NULL and fitted model.
  }
\value{ a list with elements
  \item{brier}{the overall score}
  \item{b0, b1}{per observation components of the score, under H0 and
 H1}
  \item{times}{time points of assessment}
  \item{other}
}
  \references{
    H Van Houwelingn and H Putter, Dynamic Prediction in Clinical Survival
    Analysis, CRC Press, 2012.
}

\seealso{
\code{\link{rttright}}, ~~~
}
\examples{
}
\keyword{ survival }
