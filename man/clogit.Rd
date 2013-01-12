\name{clogit}
\alias{clogit}
\title{Conditional logistic regression }
\description{
  Estimates a logistic regression model by maximising the conditional
  likelihood. Uses a model formula of the form
  \code{case.status~exposure+strata(matched.set)}.
  The default is to use the exact conditional likelihood, a commonly
  used approximate conditional likelihood is provided for compatibility
  with older software.
}
\usage{
clogit(formula, data, weights, subset, na.action,
 method=c("exact", "approximate", "efron", "breslow"),
 \dots)
}
\arguments{
  \item{formula}{Model formula}
  \item{data}{data frame }
  \item{weights}{optional, names the variable containing case weights}
  \item{subset}{optional, subset the data}
  \item{na.action}{optional na.action argument.  By default the
    global option \code{na.action} is used.}
  \item{method}{use the correct (exact) calculation in the conditional
    likelihood or one of the approximations}
  \item{\dots}{optional arguments, which will be passed to
  \code{coxph.control}}
}

\value{
  An object of class \code{"clogit"}, which is a wrapper for a
  \code{"coxph"} object.
}
\author{Thomas Lumley}
\details{
  It turns out that the logliklihood for a conditional logistic
  regresson model = loglik from a Cox model with a particular data
  structure.  Proving this is a nice homework exercise for a PhD
  statistics class; not too hard, but the fact that it is true is
  surprising.

  When a well tested Cox model routine is available many packages use
  this `trick' rather than writing a new software routine from
  scratch, and this is what the clogit routine does.
  In detail,  a stratified Cox model with each case/control group
  assigned to its own stratum, time set to a constant,
  status of 1=case 0=control,
  and using the exact partial likelihood has the same likelihood formula
  as a conditional logistic regression.  The clogit routine creates
  the necessary dummy variable of times (all 1) and the strata,
  then calls coxph.

  The computation of the exact partial likelihood can be very slow,
  however.  If a particular strata had say 10 events out of 20 subjects
  we have to add up a denominator that involves all possible ways of
  choosing 10 out of 20, which is 20!/(10! 10!) = 184756 terms. Gail et
  al describe a fast recursion method which largely ameleorates 
  this; it was incorporated into version 2.36-11 of the survival package.
  Most of the time conditional logistic modeling 
  is applied data with 1 case + k controls per set, however, 
  where the above the computational issue above does not arise.
  Thus most users will not notice the change but for others
  computation time will drop precipitously.
  The 'appoximate' option maps to the
  Breslow approximation for the Cox model, for historical reasons.

  It is not clear how case weights should be handled.  For instance if
  there are two deaths in a strata, one with weight=1 and one with
  weight=2, should the likelihood calculation consider all subsets of
  size 2 or all subsets of size 3?
  Consequently, case weights are ignored by the routine.
  }
\seealso{\code{\link{strata}},\code{\link{coxph}},\code{\link{glm}} }

\examples{
\dontrun{clogit(case ~ spontaneous + induced + strata(stratum), data=infert)}

# A multinomial response recoded to use clogit
#  The revised data set has one copy per possible outcome level, with new
#  variable tocc = target occupation for this copy, and case = whether
#  that is the actual outcome for each subject.
# See the catspec package for more details on the Logan approach.
resp <- levels(logan$occupation)
n <- nrow(logan)
indx <- rep(1:n, length(resp))
logan2 <- data.frame(logan[indx,],
                     id = indx,
                     tocc = factor(rep(resp, each=n)))
logan2$case <- (logan2$occupation == logan2$tocc)
clogit(case ~ tocc + tocc:education + strata(id), logan2)
}
\section{References}{
  Michell H Gail, Jay H Lubin and Lawrence V Rubinstein.  Likelihood
calculations for matched case-control studies and survival studies with
tied death times.  Biometrika 68:703-707, 1980.
}
\keyword{survival}
\keyword{models}
