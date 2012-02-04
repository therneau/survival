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
clogit(formula, data, method=c("exact", "approximate", "efron", "breslow"),
     na.action=getOption("na.action"), subset=NULL,control=coxph.control())
}
\arguments{
  \item{formula}{Model formula}
  \item{data}{data frame }
  \item{method}{correct or approximate conditonal likelihood?}
  \item{na.action}{missing value handling}
  \item{subset}{subset of records to use }
  \item{control}{control values}
}

\value{
  An object of class \code{"clogit"}, which is a wrapper for a
  \code{"coxph"} object.
}
\author{Thomas Lumley}
\details{
  A stratified Cox model with time=constant, status of 1=case 0=control,
  and using the exact partial likelihood has the same likelihood formula
  as a conditional logistic regression.  The clogit routine creates
  the necessary dummy variable of times (all 1) and then calls coxph.

  The computation of the exact partial likelihood can be very slow,
  however.  If a particular strata had say 10 events out of 20 subjects
  we have to add up a denominator that involves all possible ways of
  choosing 10 out of 20, which is 20!/(10! 10!) = 184756 terms. Gail et
  al describe a fast recursion method, which largely ameleorates 
  this; 
  it was incorporated into version 2.36-11 of the survival package.
  (Most of the time conditional logistic modeling 
  is applied data with 1 case + k controls per set, 
  which involves choosing 1 out of k and the computational issue
  above does not arise.)  
  The 'appoximate' option maps to the
  Breslow approximation for historical reasons.
  }
\seealso{\code{\link{strata}},\code{\link{coxph}},\code{\link{glm}} }

\examples{
\dontrun{clogit(case ~ spontaneous + induced + strata(stratum), data=infert)}

# A multinomial response recoded to use clogit
#  The data set has one copy per possible outcome level, with new
#  variable tocc = target occupation for this copy, and "case"
#  whether that is the actual outcome for each subject.
# See the catspec package for more details.
resp <- levels(logan$occupation)
n <- nrow(logan)
indx <- rep(1:n, length(nresp))
logan2 <- data.frame(logan[indx,],
                     id = indx,
                     tocc = factor(resp, each=n))
logan2$case <- (logan2$occupation == logan2$tocc)
clogit(case ~ tocc + tocc:education + strata(id), logan2)
}
\section{References}{
  Michell H Gail, Jay H Lubin and Lawrence V Rubinstein.  Likelihood
calculations for matched case-control studies and survival studies with
tied death times.  Biometrika 68:703-707, 1980.
}
\keyword{survival}%-- one or more ...
\keyword{models}
