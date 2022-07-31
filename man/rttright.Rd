\name{rttright}
\alias{rttright}
\title{Compute redistribute-to-the-right weights}
\description{
 For many survival estimands, one approach is to redistribute each
 censored observation's weight to those other observations with a longer
 survival time (think of distributing an estate to the heirs). Then
 compute on the remaining, uncensored data.
}
\usage{
rttright(formula, data, weights, subset, na.action, times, id, timefix = TRUE)
}
\arguments{
  \item{formula}{
    a formula object, which must have a 
    \code{Surv} object as the  
    response on the left of the \code{~} operator and, if desired, terms  
    separated by + operators on the right. 
    Each unique combination of predictors will define a separate strata.
  }
  \item{data}{
    a data frame in which to interpret the variables named in the formula, 
    \code{subset} and \code{weights} arguments. 
  }
  \item{weights}{
    The weights must be nonnegative and it is strongly recommended that  
    they be strictly positive, since zero weights are ambiguous, compared 
    to use of the \code{subset} argument.
  }
  \item{subset}{
    expression saying that only a subset of the rows of the data 
    should be used in the fit. 
  }
  \item{na.action}{
    a missing-data filter function, applied to the model frame, after any 
    \code{subset} argument has been used. 
    Default is \code{options()$na.action}. 
  }
  \item{times}{a vector of time points, for which to return updated
    weights.  If missing, a time after the largest time in the data is
    assumed.}
  \item{id}{optional: if the data set has multiple rows per subject, a
    a variable containing the subect identifier of each row.}
  \item{timefix}{correct for possible round-off error}
}
\details{
  The \code{formula} argument is treated exactly the same as in the
  \code{survfit} function.
  
  Redistribution is recursive: redistribute the weight of the first
  censored observation to all those with longer time, which may include
  other censored observations.  Then redistribute the next smallest and
  etc. up to the specified \code{time} value.
  After re-distributing the weight for a censored observation to other
  observations that are not censored, ordinary non-censored methods can
  often be applied.  For example, redistribution of the weights,
  followed by computation of the weighted cumulative distribution
  function, reprises the Kaplan-Meier estimator.

  A primary use of this routine is illustration of methods or
  exploration of new methods.  Methods that use RTTR directly, such as
  the Brier score, will often do these compuations internally.

  A covariate on the right hand side of the formula causes
  redistribution to occur within group; a censoring in group 1
  redistributes weights to others in group 1, etc.  This is appropriate
  when the censoring pattern depends upon group.
  
}
\value{a vector or matrix of weights, with one column for each requested time} 
\seealso{
  \code{\link{survfit}}
}
\examples{
afit <- survfit(Surv(time, status) ~1, data=aml)
rwt  <- rttright(Surv(time, status) ~1, data=aml)

# Reproduce a Kaplan-Meier
index <- order(aml$time)
cdf <- cumsum(rwt[index])  # weighted CDF
cdf <- cdf[!duplicated(aml$time[index], fromLast=TRUE)]  # remove duplicate times
cbind(time=afit$time, KM= afit$surv, RTTR= 1-cdf)

# Hormonal patients have a diffent censoring pattern
wt2 <- rttright(Surv(dtime, death) ~ hormon, rotterdam, times= 365*c(3, 5))
dim(wt2)
}
\keyword{ survival }
