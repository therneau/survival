\name{survfit.object}
\alias{survfit.object}
\alias{survfitms.object}
\title{
  Survival Curve Object 
}
\description{
This class of objects is returned by the \code{survfit} class of functions
to represent a fitted survival curve.
For a multi-state model the object has class \code{c('survfitms', 'survfit')}.

Objects of this class have methods for the functions \code{print},
\code{summary}, \code{plot}, \code{points} and \code{lines}. The
\code{\link{print.survfit}} method does more computation than is typical
for a print method and is documented on a separate page.
}
\section{Structure}{
The following components must be included in a legitimate 
\code{survfit} or \code{survfitms} object. 
}
\arguments{
  \item{n}{
    total number of observations in each curve.
  }
  \item{time}{
    the time points at which the curve has a step. 
  }
  \item{n.risk}{
    the number of subjects at risk at t. 
  }
  \item{n.event}{
    the number of events that occur at time t. 
  }
  \item{n.enter}{
    for counting process data only, and only if there was an \code{id}
    argument, the number of subjects that enter the risk set during the
    current interval.  If there are event/censoring times at 1, 3, 5 for
    instance, someone who enters at time 1 is counted in the (1, 3]
    interval, i.e., appears in the row for time 3. 
  }
  \item{n.censor}{
    for counting process data only,
    the number of subjects who exit the risk set,
    without an event,  at time t. 
    (For right censored data, this number can be computed from the successive
    values of the number at risk).
  }
  \item{surv}{
    the estimate of survival at time t+0. 
    This may be a vector or a matrix. The latter occurs when a set of
    survival curves is created from a single Cox model, in which case
    there is one column for each covariate set. 
  }
  \item{pstate}{
    a multi-state survival will have the \code{pstate} component
    instead of \code{surv}.
    It will be a matrix containing the estimated probability
    of each state at each time, one column per state.
  }
  \item{std.err}{
    for a survival curve this contains standard error of the cumulative
    hazard or -log(survival), for a multi-state curve it contains the
    standard error of prev.  This difference is a reflection of
    the fact that each is the natural calculation for that case.
  }
  \item{cumhaz}{optional.  Contains the cumulative
    hazard for each possible transition.
  }
  \item{counts}{optional. If weights were used, the \code{n.risk} etc
    elements contain weighted sums; the \code{counts} matrix will
    contain unweighted values.  Weighted values are normally more useful
    for further computation, unweighted may be preferred for labeling or
    printout.
  }
  \item{strata}{
    if there are multiple curves, this component gives the number of
    elements of the \code{time}  vector corresponding to the first curve,
    the second curve, and so on.
    The names of the elements are labels for the curves. 
  }
  \item{upper}{optional
    upper (2-sided) confidence limit for the survival curve or probability in state
  }
  \item{lower}{optional 
    lower (2-sidedconfidence limit for the survival curve or probability in state
  }
  \item{t0}{optional, the starting time for the curve}
  \item{p0, sp0}{for a multistate object, the distribution of starting
    states.  If the curve has a strata dimension, this will be a matrix
    one row per stratum.  The \code{sp0} element has the standard error
    of p0, if p0 was estimated.
  }
  \item{newdata}{for survival curves from a fitted model, this contains
    the covariate values for the curves
  }

  \item{n.id}{the total number of unique id values that contributed to
    the curve.  This is only available if the original call used the id
    option.
  }
  \item{conf.type}{
    the approximation used to compute the confidence limits. 
  }
  \item{conf.int}{
    the level of the confidence limits, e.g. 90 or 95\%. 
  }
  \item{transitions}{for multi-state data, the total number
    of transitions of each type.}
  \item{na.action}{
    the returned value from the na.action function, if any.  It will be used 
    in the printout of the curve, e.g., the number of observations deleted due 
    to missing values. 
  }
  \item{call}{
    an image of the call that produced the object. 
  }
  \item{type}{
    type of survival censoring. 
  }
  \item{influence.p, influence.c}{optional influence
    matrices for the \code{pstate} (or \code{surv}) and for the
    \code{cumhaz} estimates.
    A list with one element per stratum, each
    element of the list is an array indexed by subject, time, state.
  }
  \item{version}{the version of the object.  Will be missing, 2, or 3}
}

\section{Subscripts}{
Survfit objects can be subscripted. 
This is often used to plot a subset of the curves, for instance.
From the user's point of view the \code{survfit} object appears to be
a vector, matrix, or array of curves.
The first dimension is always the underlying number of curves or
``strata'';
for multi-state models the state is always the last dimension.
Predicted curves from a Cox model can have a second dimension
which is the number of different covariate prediction vectors.
}

\section{Details}{
The \code{survfit} object has evolved over time: when first created
there was no thought of multi-state models for instance.  This evolution
has almost entirely been accomplished by the addition of new elements.

For both plots of the curves and computation of the restricted mean
time in state (RMTS) we need the concept of a starting point t0 and
starting prevalence of the states p0 for each curve.
(Sojourn time, area under the curve and restricted mean survival time
are other labels for the RMTS).
Time 0 is not, by default, included as part of the standard
tableau of results, i.e., time, number at risk, number of events, etc.
For simple survival with a 0/1 status variable, the starting state p0 is
the obvious default of "everyone alive", and t0 is formally not discernable
from the data and so was left out. 
(A design decision made in 1986, and now far too late to change.)
However, for plots t0 is assumed to be the minimum of 0 and all observed
times. Negative survival times are unusual but not invalid.
Multi-state survival curves include \code{t0} and \code{p0} as a part of
the returned object. The first is a single value for all curves, the
second is per curve.  

The \code{survfit0} routine can be used to add these values to the main
curve data, this is done by the default print, plot, and summary methods
for survfit objects.  The methods vignette has discussion of the
rationale of how t0 and p0 are chosen in the multi-state case.
Notice that if there is an event at time t0, e.g., a death on day 0 for
competing risks, then p0 will contain the prevalence just before that
event occured.  
}
\seealso{
\code{\link{plot.survfit}},  
\code{\link{summary.survfit}},
\code{\link{print.survfit}},
\code{\link{survfit}},
\code{\link{survfit0}}
}
\keyword{survival}
