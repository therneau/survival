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
    total number of subjects in each curve.
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
    for counting process data only,
    the number of subjects that enter at time t. 
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
  \item{cumhaz hazard}{optional.  Contains the cumulative
    hazard for each possible transtion.
  }
  \item{strata}{
    if there are multiple curves, this component gives the number of
    elements of the \code{time}  vector corresponding to the first curve,
    the second curve, and so on.
    The names of the elements are labels for the curves. 
  }
  \item{upper}{optional
    upper confidence limit for the survival curve or pstate
  }
  \item{lower}{options 
    lower confidence limit for the survival curve or pstate
  }
  \item{start.time, p0, sp0}{see \code{details} below
  }
  \item{newdata}{for survival curves from a fitted model, this contains
    the covariate values for the curves
  }

  \item{n.all}{the total number of observations that were available
    For counting process data, and any time that the 
    \code{start.time} argument was used, 
    not all may have been used in creating the curve, in which case this value
    will be larger than \code{n} above.
    The \code{print} and \code{plot} routines in the package do no use this
    value, it is for information only.
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
One change in survival version 3 is the addition of a \code{survfitconf} routine
which will compute confidence intervals for a \code{survfit} object.
This allows the computation of CI intervals to be deferred until later,
if desired, rather than making them a permanent part of the object.
Later iterations of the base routines may omit the confidence intervals.

The survfit object starts at the first observation time, but survival
curves are normally plotted from time 0.
A helper routine \code{survfit23} is used to add this first time point
and align the data.
}
\seealso{
\code{\link{plot.survfit}},  
\code{\link{summary.survfit}},
\code{\link{print.survfit}},
\code{\link{survfit}},
\code{\link{survfit23}}
}
\keyword{survival}
