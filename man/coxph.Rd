\name{coxph}
\alias{coxph}
\title{
  Fit Proportional Hazards Regression Model 
}
\description{
Fits a Cox proportional hazards regression model. 
Time dependent variables, time dependent strata, multiple events per subject, 
and other extensions are incorporated using the counting process formulation 
of Andersen and Gill. 
}
\usage{
coxph(formula, data=, weights, subset, 
      na.action, init, control, 
      ties=c("efron","breslow","exact"), 
      singular.ok=TRUE, robust, 
      model=FALSE, x=FALSE, y=TRUE, tt, method=ties,
      id, cluster, istate, statedata, nocenter=c(-1, 0, 1), ...)
}
\arguments{
  \item{formula}{
    a formula object, with the response on the left of a \code{~} operator, and 
    the terms on the right.  The response must be a survival object as 
    returned by the \code{Surv} function.  For a multi-state model the
    formula may be a list of formulas.
  }
  \item{data}{
    a data.frame in which to interpret the variables named in 
    the \code{formula}, or in the \code{subset} and the \code{weights}
    argument. 
  }
  \item{weights}{
    vector of case weights, see the note below.
    For a thorough discussion of these see the
    book by Therneau and Grambsch.
  }
  \item{subset}{
    expression indicating which subset of the rows of data should be used in 
    the fit.    All observations are included by default. 
  }
  \item{na.action}{
    a missing-data filter function.  This is applied to the model.frame
    after any 
    subset argument has been used.  Default is \code{options()\$na.action}. 
  }
  \item{init}{
    vector of initial values of the iteration.  Default initial 
    value is zero for all variables. 
  }
  \item{control}{
    Object of class \code{\link{coxph.control}} specifying iteration limit
    and other control options. Default is \code{coxph.control(...)}.
  }
  \item{ties}{
    a character string specifying the method for tie handling.  If there  
    are no tied death times all the methods are equivalent.
    The Efron approximation is used as the default, it is more 
    accurate when dealing with tied death times, and is as efficient 
    computationally. (But see below for multi-state models.)
    The ``exact partial likelihood'' is 
    equivalent to a conditional logistic model, and is appropriate when
    the times are a small set of discrete values.
  }
  \item{singular.ok}{
    logical value indicating how to handle collinearity in the model matrix. 
    If \code{TRUE}, the program will automatically skip over columns of the X 
    matrix that are linear combinations of earlier columns.  In this case the 
    coefficients for such columns will be NA, and the variance matrix will 
    contain zeros. For ancillary calculations, such as the linear
    predictor, 
    the missing coefficients are treated as zeros. 
  }
  \item{robust}{should a robust variance be computed.
    The default is TRUE if: there is a \code{cluster} argument, there
    are case weights that are not 0 or 1, or there are \code{id} values
    with more than one event. 
  }

  \item{id}{optional variable name that identifies subjects.  Only
    necessary when a subject can have multiple rows in the data, and
    there is more than one event type.  This variable will normally be
    found in \code{data}.}
  \item{cluster}{optional variable which clusters the observations, for
    the purposes of a robust variance.  If present, it implies
    \code{robust}.
    This variable will normally be found in \code{data}.}
  \item{istate}{optional variable giving the current state at the start
    each interval. This variable will normally be found in \code{data}.}
  \item{statedata}{optional data set used to describe multistate models.}
  \item{model}{
    logical value: if \code{TRUE}, the model frame is returned in component
    \code{model}. 
  }
  \item{x}{
    logical value: if \code{TRUE}, the x matrix is returned in
    component \code{x}. 
  }
  \item{y}{
    logical value: if \code{TRUE}, the response vector is returned in
    component \code{y}. 
  }
  \item{tt}{optional list of time-transform functions.}
  \item{method}{alternate name for the \code{ties} argument.}
  \item{nocenter}{columns of the X matrix whose values lie strictly
    within this set are not recentered.  Remember that a factor
    variable becomes a set of 0/1 columns.}
  \item{...}{Other arguments will be passed to \code{\link{coxph.control}}
  }
}
\value{
an object of class \code{coxph} representing the fit. 
See \code{coxph.object} and \code{coxphms.object} for details.  
}
\section{Side Effects}{
  Depending on the call, the \code{predict}, \code{residuals},
and \code{survfit} routines may 
need to reconstruct the x matrix created by \code{coxph}.
It is possible for this to fail, as in the example below in
which the predict function is unable to find \code{tform}.

\preformatted{  tfun <- function(tform) coxph(tform, data=lung)
  fit <- tfun(Surv(time, status) ~ age)
  predict(fit)}

In such a case add the \code{model=TRUE} option to the
\code{coxph} call to obviate the
need for reconstruction, at the expense of a larger \code{fit}
object.  
}
\details{
The proportional hazards model is usually expressed in terms of a 
single survival time value for each person, with possible censoring. 
Andersen and Gill reformulated the same problem as a counting process; 
as time marches onward we observe the events for a subject, rather 
like watching a Geiger counter. 
The data for a subject is presented as multiple rows or "observations", 
each 
of which applies to an interval of observation (start, stop].

The routine internally scales and centers data to avoid overflow in
the argument to the exponential function.  These actions do not change
the result, but lead to more numerical stability.
Any column of the X matrix whose values lie within \code{nocenter} list
are not recentered.  The practical consequence of the default is to not
recenter dummy variables corresponding to factors.
However, arguments to offset are not scaled since there are situations
where a large offset value is a purposefully used.
In general, however, users should not avoid very large numeric values
for an offset due to possible loss of precision in the estimates.
}
\section{Case weights}{
Case weights are treated as replication weights, i.e., a case weight of
2 is equivalent to having 2 copies of that subject's observation.
When computers were much smaller grouping like subjects together was a
common trick to used to conserve memory.  Setting all weights to 2 for
instance will give the same coefficient estimate but halve the variance.
When the Efron approximation for ties (default) is employed replication
of the data will not give exactly the same coefficients as the weights option,
and in this case the weighted fit is arguably the correct one.

When the model includes a \code{cluster} term or the \code{robust=TRUE}
option the computed variance is based on an infinitesimal jackknife,
which treats any weights as sampling weights;
setting all weights to 2 will in this case give the same variance as
weights of 1.  If there are any non-integer weights, the code also
defaults \code{robust=TRUE}.
}
\section{Special terms}{
  There are a few special terms that may be used in the model equation:
\code{strata}, \code{tt}, \code{pspline}, \code{frailty} and
\code{ridge}.
Each look like an ordinary function, e.g. \code{+ strata(group)} but are
specially identifies so that they can be treated in a special way.
The term \code{+ cluster(group)} is also but is depricated, use a
a \code{cluster} arguement outside the formula instead.

A \code{strata} term identifies a stratified Cox model; separate baseline 
hazard functions are fit for each strata.

A time-transform term allows variables to vary dynamically in time.  In
this case the \code{tt} argument will be a function or a list of
functions (if there are more than one tt() term in the model) giving the
appropriate transform.   See the examples below.
If the \code{id} variable is not 
unique, it is assumed that it identifies clusters of correlated
observations.

A time-transform term allows variables to vary dynamically in time.  In
this case the \code{tt} argument will be a function or a list of
functions (if there are more than one tt() term in the model) giving the
appropriate transform.   See the examples below.

One user mistake that has recently arisen is to slavishly follow the
advice of some coding guides and prepend \code{survival::} onto
everthing, including the special terms, e.g.,
  \code{survival::coxph(survival:Surv(time, status) ~ age +
    survival::strata(inst), data=lung)}
  
First, for this actually will not fit the model that was intended, as
the :: interferes with the recognition of specials by the underlying
\code{model.frame} function; there will a coefficient per
institution rather than fitting a stratified model.
A similar issue arises from using \code{stats::offset} as a term in
a glm model.

From survival 3.8-1 onward it is also unnessary: the common formula
arguments that are part of the survival namespace will be found in that
namespace, i.e., Surv, strata, tt, pspline, and cluster.  A survival::
prefex found on any of these is actually removed before evaluation of
the formula.
This only affects the formula itself; whether to use the qualified
form \code{survival::coxph} for the call itself is a different discussion.
}

\section{Robust variance}{
The robust estimate arises from many different arguments and thus has
had many labels.  It is variously known as the
Huber sandwich estimator, White's estimate (linear models/econometrics),
the Horvitz-Thompson estimate (survey sampling), the working
independence variance (generalized estimating equations), the
infinitesimal jackknife, and the Wei, Lin, Weissfeld (WLW) estimate.

If there is an \code{id} or \code{cluster} argument in the call, or
\code{robust=TRUE}, or there are non-integer case weights, the robust
variance is the default.  If there are multiple rows per id
and no id has more than one event and the underlying proportional
hazards model is exactly correct, then the multiple rows for a given
id will represent independent increments and the variance (robust or standard)
obtained without an \code{id} statement will be an unbiased estimate.
But it will not be equal to the estimate with \code{id} specified and we
do not recommend it.
}

\section{Convergence}{
In certain data cases the actual MLE estimate of a 
coefficient is infinity, e.g., a dichotomous variable where one of the 
groups has no events.  When this happens the associated coefficient 
grows at a steady pace and a race condition will exist in the fitting 
routine: either the log likelihood converges, the information matrix 
becomes effectively singular, an argument to exp becomes too large for 
the computer hardware, or the maximum number of interactions is
exceeded.
(Most often number 1 is the first to occur.)
The routine attempts to detect when this has happened,
not always successfully.
The primary consequence for the user is that the Wald statistic =
coefficient/se(coefficient) is not valid in this case and should be
ignored; the likelihood ratio and score tests remain valid however.
}
\section{Ties}{
There are three possible choices for handling tied event times.
The Breslow approximation is the easiest to program and hence became the
first option coded for almost all computer routines. It then ended up
as the default option when other options were added in order to "maintain
backwards compatability".  The Efron option is more accurate if there are
a large number of ties, and it is the default option here.
In practice the number of ties is usually small, in which case all the
methods are statistically indistinguishable.

Using the "exact partial likelihood" approach the Cox partial likelihood
is equivalent to that for matched logistic regression.  (The
\code{clogit} function uses the \code{coxph} code to do the fit.)
It is technically appropriate when the time scale is discrete and has
only a few unique values, and some packages refer to this as the
"discrete" option.  There is also an "exact marginal likelihood" due to
Prentice which is not implemented here.

The calculation of the exact partial likelihood is numerically intense.
Say for instance 180 subjects are at risk on day 7 of which 15 had an event;
then the code needs to compute sums over all 180-choose-15 > 10^43 different
possible subsets of size 15.  There is an efficient recursive algorithm
for this task, but even with this the computation can be insufferably
long.  With (start, stop) data it is much worse since the recursion
needs to start anew for each unique start time.

Multi state models are a more difficult case.
First of all, a proper extension of the Efron argument is much more
difficult to do, and this author is not yet fully convinced that the
resulting algorithm is defensible.  
Secondly, the current code for Efron case does not consistently
compute that extended logic (and extension would require major changes
in the code).  Due to this
complexity, the default is \code{ties='breslow'} for the multistate
case. If \code{ties='efron'} is selected the current code will, in
effect, only apply to to tied transitions of the same type.

A separate issue is that of artificial ties due to floating-point
imprecision. See the vignette on this topic for a full explanation or
the \code{timefix} option in \code{coxph.control}.
Users may need to add \code{timefix=FALSE} for simulated data sets.
}
\section{Penalized regression}{
  \code{coxph} can maximise a penalised partial likelihood with
arbitrary user-defined penalty.  Supplied penalty functions include
ridge regression (\link{ridge}), smoothing splines
(\link{pspline}), and frailty models (\link{frailty}).
} 
\references{
Andersen, P. and Gill, R. (1982). 
Cox's regression model for counting processes, a large sample study. 
\emph{Annals of Statistics}
\bold{10}, 1100-1120. 

Therneau, T., Grambsch, P., Modeling Survival Data: Extending the Cox Model. 
Springer-Verlag, 2000.
}
\seealso{
  \code{\link{coxph.object}}, \code{\link{coxphms.object}},
  \code{\link{coxph.control}},
  \code{\link{cluster}},  \code{\link{strata}},  \code{\link{Surv}},
  \code{\link{survfit}}, \code{\link{pspline}}. 
}
\examples{
# Create the simplest test data set 
test1 <- list(time=c(4,3,1,1,2,2,3), 
              status=c(1,1,1,0,1,1,0), 
              x=c(0,2,1,1,1,0,0), 
              sex=c(0,0,0,0,1,1,1)) 
# Fit a stratified model 
coxph(Surv(time, status) ~ x + strata(sex), test1) 
# Create a simple data set for a time-dependent model 
test2 <- list(start=c(1,2,5,2,1,7,3,4,8,8), 
              stop=c(2,3,6,7,8,9,9,9,14,17), 
              event=c(1,1,1,1,1,1,1,0,0,0), 
              x=c(1,0,0,1,0,1,1,1,0,0)) 
summary(coxph(Surv(start, stop, event) ~ x, test2)) 

#
# Create a simple data set for a time-dependent model
#
test2 <- list(start=c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
                stop =c(2, 3, 6, 7, 8, 9, 9, 9,14,17),
                event=c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
                x    =c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0) )


summary( coxph( Surv(start, stop, event) ~ x, test2))

# Fit a stratified model, clustered on patients 

bladder1 <- bladder[bladder$enum < 5, ] 
coxph(Surv(stop, event) ~ (rx + size + number) * strata(enum),
      cluster = id, bladder1)

# Fit a time transform model using current age
coxph(Surv(time, status) ~ ph.ecog + tt(age), data=lung,
     tt=function(x,t,...) pspline(x + t/365.25))
}
\keyword{survival}

