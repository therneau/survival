\name{survfit.coxph}
\alias{survfit.coxph}
\alias{survfit.coxphms}
\title{
  Compute a Survival Curve from a Cox model
}
\description{
Computes the predicted survivor function for a Cox proportional 
hazards model. 
}
\usage{
\method{survfit}{coxph}(formula, newdata, 
        se.fit=TRUE, conf.int=.95, individual=FALSE, stype=2, ctype,
        conf.type=c("log","log-log","plain","none", "logit", "arcsin"),
        censor=TRUE, start.time, id, influence=FALSE,
        na.action=na.pass, type, ...)
\method{survfit}{coxphms}(formula, newdata, 
        se.fit=FALSE, conf.int=.95, individual=FALSE, stype=2, ctype,
        conf.type=c("log","log-log","plain","none", "logit", "arcsin"),
        censor=TRUE, start.time, id, influence=FALSE,
        na.action=na.pass, type, p0=NULL, ...)
}
\arguments{
  \item{formula}{
    A \code{coxph} object. 
  }
  \item{newdata}{
    a data frame with the same variable names as those that appear 
    in the \code{coxph} formula. One curve is produced per row.
    The curve(s) produced will be representative of a cohort whose 
    covariates correspond to the values in \code{newdata}. 
  }
  \item{se.fit}{
    a logical value indicating whether standard errors should be 
    computed.  Default is \code{TRUE} for standard models, \code{FALSE}
    for multi-state (code not yet present for that case.)
  }

  \item{conf.int}{
    the level for a two-sided confidence interval on the survival curve(s). 
    Default is 0.95. 
  }
  \item{individual}{deprecated argument, replaced by the general
    \code{id}}
  \item{stype}{computation of the survival curve, 1=direct, 2=
    exponenial of the cumulative hazard.}
  \item{ctype}{whether the cumulative hazard computation should have a
    correction for ties, 1=no, 2=yes.}
  \item{conf.type}{
    One of \code{"none"}, \code{"plain"}, \code{"log"} (the default),
    \code{"log-log"} or \code{"logit"}.  Only
    enough of the string to uniquely identify it is necessary.
    The first option causes confidence intervals not to be
    generated.  The second causes the standard intervals
    \code{curve +- k *se(curve)}, where k is determined from
    \code{conf.int}.  The log option calculates intervals based on the
    cumulative hazard or log(survival). The log-log option uses
    the log hazard or log(-log(survival)), and the logit
    log(survival/(1-survival)).
  }
  \item{censor}{if FALSE time points at which there are no events (only
    censoring) are not included in the result.}
  \item{id}{optional variable name of subject identifiers.  If this is
    present, it will be search for in the \code{newdata} data frame.
    Each group of rows in \code{newdata} with the same subject id represents
    the covariate path through time of a single subject, and the result
    will contain one curve per subject.  If the \code{coxph} fit had
    strata then that must also be specified in \code{newdata}.
    If \code{newid} is not present, then each
    individual row of \code{newdata} is presumed to represent a distinct
    subject.}

  \item{start.time}{optional starting time, a single numeric value.
    If present the returned curve contains survival after
    \code{start.time} conditional on surviving to \code{start.time}.
  }
  \item{influence}{option to return the influence values}
  \item{na.action}{the na.action to be used on the newdata argument}
  \item{type}{older argument that encompassed \code{stype} and
    \code{ctype}, now deprecated}
  \item{p0}{optional, a vector of probabilities.  The returned curve
    will be for a cohort with this mixture of starting states.  Most
    often a single state is chosen}
     \item{\dots}{for future methods}
}
\value{
an object of class \code{"survfit"}.  
See \code{survfit.object} for 
details. Methods defined for survfit objects are  
\code{print}, \code{plot}, 
\code{lines}, and \code{points}. 
}
\details{
  This routine produces Pr(state) curves based on a \code{coxph}
  model fit. For single state models it produces the single curve for
  S(t) = Pr(remain in initial state at time t), known as the survival
  curve; for multi-state models a matrix giving probabilities for all states.
  The \code{stype} argument states the type of estimate, and defaults
  to the exponential of the cumulative hazard, better known as the Breslow
  estimate.  For a multi-state Cox model this involves the exponential
    of a matrix. 
  The argument \code{stype=1} uses a non-exponential or `direct'
  estimate.  For a single endpoint coxph model the code evaluates the
  Kalbfleich-Prentice estimate, and for a multi-state model it uses an
  analog of the Aalen-Johansen estimator.  The latter approach is the
    default in the \code{mstate} package.

  The \code{ctype} option affects the estimated cumulative hazard, and
  if \code{stype=2} the estimated P(state) curves as well.  If not
  present it is chosen so as to be concordant with the 
  \code{ties} option in the \code{coxph} call. (For multistate
  \code{coxphms} objects only \code{ctype=1} is currently implemented.)
  Likewise
  the choice between a model based and robust variance estimate for the
  curve will mirror the choice made in the \code{coxph} call,
  any clustering is also inherited from the parent model.

  If the \code{newdata} argument is missing, then a curve is produced
  for a single "pseudo" subject with
  covariate values equal to the \code{means} component of the fit.
  The resulting curve(s) almost never make sense, but 
  the default remains due to an unwarranted attachment to the option shown by
  some users and by other packages.  Two particularly egregious examples
  are factor variables and interactions.  Suppose one were studying
  interspecies transmission of a virus, and the data set has a factor
  variable with levels ("pig", "chicken") and about equal numbers of
  observations for each.  The ``mean'' covariate level will be 0.5 --
  is this a flying pig?  As to interactions assume data with sex coded as 0/1,
  ages ranging from 50 to 80, and a model with age*sex.  The ``mean''
  value for the age:sex interaction term will be about 30, a value
  that does not occur in the data.
  Users are strongly advised to use the newdata argument.
  For these reasons predictions from a multistate coxph model require
  the newdata argument.

  If the \code{coxph} model contained an offset term, then the data set
  in the \code{newdata} argument should also contain that variable.

When the original model contains time-dependent covariates, then the
path of that covariate through time needs to be specified in order to
obtain a predicted curve. This requires \code{newdata} to contain
multiple lines for each hypothetical subject which gives the covariate
values, time interval, and strata for each line (a subject can change
strata), along with an \code{id} variable
 which demarks which rows belong to each subject.
The time interval must have the same (start, stop, status)
variables as the original model: although the status variable is not
used and thus can be set to a dummy value of 0 or 1, it is necessary for
the response to be recognized as a \code{Surv} object.
Last, although predictions with a time-dependent covariate path can be
useful, it is very easy to create a prediction that is senseless.  Users
are encouraged to seek out a text that discusses the issue in detail.

When a model contains strata but no time-dependent covariates the user
of this routine has a choice.
If newdata argument does not contain strata variables then the returned
object will be a matrix of survival curves with one row for each strata
in the model and one column for each row in newdata.
(This is the historical behavior of the routine.)
If newdata does contain strata variables, then the result will contain
one curve per row of newdata, based on the indicated stratum of the
original model.  In the rare case of a model with strata by covariate
interactions the strata variable must be included in newdata, the
routine does not allow it to be omitted (predictions become too confusing).
(Note that the model Surv(time, status) ~ age*strata(sex) expands internally to
strata(sex) + age:sex; the sex variable is needed for the second term
of the model.)

See \code{\link{survfit}} for more details about the counts (number of
events, number at risk, etc.)
}

\section{Notes}{
If the following pair of lines is used inside of another function then
the \code{model=TRUE} argument must be added to the coxph call:
\code{fit <- coxph(...); survfit(fit)}.
This is a consequence of the non-standard evaluation process used by the
\code{model.frame} function when a formula is involved.

Let \eqn{\log[S(t; z)]}{log(S(t; z))} be the log of the survival curve
for a fixed covariate vector \eqn{z}, then
\eqn{\log[S(t; x)]= e^{(x-z)\beta}\log[S(t; z)]}{log(S(t; x))\exp((x-z)*beta)}
is the log of the curve for any new covariate vector \eqn{x}.  
There is an unfortunate tendency to refer to the reference curve with
\eqn{z=0} as `THE' baseline hazard.  However, any \eqn{z} can be used as
the reference point, and more importantly, if \eqn{x-z} is large the
compuation can suffer severe roundoff error.  It is always safest to
provide the desired \eqn{x} values directly via \code{newdata}.
}

\section{References}{
Fleming, T. H. and Harrington, D. P. (1984).  Nonparametric estimation of the 
survival distribution in censored data.  \emph{Comm. in Statistics}  
\bold{13}, 2469-86. 

Kalbfleisch, J. D. and Prentice, R. L. (1980).
\emph{The Statistical Analysis of Failure Time Data.}
New York:Wiley. 

Link, C. L. (1984). Confidence intervals for the survival 
function using Cox's proportional hazards model with  
covariates.  \emph{Biometrics}  
\bold{40}, 601-610.

Therneau T and Grambsch P (2000), Modeling Survival Data: Extending the
Cox Model, Springer-Verlag.

Tsiatis, A. (1981). A large sample study of the estimate 
for the integrated hazard function in Cox's regression 
model for survival data. \emph{Annals of Statistics}  
\bold{9}, 93-108. 
}
\seealso{
\code{\link{print.survfit}},  
\code{\link{plot.survfit}},  
\code{\link{lines.survfit}},   
\code{\link{coxph}},  
\code{\link{Surv}},  
\code{\link{strata}}.   
}
\keyword{survival}
