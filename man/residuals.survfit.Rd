\name{residuals.survfit}
\alias{residuals.survfit}
\title{IJ residuals from a survfit object.}
\description{
  Return infinitesimal jackknife residuals from a survfit object, for
  the survival, cumulative hazard, or restricted mean time in state (RMTS).
}
\usage{
\method{residuals}{survfit}(object, times, 
    type="pstate", collapse=FALSE, weighted= collapse, data.frame=FALSE,
    extra = FALSE, ...)
}
\arguments{
  \item{object}{a \code{survfit} object}
  \item{times}{a vector of times at which the residuals are desired}
  \item{type}{the type of residual, see below}
  \item{collapse}{add the residuals for all subjects in a cluster}
  \item{weighted}{weight the residuals by each observation's weight}
  \item{data.frame}{if FALSE return a matrix or array}
  \item{extra}{return extra information when \code{data.frame=FALSE}.
    (This is used internally by the psuedo function.)
  \item{...}{arguments for other methods}
}
\details{
This function is designed to efficiently compute the per-observation
residuals for a Kaplan-Meier or Aalen-Johansen curve, also known as
infinitesimal jackknife (IJ) values, at a small number of time points.
Common usages are the creation of psuedo-values (via the \code{pseudo} function)
and IJ estimates of variance.
The residuals matrix has a value for each observation and time point
pair.
For a multi-state model the state will be a third dimension.

The residuals are the impact of each observation or cluster on the
resulting probability in state curves at the given time points,
the cumulative hazard curve at those time points,
or the expected sojourn time in each state up to the given time points.
For a simple Kaplan-Meier the \code{survfit} object contains only the
probability in the "initial" state, i.e., the survival fraction.
In this case the sojourn time, the expected amount of time spent in
the initial state, up to the specified endpoint, is commonly known as the
restricted mean survival time (RMST).
For a multistate model this same quantity is more often referred to as the
restricted mean time in state (RMTS).
It can be computed as the area under the respective probability in state curve.

The program allows any of \code{pstate}, \code{surv}, \code{cumhaz},
\code{chaz}, \code{sojourn}, \code{rmst}, \code{rmts} or \code{auc}
for the type argument, ignoring upper/lowercase, so
users can choose whichever abbreviation they like best.

When \code{collapse=TRUE} the result has the cluster identifier (which
defaults to the \code{id} variable) as the dimname for the first
dimension. 
If the \code{fit} object contains more than one curve, and the same
identifier is reused in two different curves this approach does not work
and the routine will stop with an error.
In principle this is not necessary, e.g., the result could contain two rows
with the same label, showing the separate effect on each curve,
but this was deemed too confusing. 
}

\value{A matrix or array with one row per observation or cluster, and one column
  for each value in \code{times}.  For a multi-state model the three
  dimensions are observation, state, and time.  For cumulative hazard,
  the second dimension is the set of transitions.  (A competing risks
  model for instance has 3 states and 2 transitions.)
}
\note{The first column of the data frame identifies the origin of the
  row.  If there was an \code{id} variable in the \code{survfit} call it
  will contain the values of that variable and be labeled with the
  variable name, or "(id)" if there was an expression rather than a
  name. (For example, \code{survfit(....  id= abc$def[z])}).  If there
  was no \code{id} variable the label will be "(row)", and the column
  will contain the row number of the survfit data.  For a matrix result
  the first component of dimnames has similar structure.
}
\seealso{\code{\link{survfit}}, \code{\link{survfit.formula}}
}
\examples{
fit <- survfit(Surv(time, status) ~ x, aml)
resid(fit, times=c(24, 48), type="RMTS")
}
% \keyword{ survival }

