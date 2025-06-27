\name{finegray}
\alias{finegray}
\title{Create data for a Fine-Gray model}
\description{
  The Fine-Gray model can be fit by first creating a special data set,
  and then fitting a weighted Cox model to the result.  This routine creates the
  data set.
}
\usage{
finegray(formula, data, weights, subset, na.action= na.pass, etype,
    prefix="fg", count, id, timefix=TRUE)
}
\arguments{
  \item{formula}{a standard model formula, with survival on the left and
    covariates on the right.  
}
\item{data}{an optional data frame, list or environment (or object
  coercible by as.data.frame to a data frame) containing the variables
  in the model.
}
  \item{weights}{optional vector of observation weights}
  \item{subset}{
    an optional vector specifying a subset of observations to be used in
    the fitting process.
}
  \item{na.action}{
    a function which indicates what should happen when the data contain
    NAs.  The default is set by the na.action setting of options.
}
\item{etype}{
   the event type for which a data set will be generated.  The default is
   to use whichever is listed first in the multi-state survival object.
 }

 \item{prefix}{the routine will add 4 variables to the data set: a start
   and end time for each interval, status, and a weight for the
   interval. The default names of these are "fgstart", "fgstop", "fgstatus",
   and "fgwt"; the \code{prefix} argument determines the initial portion of
   the new names.
 }

 \item{count}{a variable name in the output data set for an optional
   variable that will contain
  the the replication count for each row of the input data.  If a row is
  expanded into multiple lines it will contain 1, 2, etc.
}

\item{id}{optional, the variable name in the data set which identifies
  subjects.}
\item{timefix}{process times through the \code{aeqSurv} function to
  eliminate potential roundoff issues.}
}

\details{
  The function expects a multi-state survival expression or variable as
  the left hand side of the formula, e.g. \code{Surv(atime, astat)}
  where \code{astat} is a factor whose first level represents censoring
  and remaining levels are states.  The output data set will contain simple
  survival data (status = 0 or 1) for a single endpoint of interest.
  For exposition
  call this endpoint A and lump all others as endpoint B.
  In the output data set subjects who experience endpoint B become
  censored observations 
  whose times are artificially extended to the right, with a
  decreasing case weight from interval to interval.  
  The output data set will normally contain many more rows than the
  input.

  The algorithm allows for delayed entry, and only a limited form of
  time-dependent covariates.  That is, when subjects with endpoint B are
  extended, those future covariate values stay constant; so there is an
  implicit assumption that no more changes would have occurred if
  the event had not intervened and follow-up had been longer. 
  For predictable time-dependent covariates the final data set could be
  further processed to fix this, but this is not included in the
  function.  Geskus for example considers an example with different
  calendar epochs, corresponding to a change in standard medical
  practice for the disese, as a covariate.
  dependent covariates.  
  If there are time dependent covariates or delayed entry, e.g.., the input data
  set had \code{Surv(entry, exit, stat)} as the left hand side, then
  an \code{id} statement is required.  The program does data checks
  in this case, and needs to know which rows belong to each subject.

  The output data set will often have gaps. Say that there were events
  at time 50 and 100 (and none between) and censoring at 60, 70, and 80.
  Formally, a non event subjects at risk from 50 to 100 will have
  different weights in each of
  the 3 intervals 50-60, 60-70, and 80-100, but because the middle
  interval does not span any event times the subsequent Cox model will
  never use that row.  The \code{finegray} output omits such rows.
  
  In the simulation results of Geskus, he found that the robust sandwich
  estimate of the standard error was equal to the Fine and Gray estimate
  up to 3 digits, we thus recomment including the \code{cluster} option
  in the subsequent \code{coxph} call.
}
\value{a data frame}
\references{
  Fine JP and Gray RJ (1999) A proportional hazards model for the
  subdistribution of a competing risk. JASA 94:496-509.

  Geskus RB (2011). Cause-Specific Cumulative Incidence Estimation and the
  Fine and Gray Model Under Both Left Truncation and Right Censoring.
  Biometrics 67, 39-49.

}
\author{Terry Therneau}
\seealso{\code{\link{coxph}}, \code{\link{aeqSurv}}}

\examples{
# Treat time to death and plasma cell malignancy as competing risks
etime <- with(mgus2, ifelse(pstat==0, futime, ptime))
event <- with(mgus2, ifelse(pstat==0, 2*death, 1))
event <- factor(event, 0:2, labels=c("censor", "pcm", "death"))

# FG model for PCM
pdata <- finegray(Surv(etime, event) ~ ., data=mgus2)
fgfit <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex,
                     weight=fgwt, data=pdata, cluster=id)

# Compute the weights separately by sex
adata <- finegray(Surv(etime, event) ~ . + strata(sex),
             data=mgus2, na.action=na.pass)
}
\keyword{survival}

