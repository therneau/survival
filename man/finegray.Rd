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
  In the output data set subjects who did not experience the event of
  interest become censored subjects
  whose times are artificially extended over multiple intervals, with a
  decreasing case weight from interval to interval.  
  The output data set will normally contain many more rows than the
  input.

  Time dependent covariates are allowed, but not (currently) delayed
  entry.  If there are time dependent covariates, e.g.., the input data
  set had \code{Surv(entry, exit, stat)} as the left hand side, then
  an \code{id} statement is required.  The program does data checks
  in this case, and needs to know which rows belong to each subject.
  
  See the competing risks vignette for more details.
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
                     weight=fgwt, data=pdata)

# Compute the weights separately by sex
adata <- finegray(Surv(etime, event) ~ . + strata(sex),
             data=mgus2, na.action=na.pass)
}
\keyword{survival}

