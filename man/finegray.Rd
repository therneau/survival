\name{finegray}
\alias{finegray}
\title{Create data for a Fine-Gray model}
\description{
  The Fine-Gray model can be fit by first creating a special data set,
  and then fitting a Cox model to the result.  This routine creates the
  data.
}
\usage{
finegray(formula, data, subset, na.action= na.pass, endpoint,
    count="", tstart="fgstart", end="fgend", status="fgstat", wt="fgwt")
}
\arguments{
  \item{formula}{a standard model formula, with survival on the left and
    covariates on the right.  
}
\item{data}{an optional data frame, list or environment (or object
  coercible by as.data.frame to a data frame) containing the variables
  in the model.
}
  \item{subset}{
    an optional vector specifying a subset of observations to be used in
    the fitting process.
}
  \item{na.action}{
    a function which indicates what should happen when the data contain
    NAs.  The default is set by the na.action setting of options.
}
\item{endpoint}{
   the endpoint for which a data set will be generated.  The default is
   to use the first endpoint state.
}
\item{count}{a variable name in the output data set that will contain
  the the replication count for each row of the input data.  If a row is
  expanded into multiple lines it will contain 1, 2, etc.
}

\item{tstart, end, status}{output variable names for the survival time
  variables.  See the details below.}
\item{wt}{output variable name for the created weight variable.}
}
\details{
  The function expects a multi-state survival expression or variable as
  the left hand side of the formula, e.g. \code{Surv(atime, astat)}
  where \code{astat} is a factor whose first level represents censoring
  and remaining levels are states.
  In the output data set subjects who did not experience the event of
  interest, i.e., the \code{endpoint} argument, become censored subjects
  whose times are artificially extended over multiple intervals, with a
  decreasing case weight from interval to interval.  
  The output data set will normally contain many more rows than the
  input.

  The output data set contains a simple (start, end, event) survival
  times where the event is 0 or 1.  
  To the extent possible, the created output data set will reuse variable
  names from the input, i.e., in the example just given the output data set
  would contain the variables "tstart", "atime", and "status".
  However, the parsing portion is (purposefully) unsophisticated, and if
  confused the function will fall back to the default names given in the
  \code{tstart}, \code{end} and \code{status} arguments.  The expression
  \code{Surv( (atime), astat)} is sufficient, for example.

  Time dependent covariates are allowed, but not (currently) delayed
  entry.
  See the competing risks vignette for more details.
}
\value{a data frame}
\references{
  Fine JP and Gray RJ (1999) A proportional hazards model for the
  subdistribution of a competing risk. JASA 94:496-509.

  Geskus RB (2011). Cause-Specific Cumulative Incidence Estimation and the
  Fine and Gray Model Under Both Left Truncation and Right Censoring.
  Biometrics 67, 39–49.

}
\author{Terry Therneau}
\seealso{\code{\link{coxph}}}
}
\examples{
# Treat time to death and plasma cell malignancy as competing risks
etime <- with(mgus2, ifelse(pstat==0, futime, ptime))
event <- with(mgus2, ifelse(pstat==0, 2*death, 1))
event <- factor(event, 0:2, labels=c("censor", "pcm", "death"))

# multi-state curves
curves <- survfit(Surv(etime, event) ~ sex, data=mgus2)

# FG model for PCM
pdata <- finegray(Surv(etime, event) ~ ., data=mgus2)
fgfit <- coxph(Surv(etime, event) ~ age + sex, weight=fgwt, data=pdata)
}
\keyword{survival}

