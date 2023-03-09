\name{pseudo}
\alias{pseudo}
\title{
Pseudo values for survival.
}
\description{
Produce pseudo values from a survival curve.
}
\usage{
pseudo(fit, times, type, addNA=TRUE, data.frame=FALSE, minus1=FALSE, ...)
}
\arguments{
  \item{fit}{a \code{survfit} object, or one that inherits that class.   
}
  \item{times}{
    a vector of time points, at which to evaluate the pseudo values.
    }
  \item{type}{
    the type of value, either the probabilty in state \code{pstate},
    the cumulative hazard \code{cumhaz} or the expected sojourn time in
    the state \code{sojourn}.  
  }
  \item{addNA}{If any observations were removed due to missing values
    in the \code{fit} object, add those rows (as NA) into the return.
    This causes the result of pseudo to match the original dataframe.
  }
  \item{data.frame}{if TRUE, return the data in "long" form as a
    data.frame with id, time, and pseudo as variables.}
  \item{minus1}{use n-1 as the multiplier rather than n}.
  \item{\dots}{
    other arguments to the \code{residuals.survfit} function, which does
 the majority of the work, e.g., \code{collapse} and \code{weighted}.
}
}
\details{
This function computes pseudo values based on a first order Taylor
series, also known as the "infinitesimal jackknife" (IJ) or "dfbeta"
residuals.  To be completely correct these results could perhaps be
called `IJ pseudo values' or even pseudo psuedo-values.
For moderate to large data, however, the resulting values will
be almost identical, numerically, to the ordinary jackknife.

A primary advantage of this approach is computational speed.
Other features, neither good nor bad, are that they will agree with
robust standard errors of other survival package estimates,
which are based on the IJ, and that the mean of the estimates, over
subjects, is exactly the underlying survival estimate.

For the \code{type} variable, \code{surv} is an acceptable synonym for
\code{pstate}, and \code{rmst} and \code{rmts} are equivalent to \code{sojourn}.
All of these are case insensitive.

The result from this routine is simply n times the IJ value, where n is
the number of subjects.
(If the the \code{survfit} call included and \code{id} option, n is
the number of unique id values, otherwise the number of rows in the data set.) 
IJ values are well defined for all variants of the Aalen-Johansen
estimate, as computed by the \code{survfit} function; indeed, they are
the basis for standard errors of the result.
Understanding of the properties of the pseudo-values, however, is still
evolving.  Validity has been shown for the simplest case (Kaplan-Meier),
for competing risks, and for the corresponding sojourn times.
On the other hand, one must be careful when the data includes
left-truncation (P. K. Andersen, personal communication), and also with
pseudo-values for the cumulative hazard.
As understanding evolves, treat this routine's results as a reseach
tool, not production, for the more complex models.
}
\value{
  A vector, matrix, or array.  The first dimension is always the number of
  observations in \code{fit} object, in the same order as the original
  data set (less any missing values that were removed when creating the
  survfit object);
  the second, if applicable, corresponds to \code{fit$states}, e.g.,
  multi-state
  survival, and the last dimension to the selected time points.
  (If there are multiple rows for a given id, there is only one
  pseudovalue per unique id.)

  For the data.frame option, a data frame containing values for id,
  time, and pseudo.  If the original \code{survfit} call contained an
  \code{id} statement, then the values in the \code{id} column will be
  taken from that variable.  If the \code{id} statement has a simple
  form, e.g., \code{id = patno}, then the name of the id column will
  be `patno', otherwise it will be named `(id)'.
}
\note{
  The code will be slightly faster if the \code{model=TRUE} option is
  used in the \code{survfit} call.  It may be essential if the
  survfit/pseudo pair is used inside another function.
}
\references{
PK Andersen and M Pohar-Perme, Pseudo-observations in surivival
analysis, Stat Methods Medical Res, 2010; 19:71-99
}
\seealso{
\code{\link{residuals.survfit}}
}
\examples{
fit1 <- survfit(Surv(time, status) ~ 1, data=lung)
yhat <- pseudo(fit1, times=c(365, 730))
dim(yhat)
lfit <- lm(yhat[,1] ~ ph.ecog + age + sex, data=lung)

# Restricted Mean Time in State (RMST) 
rms <- pseudo(fit1, times= 730, type='RMST') # 2 years
rfit <- lm(rms ~ ph.ecog + sex, data=lung)
rhat <- predict(rfit, newdata=expand.grid(ph.ecog=0:3, sex=1:2), se.fit=TRUE)
# print it out nicely
temp1 <- cbind(matrix(rhat$fit, 4,2))
temp2 <- cbind(matrix(rhat$se.fit, 4, 2))
temp3 <- cbind(temp1[,1], temp2[,1], temp1[,2], temp2[,2])
dimnames(temp3) <- list(paste("ph.ecog", 0:3), 
                        c("Male RMST", "(se)", "Female RMST", "(se)"))

round(temp3, 1)
# compare this to the fully non-parametric estimate
fit2 <- survfit(Surv(time, status) ~ ph.ecog, data=lung)
print(fit2, rmean=730)
# the estimate for ph.ecog=3 is very unstable (n=1), pseudovalues smooth it.
#
# In all the above we should be using the robust variance, e.g., svyglm, but
#  a recommended package can't depend on external libraries.
# See the vignette for a more complete exposition.
}
\keyword{ survival }

