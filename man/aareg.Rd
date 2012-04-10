\name{aareg}
\alias{aareg}
\title{
  Aalen's additive regression model for censored data
}
\description{
Returns an object of class \code{"aareg"} that
represents an Aalen model.
}
\usage{
aareg(formula, data, weights, subset, na.action,
   qrtol=1e-07, nmin, dfbeta=FALSE, taper=1,
   test = c('aalen', 'variance', 'nrisk'),
    model=FALSE, x=FALSE, y=FALSE) 
}
\arguments{
\item{formula}{
a formula object, with the response on the left of a `~' operator and 
the terms,
 separated by \code{+} operators, on the right. 
The response must be a \code{Surv} object.
Due to a particular computational approach that is used, the model
MUST include an intercept term.  If "-1" is used in the model
formula the program will ignore it.
}
\item{data}{
data frame in which to interpret the variables named in the
\code{formula},
\code{subset}, and \code{weights} arguments.
This may also be a single number to handle
some speci al cases -- see below for details. If
\code{data} is missi ng, the variables in the
model formula should be in the search path.
}
\item{weights}{
vector of observation weights. If supplied, the fitting algorithm
minimizes the sum of the weights multiplied by the squared residuals
(see below for additional technical details). The length of
\code{weights} must be the same as the number of
observations. The weights must be nonnegative and it i s recommended
that they be strictly positive, since zero weights are ambiguous.  To
exclude particular observations from the model, use the
\code{subset} argument instead of zero weights.
}
\item{subset}{
expression specifying which subset of observations should be used in
the fit. Th is can be a logical vector (which is replicated to have
length equal to the numb er of observations), a numeric vector
indicating the observation numbers to be i ncluded, or a character
vector of the observation names that should be included.  All
observations are included by default.
}
\item{na.action}{
a function to filter missing data. This is applied to the
\code{model.fr ame} after any
\code{subset} argument has be en applied. The
default is \code{na.fail}, which returns a n
error if any missing values are found. An alternative is
\code{na.excl ude}, which deletes observations
that contain one or more missing values.
}
\item{qrtol}{
tolerance for detection of singularity in the QR decomposition
}
\item{nmin}{
minimum number of observations for an estimate; defaults to 3 times the
number of covariates.
This essentially truncates the computations near the tail of the data
set, when n is small and the calcualtions can become numerically
unstable.
}
\item{dfbeta}{
should the array of dfbeta residuals be computed.  This implies computation
of the sandwich variance estimate.
The residuals will always be computed if there is a
\code{cluster} term in the model formula.
}
\item{taper}{
allows for a smoothed variance estimate.
Var(x), where x is the set of covariates, is an important component of the
calculations for the Aalen regression model.  
At any given time point t, it is computed over all subjects who are still
at risk at time t.
The tape argument allows smoothing these estimates,
for example \code{taper=(1:4)/4} would cause
the variance estimate used at any event time to be a weighted average
of the estimated variance matrices at the last 4 death times,
with a weight of 1 for the current death time and decreasing to 1/4 for
prior event times.
The default value gives the standard Aalen model.
}
\item{test}{
selects the weighting to be used, for computing an overall ``average''
coefficient vector over time and the subsequent test for equality to zero.
}
\item{model, x, y }{
should copies of the model frame, the x matrix of predictors, or the
response vector y be included in the saved result.
}
}
\value{
  an object of class \code{"aareg"} 
  representing the fit, with the following components:
  \item{n}{vector containing the number of observations in the data set,
    the number of event times, and the number of event times used in the
    computation}
  \item{times}{vector of sorted event times, which may contain
    duplicates}
  \item{nrisk}{vector containing the number of subjects at risk, of the
    same length as \code{times}}
  \item{coefficient}{matrix of coefficients, with one row per event and
    one column per covariate}
  \item{test.statistic}{the value of the test statistic, a vector with
    one element per covariate}
  \item{test.var}{variance-covariance matrix for the test}
  \item{test}{the type of test; a copy of the \code{test} argument
    above}
  \item{tweight}{matrix of weights used in the computation, one row per
    event}
  \item{call}{a copy of the call that produced this result}
}
\details{
The Aalen model assumes that the cumulative hazard H(t) for a subject
can be expressed as a(t) + X B(t), where a(t) is a time-dependent intercept
term, X is the vector of covariates for the subject (possibly time-dependent),
and B(t) is a time-dependent matrix of coefficients.
The estimates are inheritly non-parametric; a fit of the model will normally
be followed by one or more plots of the estimates.

The estimates may become unstable near the tail of a data set, since the
increment to B at time t is based on the subjects still at risk at time
t.  The tolerance and/or nmin parameters may act to truncate the estimate
before the last death.
The \code{taper} argument can also be used to smooth
out the tail of the curve.
In practice, the addition of a taper such as 1:10 appears to have little
effect on death times when n is still reasonably large, but can considerably
dampen wild occilations in the tail of the plot.  
}
\section{References}{
Aalen, O.O. (1989). A linear regression model for the analysis of life times.
Statistics in Medicine, 8:907-925.

Aalen, O.O (1993). Further results on the non-parametric linear model in
survival analysis.  Statistics in Medicine. 12:1569-1588.
}
\seealso{
print.aareg, summary.aareg, plot.aareg
}
\examples{
# Fit a model to the lung cancer data set
lfit <- aareg(Surv(time, status) ~ age + sex + ph.ecog, data=lung,
                     nmin=1)
\dontrun{
lfit
Call:
aareg(formula = Surv(time, status) ~ age + sex + ph.ecog, data = lung, nmin = 1
        )

  n=227 (1 observations deleted due to missing values)
    138 out of 138 unique event times used

              slope      coef se(coef)     z        p 
Intercept  5.26e-03  5.99e-03 4.74e-03  1.26 0.207000
      age  4.26e-05  7.02e-05 7.23e-05  0.97 0.332000
      sex -3.29e-03 -4.02e-03 1.22e-03 -3.30 0.000976
  ph.ecog  3.14e-03  3.80e-03 1.03e-03  3.70 0.000214

Chisq=26.73 on 3 df, p=6.7e-06; test weights=aalen

plot(lfit[4], ylim=c(-4,4))  # Draw a plot of the function for ph.ecog
}
lfit2 <- aareg(Surv(time, status) ~ age + sex + ph.ecog, data=lung,
                  nmin=1, taper=1:10)
\dontrun{lines(lfit2[4], col=2)  # Nearly the same, until the last point}

# A fit to the mulitple-infection data set of children with
# Chronic Granuomatous Disease.  See section 8.5 of Therneau and Grambsch.
fita2 <- aareg(Surv(tstart, tstop, status) ~ treat + age + inherit +
                         steroids + cluster(id), data=cgd)
\dontrun{
  n= 203 
    69 out of 70 unique event times used

                     slope      coef se(coef) robust se     z        p
Intercept         0.004670  0.017800 0.002780  0.003910  4.55 5.30e-06
treatrIFN-g      -0.002520 -0.010100 0.002290  0.003020 -3.36 7.87e-04
age              -0.000101 -0.000317 0.000115  0.000117 -2.70 6.84e-03
inheritautosomal  0.001330  0.003830 0.002800  0.002420  1.58 1.14e-01
steroids          0.004620  0.013200 0.010600  0.009700  1.36 1.73e-01

Chisq=16.74 on 4 df, p=0.0022; test weights=aalen
}
}
\keyword{survival}
% docclass is function
