\name{survdiff}
\alias{survdiff}
\alias{print.survdiff}
\title{
Test Survival Curve Differences
}
\description{
Tests if there is a difference between two or more survival curves using
the \eqn{G^\rho}{G-rho} family of tests, or for a single curve against a known alternative.
}
\usage{
survdiff(formula, data, subset, na.action, rho=0, timefix=TRUE)
}
\arguments{
\item{formula}{
a formula expression as for other survival models, of the form
\code{Surv(time, status) ~ predictors}.  For a one-sample test, the predictors
must consist of a single \code{offset(sp)} term, where \code{sp} is a vector giving the
survival probability of each subject.  For a k-sample test, each unique
combination of predictors defines a subgroup.
A \code{strata} term may be used to produce a stratified test.
To cause missing values in the predictors to be treated as a separate
group, rather than being omitted, use the \code{strata} function with its
\code{na.group=T} argument.
}
\item{data}{
an optional data frame in which to interpret the variables occurring in the
formula.
}
\item{subset}{
expression indicating which subset of the rows of data should be used in
the fit.  This can be a logical vector (which is replicated to have
length equal to the number of observations), a numeric vector indicating
which observation numbers are to be included (or excluded if negative),
or a character vector of row names to be included.  All observations are
included by default.
}
\item{na.action}{
a missing-data filter function.  This is applied to the \code{model.frame} after any
subset argument has been used.  Default is \code{options()$na.action}.
}
\item{rho}{
a scalar parameter that controls the type of test.
}
\item{timefix}{process times through the \code{aeqSurv} function to
  eliminate potential roundoff issues.}

}
\value{
a list with components:

\item{n}{
the number of subjects in each group.
}
\item{obs}{
the weighted observed number of events in each group.
If there are strata, this will be a matrix with one column per stratum.
}
\item{exp}{
the weighted expected number of events in each group.
If there are strata, this will be a matrix with one column per stratum.
}
\item{chisq}{
the chisquare statistic for a test of equality.
}
\item{var}{
the variance matrix of the test.
}
\item{strata}{
optionally, the number of subjects contained in each stratum.
}
\item{pvalue}{the p-value corresponding to the Chisquare statistic}
}
\section{Description}{
This function implements the G-rho family of
Harrington and Fleming (1982), with weights on each death of \eqn{S(t)^\rho}{S(t)^rho},
where \eqn{S(t)}{S} is the Kaplan-Meier estimate of survival.
With \code{rho = 0} this is the log-rank or Mantel-Haenszel test,
and with \code{rho = 1} it is equivalent to the Peto & Peto modification
of the Gehan-Wilcoxon test.

Peto and Peto show that the Gehan-Wilcoxon test can be badly biased if
the two groups have different censoring patterns, and proposed an
alternative.  Prentice and Marek later showed an actual example where this
issue occurs.  For most data sets the Gehan-Wilcoxon and
Peto-Peto-Prentice variant will hardly differ, however.

If the right hand side of the formula consists only of an offset term,
then a one sample test is done.
To cause missing values in the predictors to be treated as a separate
group, rather than being omitted, use the \code{factor} function with its
\code{exclude} argument to recode the righ-hand-side covariate.
}
\references{
  Harrington, D. P. and Fleming, T. R. (1982).
A class of rank test procedures for censored survival data.
 Biometrika, 553-566.

  Peto R. Peto and Peto, J. (1972) Asymptotically efficient rank invariant
  test procedures (with discussion), JRSSA, 185-206.

  Prentice, R. and Marek, P. (1979)  A qualitative discrepancy between
  censored data rank tests, Biometics, 861--867.
}
\examples{
## Two-sample test
survdiff(Surv(futime, fustat) ~ rx,data=ovarian)

## Stratified 7-sample test

survdiff(Surv(time, status) ~ pat.karno + strata(inst), data=lung)

## Expected survival for heart transplant patients based on
## US mortality tables
expect <- survexp(futime ~ 1, data=jasa, cohort=FALSE,
                  rmap= list(age=(accept.dt - birth.dt), sex=1, year=accept.dt),
                  ratetable=survexp.us)
## actual survival is much worse (no surprise)
survdiff(Surv(jasa$futime, jasa$fustat) ~ offset(expect))

# The free light chain data set is close to the population.
e2 <- survexp(futime ~ 1, data=flchain, cohort=FALSE,
              rmap= list(age= age*365.25, sex=sex, 
                         year=as.Date(paste0(sample.yr, "-07-01"))),
              ratetable= survexp.mn)
survdiff(Surv(futime, death) ~ offset(e2), flchain)
}
\keyword{survival}

