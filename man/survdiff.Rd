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
survdiff(formula, data, subset, na.action, rho=0)
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
}}
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
}}
\section{METHOD}{
This function implements the G-rho family of
Harrington and Fleming (1982), with weights on each death of \eqn{S(t)^\rho}{S(t)^rho},
where \eqn{S(t)}{S} is the Kaplan-Meier estimate of survival.
With \code{rho = 0} this is the log-rank or Mantel-Haenszel test,
and with \code{rho = 1} it is equivalent to the Peto & Peto modification
of the Gehan-Wilcoxon test.


If the right hand side of the formula consists only of an offset term,
then a one sample test is done.
To cause missing values in the predictors to be treated as a separate
group, rather than being omitted, use the \code{factor} function with its
\code{exclude} argument.
}
\references{
Harrington, D. P. and Fleming, T. R. (1982).
A class of rank test procedures for censored survival data.
\emph{Biometrika}
\bold{69}, 553-566.}
\examples{
## Two-sample test
survdiff(Surv(futime, fustat) ~ rx,data=ovarian)

## Stratified 7-sample test

survdiff(Surv(time, status) ~ pat.karno + strata(inst), data=lung)

## Expected survival for heart transplant patients based on
## US mortality tables
expect <- survexp(futime ~ ratetable(age=(accept.dt - birth.dt),
     sex=1,year=accept.dt,race="white"), jasa, cohort=FALSE,
     ratetable=survexp.usr)
## actual survival is much worse (no surprise)
survdiff(Surv(jasa$futime, jasa$fustat) ~ offset(expect))
}
\keyword{survival}
% Converted by Sd2Rd version 0.3-2. and hand edited.
