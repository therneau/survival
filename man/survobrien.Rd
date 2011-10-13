\name{survobrien}
\alias{survobrien}
\title{
 O'Brien's Test for Association of a Single Variable with Survival
}
\description{
Peter O'Brien's test for association of a single variable with survival 
This test is proposed in Biometrics, June 1978. 
}
\usage{
survobrien(formula, data, subset, na.action, transform)
}
\arguments{
\item{formula}{
a valid formula for a cox model.
}
  \item{data}{
    a data.frame in which to interpret the variables named in 
    the \code{formula}, or in the \code{subset} and the \code{weights}
    argument. 
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

\item{transform}{the transformation function to be applied at each
  time point. The default is O'Brien's suggestion logit(tr) where
  tr = (rank(x)- 1/2)/ length(x) is the rank shifted to the range
  0-1 and logit(x) = log(x/(1-x)) is the logit transform.
}}
\value{
  a new data frame.  The response variables will be column names
  returned by the \code{Surv} function, i.e., "time" and "status" for
  simple survival data, or "start", "stop", "status" for counting
  process data.  Each individual event time is identified by the
  value of the variable \code{.strata.}.  Other variables retain
  their original names.  If a
predictor variable is a factor or is protected with \code{I()}, it is 
retained as is.  Other predictor variables have been replaced with 
time-dependent logit scores. 

The new data frame will have many 
more rows that the original data, approximately the original number
of rows * number of deaths/2. 
}
\section{Method}{
A time-dependent cox model can now be fit to the new data. 
The univariate statistic, as originally proposed, is equivalent to 
single variable score tests from the time-dependent model. 
This equivalence is the rationale for using the time dependent model as a 
multivariate extension of the original paper. 

In O'Brien's method, the x variables are re-ranked at each death time.  A 
simpler method, proposed by Prentice, ranks the data only once at the 
start. The results are usually similar. 
}
\references{
O'Brien, Peter, "A Nonparametric Test for Association with Censored Data", 
\emph{Biometrics} 34: 243-250, 1978. 
}

\note{
  A prior version of the routine returned new time variables rather than
  a strata.  Unfortunately, that strategy does not work if the original
  formula has a strata statement.  This new data set will be the same
  size, but the \code{coxph} routine will process it slightly faster.
  }
\seealso{
\code{\link{survdiff}} 
}
\keyword{survival}
\examples{
xx <- survobrien(Surv(futime, fustat) ~ age + factor(rx) + I(ecog.ps), 
			       data=ovarian) 
coxph(Surv(time, status) ~ age + strata(.strata.), data=xx) 
}
