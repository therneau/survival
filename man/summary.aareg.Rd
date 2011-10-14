\name{summary.aareg}
\alias{summary.aareg}
\title{
  Summarize an aareg fit
}
\description{
Creates the overall test statistics for an Aalen additive regression model
}
\usage{
\method{summary}{aareg}(object, maxtime, test=c("aalen", "nrisk"), scale=1,...)
}
\arguments{
\item{object}{
the result of a call to the \code{aareg} function
}
\item{maxtime}{
truncate the input to the model at time "maxtime"
}
\item{test}{
the relative time weights that will be used to compute the test
}
\item{scale}{
scales the coefficients.  
For some data sets, the coefficients of the Aalen model will be very
small (10-4); this simply multiplies the printed values by a constant,
say 1e6, to make the printout easier to read.
}
\item{\dots}{for future methods}
}
\value{
a list is returned with the following components

\item{ table }{
a matrix with rows for the intercept and each covariate, and columns
giving a slope estimate, the test statistic, it's standard error,
the z-score and a p-value
}
\item{ test }{
the time weighting used for computing the test statistics
}
\item{ test.statistic }{
the vector of test statistics
}
\item{ test.var }{
the model based variance matrix for the test statistic
}
\item{ test.var2 }{
optionally, a robust variance matrix for the test statistic
}
\item{ chisq }{
the overall test (ignoring the intercept term) for significance of
any variable
}
\item{ n }{
a vector containing the number of observations, the number of unique death
times used in the computation, and the total number of unique death times
}
}
\details{
It is not uncommon for the very right-hand tail of the plot to have
large outlying values, particularly for the standard error.  
The \code{maxtime} parameter can then be used to
truncate the range so as to avoid these.
This gives an updated value for the test statistics, without refitting the
model.

The slope is based on a weighted linear regression to the cumulative
coefficient plot, and may be a useful measure of the overall size
of the effect.  For instance when two models include a common variable,
"age" for instance, this may help to assess how much the fit changed due
to the other variables, in leiu of overlaying the two plots.  (Of course
the plots are often highly non-linear, so it is only a rough substitute).
The slope is not directly related to the test statistic, as the latter
is invariant to any monotone transformation of time.
}
\seealso{
aareg, plot.aareg
}
\examples{
afit <- aareg(Surv(time, status) ~ age + sex + ph.ecog, data=lung,
     dfbeta=TRUE)
summary(afit)
\dontrun{
              slope   test se(test) robust se     z        p 
Intercept  5.05e-03    1.9     1.54      1.55  1.23 0.219000
      age  4.01e-05  108.0   109.00    106.00  1.02 0.307000
      sex -3.16e-03  -19.5     5.90      5.95 -3.28 0.001030
  ph.ecog  3.01e-03   33.2     9.18      9.17  3.62 0.000299

Chisq=22.84 on 3 df, p=4.4e-05; test weights=aalen
}

summary(afit, maxtime=600)
\dontrun{
              slope   test se(test) robust se      z        p 
Intercept  4.16e-03   2.13     1.48      1.47  1.450 0.146000
      age  2.82e-05  85.80   106.00    100.00  0.857 0.392000
      sex -2.54e-03 -20.60     5.61      5.63 -3.660 0.000256
  ph.ecog  2.47e-03  31.60     8.91      8.67  3.640 0.000271

Chisq=27.08 on 3 df, p=5.7e-06; test weights=aalen
}}
\keyword{survival}

