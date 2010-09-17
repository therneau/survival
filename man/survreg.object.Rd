\name{survreg.object}
\alias{survreg.object}
\alias{print.survreg}
\alias{summary.survreg}

\title{
Parametric Survival Model Object
}
\description{
This class of objects is returned by the \code{survreg} function
to represent a fitted parametric survival model.
Objects of this class have methods for the functions \code{print},
\code{summary}, \code{predict}, and \code{residuals}.
}
\section{COMPONENTS}{
The following components must be included in a legitimate \code{survreg} object.
\describe{
\item{coefficients}{
the coefficients of the \code{linear.predictors}, which multiply  the
columns of the model
matrix.
It does not include the estimate of error (sigma).
The names of the coefficients are the names of the
single-degree-of-freedom effects (the columns of the
model matrix).
If the model is over-determined there will
be missing values in the coefficients corresponding to non-estimable
coefficients.
}
\item{icoef}{
coefficients of the baseline model, which will contain the intercept
and log(scale), or mulitple scale factors for a stratified model.
}
\item{var}{
the variance-covariance matrix for the parameters, including the log(scale)
parameter(s).
}
\item{loglik}{
a vector of length 2, containing the log-likelihood for the baseline and
full models. 
}
\item{iter}{
the number of iterations required
}
\item{linear.predictors}{
the linear predictor for each subject.
}
\item{df}{
the degrees of freedom for the final model.  For a penalized model
this will be a vector with one element per term.
}
\item{scale}{
the scale factor(s), with length equal to the number of strata.
}
\item{idf}{
degrees of freedom for the initial model.
}
\item{means}{
a vector of the column means of the coefficient matrix.
}
\item{dist}{
the distribution used in the fit.


The object will also have the following components found in 
other model results (some are optional):
\code{linear predictors}, \code{weights}, \code{x}, \code{y}, \code{model}, 
\code{call}, \code{terms} and \code{formula}.
See \code{lm}.
}}}
\seealso{
\code{\link{survreg}}, \code{\link{lm}}
}
\keyword{regression}
\keyword{survival}
% Converted by Sd2Rd version 0.3-2.
