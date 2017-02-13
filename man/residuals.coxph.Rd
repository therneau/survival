\name{residuals.coxph}
\alias{residuals.coxph.penal}
\alias{residuals.coxph.null}
\alias{residuals.coxph}
\title{
Calculate Residuals for a `coxph' Fit
}
\description{
Calculates martingale, deviance, score or Schoenfeld residuals for a
Cox proportional hazards model.
}
\usage{
\method{residuals}{coxph}(object,
       type=c("martingale", "deviance", "score", "schoenfeld",
	      "dfbeta", "dfbetas", "scaledsch","partial"),
       collapse=FALSE, weighted=FALSE, ...)
\method{residuals}{coxph.null}(object,
       type=c("martingale", "deviance","score","schoenfeld"),
       collapse=FALSE, weighted=FALSE, ...)
}
\arguments{
\item{object}{
an object inheriting from class \code{coxph}, representing a fitted Cox
regression model.
Typically this is the output from the \code{coxph} function.
}
\item{type}{
character string indicating the type of residual desired.
Possible values are \code{"martingale"}, \code{"deviance"}, \code{"score"}, \code{"schoenfeld"},
"dfbeta"', \code{"dfbetas"}, and \code{"scaledsch"}.
Only enough of the string to determine a unique match is required.
}
\item{collapse}{
vector indicating which rows to collapse (sum) over.
In time-dependent models more than one row data can pertain
to a single individual.
If there were 4 individuals represented by 3, 1, 2 and 4 rows of data
respectively, then \code{collapse=c(1,1,1, 2, 3,3, 4,4,4,4)} could be used to
obtain per subject rather than per observation residuals.
}
\item{weighted}{
if \code{TRUE} and the model was fit with case weights, then the weighted
residuals are returned.
}\item{...}{other unused arguments}}
\value{
For martingale and deviance residuals, the returned object is a vector
with one element for each subject (without \code{collapse}).
For score residuals it is a matrix
with one row per subject and one column per variable.
The row order will match the input data for the original fit.
For Schoenfeld residuals, the returned object is a matrix with one row
for each event and one column per variable.  The rows are ordered by time
within strata, and an attribute \code{strata} is attached that contains the
number of observations in each strata.
The scaled Schoenfeld residuals are used in the \code{cox.zph} function.


The score residuals are each individual's contribution to the score vector.
Two transformations of
this are often more useful: \code{dfbeta} is the approximate change in the
coefficient vector if that observation were dropped,
and \code{dfbetas} is the approximate change in the coefficients, scaled by
the standard error for the coefficients.
}
\section{NOTE}{
For deviance residuals, the status variable may need to be reconstructed.
For score and Schoenfeld residuals, the X matrix will need to be reconstructed.
}
\references{
T. Therneau, P. Grambsch, and T. Fleming. "Martingale based residuals
for survival models", \emph{Biometrika}, March 1990.
}
\seealso{
\code{\link{coxph}}}
\examples{

 fit <- coxph(Surv(start, stop, event) ~ (age + surgery)* transplant,
               data=heart)
 mresid <- resid(fit, collapse=heart$id)
}
\keyword{survival}
% Converted by Sd2Rd version 0.3-2.
