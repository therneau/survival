\name{cox.zph}
\alias{cox.zph}
\alias{[.cox.zph}
\alias{print.cox.zph}
\title{
  Test the Proportional Hazards Assumption of a Cox Regression
}
\description{
Test the proportional hazards assumption for a Cox regression model fit 
(\code{coxph}). 
}
\usage{
cox.zph(fit, transform="km", terms=TRUE, singledf=FALSE, global=TRUE)
}
\arguments{
  \item{fit}{
    the result of fitting a Cox regression model, using the
    \code{coxph} or \code{coxme} functions. 
  }
  \item{transform}{
    a character string specifying how the survival times should be transformed 
    before the test is performed. 
    Possible values are \code{"km"}, \code{"rank"}, \code{"identity"} or a 
    function of one argument. 
  }
  \item{terms}{if TRUE, do a test for each term in the model rather than
    for each separate covariate.  For a factor variable with k levels,
    for instance, this would lead to a k-1 degree of freedom test.  The
    plot for such variables will be a single curve evaluating the linear
    predictor over time.}
  \item{singledf}{use a single degree of freedom test for terms that
    have multiple coefficients, i.e., the test that corresponds most
    closely to the plot.  If \code{terms=FALSE} this argument has no
    effect.}
  \item{global}{
    should a global chi-square test be done, in addition to the 
    per-variable or per-term tests tests. 
  }
}
\value{
an object of class \code{"cox.zph"}, with components: 

\item{table}{
a matrix with one row for each variable, and optionally a last row for 
the global test. 
Columns of the matrix contain a score test of for addition of the
time-dependent term, the degrees of freedom,
and the two-sided p-value. 
}
\item{x}{
the transformed time axis. 
}
\item{time}{the untransformed time values}
\item{y}{
the matrix of scaled Schoenfeld residuals.  There will be one column per 
term or per variable (depending on the \code{terms} option above),
and one row per event.  The row labels are a rounded form of the
original times.
}
\item{call}{
the calling sequence for the routine. 
}}

\details{
The computations require the original \code{x} matrix of the Cox model fit. 
Thus it saves time if the \code{x=TRUE} option is used in \code{coxph}. 
This function would usually be followed by both a plot and a print of 
the result. 
The plot gives an estimate of the time-dependent coefficient
\eqn{\beta(t)}{beta(t)}. 
If the proportional hazards assumption is true, the true
\eqn{\beta(t)}{beta(t)} function would be a  horizontal line.
The printout shows the score test for a time dependent effect of 0. 

Random effects terms such a \code{frailty} or random effects in a
\code{coxme} model are not checked for proportional hazards.
}
\note{
  In versions of the package before survival3.0 the function
  computed a fast approximation to the score test.  Later versions
  compute the actual score test.
}
\references{
P. Grambsch and T. Therneau (1994), 
Proportional hazards tests and diagnostics based on weighted residuals. 
\emph{Biometrika,}
\bold{81}, 515-26. 
}
\seealso{
\code{\link{coxph}},  \code{\link{Surv}}.   
}
\examples{
fit <- coxph(Surv(futime, fustat) ~ age + ecog.ps,  
             data=ovarian) 
temp <- cox.zph(fit) 
print(temp)                  # display the results 
plot(temp)                   # plot curves 
}
\keyword{survival}
