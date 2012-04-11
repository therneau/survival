\name{coxph}
\alias{coxph}
\alias{vcov.coxph}
\alias{print.coxph.null}
\alias{print.coxph.penal}
\alias{coxph.penalty}
\alias{[.coxph.penalty}
\alias{coxph.getdata}
\alias{summary.coxph.penal}
\title{
  Fit Proportional Hazards Regression Model 
}
\description{
Fits a Cox proportional hazards regression model. 
Time dependent variables, time dependent strata, multiple events per subject, 
and other extensions are incorporated using the counting process formulation 
of Andersen and Gill. 
}
\usage{
coxph(formula, data=, weights, subset, 
      na.action, init, control, 
      ties=c("efron","breslow","exact"), 
      singular.ok=TRUE, robust=FALSE, 
      model=FALSE, x=FALSE, y=TRUE, tt, method, ...)
}
\arguments{
  \item{formula}{
    a formula object, with the response on the left of a \code{~} operator, and 
    the terms on the right.  The response must be a survival object as 
    returned by the \code{Surv} function. 
  }
  \item{data}{
    a data.frame in which to interpret the variables named in 
    the \code{formula}, or in the \code{subset} and the \code{weights}
    argument. 
  }
  \item{weights}{
    vector of case weights.  If \code{weights} is a vector of integers,
    then the 
    estimated coefficients are equivalent to estimating the model from data 
    with the individual cases replicated as many times as indicated by 
    \code{weights}. 
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
  \item{init}{
    vector of initial values of the iteration.  Default initial 
    value is zero for all variables. 
  }
  \item{control}{
    Object of class \code{\link{coxph.control}} specifying iteration limit
    and other control options. Default is \code{coxph.control(...)}.
  }
  \item{ties}{
    a character string specifying the method for tie handling.  If there  
    are no tied death times all the methods are equivalent. 
    Nearly all Cox regression programs use the Breslow method by default, 
    but not this one. 
    The Efron approximation is used as the default here, it is more 
    accurate when dealing with tied death times, and is as efficient 
    computationally. 
    The ``exact partial likelihood'' is 
    equivalent to a conditional logistic model, and is appropriate when
    the times are a small set of discrete values.
    If there are a large number of ties and (start, stop) style survival
    data the computational time will be excessive. 
  }
  \item{singular.ok}{
    logical value indicating how to handle collinearity in the model matrix. 
    If \code{TRUE}, the program will automatically skip over columns of the X 
    matrix that are linear combinations of earlier columns.  In this case the 
    coefficients for such columns will be NA, and the variance matrix will 
    contain zeros. For ancillary calculations, such as the linear
    predictor, 
    the missing coefficients are treated as zeros. 
  }
  \item{robust}{
    this argument has been deprecated, use a cluster term in the model
    instead.
  }
  \item{model}{
    logical value: if \code{TRUE}, the model frame is returned in component
    \code{model}. 
  }
  \item{x}{
    logical value: if \code{TRUE}, the x matrix is returned in
    component \code{x}. 
  }
  \item{y}{
    logical value: if \code{TRUE}, the response vector is returned in
    component \code{y}. 
  }
  \item{tt}{optional list of time-transform functions.}
  \item{method}{alternate name for the \code{ties} argument.}
  \item{...}{Other arguments will be passed to \code{\link{coxph.control}}
  }
}
\value{
an object of class \code{coxph} representing the fit. 
See \code{coxph.object} for details.  
}
\section{Side Effects}{
  Depending on the call, the \code{predict}, \code{residuals},
and \code{survfit} routines may 
need to reconstruct the x matrix created by \code{coxph}.
It is possible for this to fail, as in the example below in
which the predict function is unable to find \code{tform}.

\preformatted{  tfun <- function(tform) coxph(tform, data=lung)
  fit <- tfun(Surv(time, status) ~ age)
  predict(fit)}

In such a case add the \code{model=TRUE} option to the
\code{coxph} call to obviate the
need for reconstruction, at the expense of a larger \code{fit}
object.  
}
\details{
The proportional hazards model is usually expressed in terms of a 
single survival time value for each person, with possible censoring. 
Andersen and Gill reformulated the same problem as a counting process; 
as time marches onward we observe the events for a subject, rather 
like watching a Geiger counter. 
The data for a subject is presented as multiple rows or "observations", 
each 
of which applies to an interval of observation (start, stop].

The routine internally scales and centers data to avoid overflow in
the argument to the exponential function.  These actions do not change
the result, but lead to more numerical stability.
However, arguments to offset are not scaled since there are situations
where a large offset value is a purposefully used.
Users should not use normally allow large numeric offset values.
}
\section{Special terms}{
There are three special terms that may be used in the model equation. 
A \code{strata} term identifies a stratified Cox model; separate baseline 
hazard functions are fit for each strata. 
The \code{cluster} term is used to compute a robust variance for the model. 
The term \code{+ cluster(id)} where each value of \code{id} is unique is equivalent to 
specifying the \code{robust=T} argument, and produces an approximate 
jackknife estimate of the variance.  If the \code{id} variable were not 
unique, but instead 
identifies clusters of correlated observations, then the variance 
estimate is based on a grouped jackknife.

A time-transform term allows variables to vary dynamically in time.  In
this case the \code{tt} argument will be a function or a list of
functions (if there are more than one tt() term in the model) giving the
appropriate transform.   See the examples below.
}
\section{Convergence}{
In certain data cases the actual MLE estimate of a 
coefficient is infinity, e.g., a dichotomous variable where one of the 
groups has no events.  When this happens the associated coefficient 
grows at a steady pace and a race condition will exist in the fitting 
routine: either the log likelihood converges, the information matrix 
becomes effectively singular, an argument to exp becomes too large for 
the computer hardware, or the maximum number of interactions is
exceeded.
(Nearly always the first occurs.)
The routine attempts to detect when this has happened,
not always successfully.
The primary consequence for he user is that the Wald statistic =
coefficient/se(coefficient) is not valid in this case and should be
ignored; the likelihood ratio and score tests remain valid however.
}
\section{Penalized regression}{
  \code{coxph} can now maximise a penalised partial likelihood with
arbitrary user-defined penalty.  Supplied penalty functions include
ridge regression (\link{ridge}), smoothing splines
(\link{pspline}), and frailty models (\link{frailty}).
} 
\references{
Andersen, P. and Gill, R. (1982). 
Cox's regression model for counting processes, a large sample study. 
\emph{Annals of Statistics}
\bold{10}, 1100-1120. 

Therneau, T., Grambsch, P., Modeling Survival Data: Extending the Cox Model. 
Springer-Verlag, 2000.
}
\seealso{
\code{\link{cluster}},  \code{\link{strata}},  \code{\link{Surv}},
\code{\link{survfit}}, \code{\link{pspline}}, \code{\link{frailty}},
\code{\link{ridge}}.   
}
\examples{
# Create the simplest test data set 
test1 <- list(time=c(4,3,1,1,2,2,3), 
              status=c(1,1,1,0,1,1,0), 
              x=c(0,2,1,1,1,0,0), 
              sex=c(0,0,0,0,1,1,1)) 
# Fit a stratified model 
coxph(Surv(time, status) ~ x + strata(sex), test1) 
# Create a simple data set for a time-dependent model 
test2 <- list(start=c(1,2,5,2,1,7,3,4,8,8), 
              stop=c(2,3,6,7,8,9,9,9,14,17), 
              event=c(1,1,1,1,1,1,1,0,0,0), 
              x=c(1,0,0,1,0,1,1,1,0,0)) 
summary(coxph(Surv(start, stop, event) ~ x, test2)) 

#
# Create a simple data set for a time-dependent model
#
test2 <- list(start=c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
                stop =c(2, 3, 6, 7, 8, 9, 9, 9,14,17),
                event=c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
                x    =c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0) )


summary( coxph( Surv(start, stop, event) ~ x, test2))

# Fit a stratified model, clustered on patients 

bladder1 <- bladder[bladder$enum < 5, ] 
coxph(Surv(stop, event) ~ (rx + size + number) * strata(enum) + 
      cluster(id), bladder1)

# Fit a time transform model using current age
coxph(Surv(time, status) ~ ph.ecog + tt(age), data=lung,
     tt=function(x,t,...) pspline(x + t/365.25))
}
\keyword{survival}

