\name{coxph.detail}
\alias{coxph.detail}
\title{
Details of a Cox Model Fit
}
\description{
Returns the individual contributions to the first and second derivative
matrix, at each unique event time.
}
\usage{
coxph.detail(object, riskmat=FALSE)
}
\arguments{
 \item{object}{
   a Cox model object, i.e., the result of \code{coxph}.
 }
 \item{riskmat}{
   include the at-risk indicator matrix in the output?
 }
}
\value{
   a list with components
   \item{time}{
   the vector of unique event times
   }
    \item{nevent}{
    the number of events at each of these time points.
    }
    \item{means}{
    a matrix with one row for each event time and one column for each variable
    in the Cox model, containing the weighted mean of the variable at that time,
    over all subjects still at risk at that time.  The weights are the risk
    weights \code{exp(x \%*\% fit$coef)}.
    }
    \item{nrisk}{
    number of subjects at risk.
    }
   \item{score}{
    the contribution to the score vector (first derivative of the log
    partial likelihood) at each time point.
    }
    \item{imat}{
    the contribution to the information matrix (second derivative of the
    log partial likelihood) at each time point.
    }
    \item{hazard}{
    the hazard increment.  Note that the hazard and variance of the
    hazard are always for some particular future subject.  This routine
    uses \code{object$mean} as the future subject.
    }
     \item{varhaz}{
    the variance of the hazard increment.
    }
    \item{x,y}{
    copies of the input data.
    }
    \item{strata}{
    only present for a stratified Cox model, this is
    a table giving the number of time points of component \code{time} that
    were contributed by each of the strata.
    }
    \item{riskmat}{
      a matrix with one row for each time and one column for each observation
    containing a 0/1 value to indicate whether that observation was (1) or
    was not (0) at risk at the given time point.
  }
}
\details{
This function may be useful for those who wish to investigate new methods or
extensions to the Cox model.  The example below shows one way to calculate
the Schoenfeld residuals.
}
\seealso{
  \code{\link{coxph}}, \code{\link{residuals.coxph}}
  }
\examples{
fit   <- coxph(Surv(futime,fustat) ~ age + rx + ecog.ps, ovarian, x=TRUE)
fitd  <- coxph.detail(fit)
#  There is one Schoenfeld residual for each unique death.  It is a
# vector (covariates for the subject who died) - (weighted mean covariate
# vector at that time).  The weighted mean is defined over the subjects
# still at risk, with exp(X beta) as the weight.

events <- fit$y[,2]==1
etime  <- fit$y[events,1]   #the event times --- may have duplicates
indx   <- match(etime, fitd$time)
schoen <- fit$x[events,] - fitd$means[indx,]
}
\keyword{survival}

