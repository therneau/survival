\name{survfit.matrix}
\alias{survfit.matrix}
\title{Create Aalen-Johansen estimates of multi-state survival from
  a matrix of hazards.}
\description{
  This allows one to create the Aalen-Johansen estimate of P, a matrix
  with one column per state and one row per time, starting with the
  individual hazard estimates.  Each row of P will sum to 1.
  Note that this routine has been superseded by the use of multi-state
  Cox models, and will eventually be removed.
}
\usage{
\method{survfit}{matrix}(formula, p0, method = c("discrete", "matexp"),
        start.time, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a matrix of lists, each element of which is either NULL or
    a survival curve object. 
}
  \item{p0}{the initial state vector.  The names of this vector are used
  as the names of the states in the output object.  If there are
  multiple curves then \code{p0} can be a matrix with one row per curve.
}
  \item{method}{
    use a product of discrete hazards, or a product of matrix exponentials.
    See details below.
  }
  \item{start.time}{optional; start the calculations at a given starting
    point}
  \item{...}{further arguments used by other survfit methods}
}
\details{
  On input the matrix should contain a set of predicted curves for each
  possible transition, and NULL in other positions.  Each of the
  predictions will have been obtained from the relevant Cox model. 
  This approach for multistate curves is easy to use but has some
  caveats.  
  First, the input curves must be consistent.
  The routine checks as best it can, but can easy be fooled.
  For instance, if one were to fit two Cox models, obtain
  predictions for males and females from one, and for treatment A and B
  from the other, this routine will create two curves but they are
  not meaningful.
  A second issue is that standard errors are not produced.

  The names of the resulting states are taken from the names of the
  vector of initial state probabilities.  If they are missing, then the
  dimnames of the input matrix are used, and lacking that the labels
  '1', '2', etc. are used.

  For the usual Aalen-Johansen estimator the multiplier at each event
  time is the matrix of hazards H (also written as I + dA).
  When using predicted survival curves from a Cox model, however,
  it is possible to get predicted hazards that are greater than 1, which
  leads to probabilities less than 0.
  If the \code{method} argument is not supplied and the input
  curves are derived from a Cox model this routine instead uses the
  approximation expm(H-I) as the multiplier, which always gives valid
  probabilities.
  (This is also the standard approach for ordinary
  survival curves from a Cox model.) 
}
\value{a survfitms object}
\author{Terry Therneau}

\note{The R syntax for creating a matrix of lists is very fussy.}

\seealso{\code{\link{survfit}}}

\examples{
etime <- with(mgus2, ifelse(pstat==0, futime, ptime))
event <- with(mgus2, ifelse(pstat==0, 2*death, 1))
event <- factor(event, 0:2, labels=c("censor", "pcm", "death"))

cfit1 <- coxph(Surv(etime, event=="pcm") ~ age + sex, mgus2)
cfit2 <- coxph(Surv(etime, event=="death") ~ age + sex, mgus2)

# predicted competing risk curves for a 72 year old with mspike of 1.2
# (median values), male and female.
# The survfit call is a bit faster without standard errors.
newdata <- expand.grid(sex=c("F", "M"), age=72, mspike=1.2)

AJmat <- matrix(list(), 3,3)
AJmat[1,2] <- list(survfit(cfit1, newdata, std.err=FALSE))
AJmat[1,3] <- list(survfit(cfit2, newdata, std.err=FALSE))
csurv  <- survfit(AJmat, p0 =c(entry=1, PCM=0, death=0))
}

\keyword{survival }
