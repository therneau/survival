\name{coxphms.object}
\alias{coxphms.object}
\title{
  Multi-state Proportional Hazards Regression Object 
}
\description{
This class of objects is returned by the \code{coxph} class of functions 
to represent a fitted hazards model, when the model has
multiple states.  The object inherits from the \code{coxph} class.
}
\section{Components}{
  The object has all the components of a \code{coxph} object, with the
  following additions and variations.
}
\arguments{
  \item{states}{a character vector listing the states in the model}
\item{cmap}{the coefficient map. A matrix containing
  a column for each transition and a row for each coefficient, the value
  maps that transition/coefficient pair to a position in the coefficient
  vector.
  If a particular covariate is not used by a transition the matrix
  will contain a zero in that position, if two transitions share a
  coefficient the matrix will contain repeats.}
\item{smap}{the stratum map.
  The row labeled `(Baseline)' identifies transitions that do or do not
  share a baseline hazard.  Further rows correspond to strata() terms
  in the model, each of which may apply to some transitions and not others.
}
\item{rmap}{mapping for the residuals and linear predictors.  A two
  column matrix with one row for each element of the vectors
  and two columns, the first contains the data row and the second the
  transition.}
}
\details{
  In a multi-state model a set of intermediate observations is created
  during the computation, with a separate set of data rows for each
  transition.  An observation (id and time interval) that is at risk for
  more than one transition will for instance have a linear predictor and
  residual for each of the potential transitions.  As a result the vector
  of linear predictors will be longer than the number of observations.
  The \code{rmap} matrix shows the mapping.
  }
\seealso{
\code{\link{coxph}},  \code{\link{coxph.object}}
}
\keyword{survival}
