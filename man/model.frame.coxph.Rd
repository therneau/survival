\name{model.frame.coxph}
\Rdversion{1.1}
\alias{model.frame.coxph}
\title{Model.frame method for coxph objects}
\description{ Recreate the model frame of a coxph fit. }
\usage{
\method{model.frame}{coxph}(formula, ...)
}
\arguments{
  \item{formula}{the result of a \code{coxph} fit}
  \item{\dots}{other arguments to \code{model.frame}}
}
\details{
  For details, see the manual page for the generic function.
  This function would rarely be called by a user, it is mostly used
  inside functions like \code{residual} that need to recreate the data
  set from a model in order to do further calculations.
}
\value{the model frame used in the original fit, or a parallel one for
  new data.
}
\author{Terry Therneau}
\seealso{\code{\link{model.frame}}}
\keyword{ survival }
