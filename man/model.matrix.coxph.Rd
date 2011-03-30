\name{model.matrix.coxph}
\Rdversion{1.1}
\alias{model.matrix.coxph}
\title{
Model.matrix method for coxph models
}
\description{
Reconstruct the model matrix for a cox model.
}
\usage{
\method{model.matrix}{coxph}(object, data=NULL, contrast.arg =
 object$contrasts, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{the result of a \code{coxph} model}
  \item{data}{optional, a data frame from which to obtain the data}
  \item{contrast.arg}{optional, a contrasts object describing how
    factors should be coded}
  \item{\dots}{other possible argument to \code{model.frame}}
}
\details{
  When there is a \code{data} argument this function differs from most
  of the other \code{model.matrix} methods in that the response variable
  for the original formual is \emph{not} required to be in the data.

  If the data frame contains a \code{terms} attribute then it is
  assumed to be the result of a call to \code{model.frame}, otherwise
  a call to \code{model.frame} is applied with the data as an argument.
}
\value{
The model matrix for the fit
}
\author{Terry Therneau}
\seealso{\code{\link{model.matrix}}}
\examples{
fit1 <- coxph(Surv(time, status) ~ age + factor(ph.ecog), data=lung)
xfit <- model.matrix(fit1)

fit2 <- coxph(Surv(time, status) ~ age + factor(ph.ecog), data=lung,
                                 x=TRUE)
all.equal(model.matrix(fit1), fit2$x)
}
\keyword{ survival }
