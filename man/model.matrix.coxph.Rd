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
\method{model.matrix}{coxph}(object, data, contrast.arg =
object$contrasts, mf, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{the result of a \code{coxph} model}
  \item{data}{optional, a data frame from which to obtain the data}
  \item{contrast.arg}{optional, a contrasts object describing how
    factors should be coded}
  \item{mf}{an optional model frame.  Used interally when programs have
    already obtained a model frame.}
  \item{\dots}{other possible argument to \code{model.frame}}
}
\details{
for details see the documentation of the generic method.
}
\value{
The model matrix for the fit
}
\author{Terry Therneau}
\seealso{\code{\link{model.matrix}}}
\examples{
fit1 <- coxph(Surv(time, status) ~ age + factor(ph.ecog), data=lung)
linear.predictor <- model.matrix(fit1) %*% coef(fit1)

fit2 <- coxph(Surv(time, status) ~ age + factor(ph.ecog), data=lung,
                                 x=TRUE)
all.equal(model.matrix(fit1), fit2$x)
}
\keyword{ survival }
