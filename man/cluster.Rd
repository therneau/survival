\name{cluster}
\alias{cluster}
\title{
  Identify clusters. 
}
\description{
This is a special function used in the context of survival models.  It 
identifies correlated groups of observations, and is used on the right hand 
side of a formula.
Using \code{cluster()} in a formula implies that
robust sandwich variance estimators are desired.}
\usage{
cluster(x)
}
\arguments{
  \item{x}{
    A character, factor, or numeric variable. 
  }
}
\value{
  \code{x} 
}
\details{
The function's only action is semantic, to mark a variable as the 
cluster indicator.
The resulting variance is what is known as the ``working independence''
variance  in a GEE model.
Note that one cannot use both a frailty term and a cluster term in the
same model, the first is a mixed-effects approach to correlation and the
second a GEE approach, and these don't mix.
}
\seealso{
\code{\link{coxph}},  \code{\link{survreg}}   
}
\examples{
marginal.model <- coxph(Surv(time, status) ~ rx + cluster(litter), rats)
frailty.model  <- coxph(Surv(time, status) ~ rx + frailty(litter), rats)
}
\keyword{survival}


