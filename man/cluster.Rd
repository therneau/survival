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
}
\seealso{
\code{\link{coxph}},  \code{\link{survreg}}   
}
\examples{
marginal.model <- coxph(Surv(time, status) ~ rx + cluster(litter), rats)
frailty.model  <- coxph(Surv(time, status) ~ rx + frailty(litter), rats)
}
\keyword{survival}


