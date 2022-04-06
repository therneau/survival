\name{levels.Surv}
\alias{levels.Surv}
\title{Return the states of a multi-state Surv object
}
\description{
  For a multi-state \code{Surv} object, this will return the names
  of the states.  
}
\usage{
\method{levels}{Surv}(x)
}
\arguments{
  \item{x}{a \code{Surv} object}
}
\value{
  for a multi-state \code{Surv} object, the vector of state names
  (excluding censoring); or NULL for an ordinary \code{Surv} object
}
\examples{
y1 <- Surv(c(1,5, 9, 17,21, 30),
           factor(c(0, 1, 2,1,0,2), 0:2, c("censored", "progression", "death")))
levels(y1)

y2 <- Surv(1:6, rep(0:1, 3))
y2
levels(y2)
}
\keyword{ survival }
