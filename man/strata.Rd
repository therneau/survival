\name{strata}
\alias{strata}
\title{
  Identify Stratification Variables 
}
\description{
This is a special function used in the context of the Cox survival model. 
It identifies stratification variables when they appear on the right hand 
side of a formula. 
}
\usage{
strata(..., na.group=FALSE, shortlabel=FALSE, sep=', ')
}
\arguments{
\item{\dots}{
any number of variables.  All must be the same length. 
}
\item{na.group}{
a logical variable, if \code{TRUE}, then missing values are treated as a 
distinct level of each variable. 
}
\item{shortlabel}{if \code{TRUE} omit variable names from resulting
  factor labels}
\item{sep}{
the character used to separate groups, in the created label
}
}
\value{
a new factor, whose levels are all possible combinations of the factors 
supplied as arguments. 
}
\details{
The result is identical to the \code{interaction} function, 
but for the labeling of the factors (\code{strata} is more verbose). 
}
\seealso{
\code{\link{coxph}}, \code{\link{interaction}}
}
\examples{
a<-factor(rep(1:3,4))
b<-factor(rep(1:4,3))
levels(strata(a))
levels(strata(a,b,shortlabel=TRUE))

coxph(Surv(futime, fustat) ~ age + strata(rx), data=ovarian) 
}
\keyword{survival}
