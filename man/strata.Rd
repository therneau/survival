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
strata(..., na.group=FALSE, shortlabel, sep=', ')
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
  factor labels.  The default action is to omit the names if all of the
  arguments are factors, and none of them was named.}
\item{sep}{
the character used to separate groups, in the created label
}
}
\value{
a new factor, whose levels are all possible combinations of the factors 
supplied as arguments. 
}
\details{
When used outside of a \code{coxph} formula the result of the function
is essentially identical to the \code{interaction} function, 
though the labels from \code{strata} are often more verbose. 
}
\seealso{
\code{\link{coxph}}, \code{\link{interaction}}
}
\examples{
a <- factor(rep(1:3,4), labels=c("low", "medium", "high"))
b <- factor(rep(1:4,3))
levels(strata(b))
levels(strata(a,b,shortlabel=TRUE))

coxph(Surv(futime, fustat) ~ age + strata(rx), data=ovarian) 
}
\keyword{survival}
