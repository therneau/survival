\name{clogit}
\alias{clogit}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Conditional logistic regression }
\description{
  Estimates a logistic regression model by maximising the conditional
  likelihood. Uses a model formula of the form
  \code{case.status~exposure+strata(matched.set)}.
  The default is to use the exact conditional likelihood, a commonly
  used approximate conditional likelihood is provided for compatibility
  with older software.
}
\usage{
clogit(formula, data, method=c("exact", "approximate"), na.action=getOption("na.action"), subset=NULL,control=coxph.control())
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{formula}{Model formula}
  \item{data}{data frame }
  \item{method}{correct or approximate conditonal likelihood?}
  \item{na.action}{missing value handling}
  \item{subset}{subset of records to use }
  \item{control}{control values}
}

\value{
  An object of class \code{"clogit"}, which is a wrapper for a
  \code{"coxph"} object.
}
\author{Thomas Lumley}


\seealso{\code{\link{strata}},\code{\link{coxph}},\code{\link{glm}} }

\examples{
\dontrun{clogit(case~spontaneous+induced+strata(stratum),data=infert)}
}
\keyword{survival}%-- one or more ...
\keyword{models}
