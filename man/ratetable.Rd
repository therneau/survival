\name{ratetable}
\alias{ratetable}
\alias{[.ratetable}
\alias{[.ratetable2}
\alias{print.ratetable}
\alias{is.na.ratetable}
\alias{summary.ratetable}
\title{Ratetable reference in formula}
\description{
This function matches variable names in data to those in a ratetable for
\code{\link{survexp}}
}
\usage{
ratetable(...)
}

\arguments{
  \item{\dots}{tags matching dimensions of the ratetable and variables
    in the data frame (see example)}
}

\value{
  A data frame
}

\seealso{\code{\link{survexp}},\code{\link{survexp.us}},\code{\link{is.ratetable}}}

\examples{

fit <- survfit(Surv(time, status) ~ sex, pbc,subset=1:312)

# The data set does not have entry date, use the midpoint of the study
efit <- survexp(~ ratetable(sex=sex,age=age*365.35,year=as.Date('1979/1/1')) +
            sex, data=pbc,  times=(0:24)*182)
\dontrun{
plot(fit, mark.time=F, xscale=365.25, xlab="Years post diagnosis",
        ylab="Survival")
lines(efit, col=2, xscale=365.25) # Add the expected survival line
}
}
\keyword{survival}%-- one or more ...
