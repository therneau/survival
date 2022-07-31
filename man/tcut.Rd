\name{tcut}
\alias{tcut}
\alias{[.tcut}
\alias{levels.tcut}
\title{Factors for person-year calculations}
\description{
Attaches categories for person-year calculations to a variable without
losing the underlying continuous representation
}
\usage{
tcut(x, breaks, labels, scale=1)
\method{levels}{tcut}(x)
}
\arguments{
  \item{x}{numeric/date variable }
  \item{breaks}{breaks between categories, which are right-continuous }
  \item{labels}{labels for categories }
  \item{scale}{Multiply \code{x} and \code{breaks} by this.}
}

\value{
  An object of class \code{tcut}
}

\seealso{ \code{\link{cut}}, \code{\link{pyears}} }

\examples{
# For pyears, all time variable need to be on the same scale; but
# futime is in months and age is in years
test <- mgus2
test$years <- test$futime/30.5   # follow-up in years

# first grouping based on years from starting age (= current age)
# second based on years since enrollment (all start at 0)
test$agegrp <- tcut(test$age, c(0,60, 70, 80, 100), 
                     c("<=60", "60-70", "70-80", ">80"))
test$fgrp  <- tcut(rep(0, nrow(test)), c(0, 1, 5, 10, 100),
                   c("0-1yr", "1-5yr", "5-10yr", ">10yr"))

# death rates per 1000, by age group
pfit1 <- pyears(Surv(years, death) ~ agegrp, scale =1000, data=test)
round(pfit1$event/ pfit1$pyears) 

#death rates per 100, by follow-up year and age
# there are excess deaths in the first year, within each age stratum
pfit2 <- pyears(Surv(years, death) ~ fgrp + agegrp, scale =1000, data=test)
round(pfit2$event/ pfit2$pyears)  
}
\keyword{survival}

