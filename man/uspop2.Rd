\name{uspop2}
\alias{uspop2}
\docType{data}
\title{Projected US Population}
\description{US population by age and sex, for 2000 through 2020}
\usage{data(uspop2)}
\format{
  The data is a matrix with dimensions age, sex, and calendar year.
  Age goes from 0 through 100, where the value for age 100 is the total
  for all ages of 100 or greater.
}
\details{ This data is often used as a "standardized" population for
  epidemiolgy studies.}

\source{
NP2008_D1:  Projected Population by Single Year of Age, Sex, Race, and Hispanic 
Origin for the United States: July 1, 2000 to July 1, 2050,
www.census.gov/population/projections.
}
\examples{
us50 <- uspop2[51:101,, "2000"]  #US 2000 population, 50 and over
age <- as.integer(dimnames(us50)[[1]])
smat <- model.matrix( ~ factor(floor(age/5)) -1)
ustot <- t(smat) \%*\% us50  #totals by 5 year age groups
temp <- c(50,55, 60, 65, 70, 75, 80, 85, 90, 95)
dimnames(ustot) <- list(c(paste(temp, temp+4, sep="-"), "100+"),
                         c("male", "female"))
}
\seealso{\code{\link{uspop}}}
\keyword{datasets}
