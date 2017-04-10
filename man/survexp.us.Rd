\name{ratetables}
\alias{survexp.us}
\alias{survexp.usr}
\alias{survexp.mn}
\title{
  Census Data Sets for the Expected Survival and Person Years Functions
}
\description{
  Census data sets for the expected survival and person years functions.
}
\details{
  \describe{
    \item{us}{
      total United States population, by age and sex, 1940 to 2004.
    }
    \item{usr}{
      United States population, by age, sex and race, 1940 to 2004.
      Race is white, nonwhite, or black.  For 1960 and 1970 the black
      population values were not reported separately, so the nonwhite
      values were used.
    }
    \item{mn}{
      total Minnesota population, by age and sex, 1970 to 2004.
    }
   }
  
  Each of these tables contains the daily hazard rate for a matched
  subject from the population, defined as \eqn{-\log(1-q)/365.25} where
  \eqn{q} is the 1 year probability of death as reported in the original
  tables.  For age 25 in 1970, for instance, \eqn{p = 1-q} is is the
  probability that a subject who becomes 25 years of age in 1970 will
  achieve his/her 26th birthday.  The tables are recast in terms of
  hazard per day entirely for computational convenience.
  
  Each table is stored as an array, with additional attributes, and
  can be subset and manipulated as standard R arrays. 
}
\examples{
survexp.uswhite <- survexp.usr[,,"white",]
}
\keyword{survival}
\keyword{datasets}
