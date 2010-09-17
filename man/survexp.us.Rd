\name{ratetables}
\alias{survexp.us}
\alias{survexp.usr}
\alias{survexp.az}
\alias{survexp.azr}
\alias{survexp.fl}
\alias{survexp.flr}
\alias{survexp.mn}
\alias{survexp.mnwhite}
\alias{survexp.wnc}
\title{
  Census Data Sets for the Expected Survival and Person Years Functions
}
\description{
  Census data sets for the expected survival and person years functions.
}
\details{
  \describe{
    \item{us}{
      total United States population, by age and sex, 1960 to 1980.
    }
    \item{uswhite}{
      United States white population, by age and sex, 1950 to 1980. This
      is no longer included, but can be extracted from
      \code{survexp.usr} as shown in the examples.
    }
    \item{usr}{
      United States population, by age, sex and race, 1960 to 1980.
      Race is white, nonwhite, or black.  For 1960 and 1970 the black
      population values were not reported separately, so the nonwhite
      values were used.
    }
    \item{mn}{
      total Minnesota population, by age and sex, 1970 and 1980.
    }
    \item{mnwhite}{
      Minnesota white population, by age and sex, 1960 to 1980.
    }
    \item{fl}{
      total Florida population, by age and sex, 1970 and 1980.
    }
    \item{flr}{
      Florida population, by age, sex and race, 1970-1980.  Race is
      white, nonwhite, or black.  For 1970 the black population values
      were not reported separately, so the nonwhite values were used.
    }
    \item{az}{
      total Arizona population, by age and sex, 1970 and 1980.
    }
    \item{azr}{
      Arizona population, by age, sex and race, 1970-1980.  Race is
      white versus nonwhite.  For 1970 the nonwhite population values
      were not reported separately.
      In order to make the rate table be a matrix, the 1980 values were
      repeated.  (White and non-white values are quite different).
    }
  }
  
  Each of these tables contains the daily hazard rate for a matched
  subject from the population, defined as \eqn{-\log(1-q)/365.24} where
  \eqn{q} is the 1 year probability of death as reported in the original
  tables.  For age 25 in 1970, for instance, \eqn{p = 1-q} is is the
  probability that a subject who becomes 25 years of age in 1970 will
  achieve his/her 26th birthday.  The tables are recast in terms of
  hazard per day entirely for computational convenience.  (The fraction
  .24 in the denominator is based on 24 leap years per century.)

  Each table is stored as an array, with additional attributes, and
  can be subset and manipulated as standard S arrays.  Interpolation
  between calendar years is done \dQuote{on the fly} by the
  \code{\link{survexp}} routine. 

  Some of the deficiencies, e.g., 1970 Arizona non-white, are a result
  of local (Mayo Clinic) conditions.  The data probably exists, but we
  don't have a copy it in the library.

  The tables have been augmented to contain extrapolated values for 1990
  and 2000.  The details can be found in Mayo Clinic Biostatistics
  technical report 63 at \url{http://www.mayo.edu/hsr/techrpt.html}.
}
\examples{
survexp.uswhite <- survexp.usr[,,"white",]
}
\keyword{survival}
\keyword{datasets}
