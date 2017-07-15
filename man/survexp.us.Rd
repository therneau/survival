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
      total United States population, by age and sex, 1940 to 2012.
    }
    \item{usr}{
      United States population, by age, sex and race, 1940 to 2012.
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
  tables from the US Census.  For age 25 in 1970, for instance,
  \eqn{p = 1-q} is is the
  probability that a subject who becomes 25 years of age in 1970 will
  achieve his/her 26th birthday.  The tables are recast in terms of
  hazard per day entirely for computational convenience.
  
  Each table is stored as an array, with additional attributes, and
  can be subset and manipulated as standard R arrays.  The attributes
  are
  \itemize{
    \item dimnames: a standard list of R dimnames.  For an categorical
    dimension such as \code{sex}, user's data can be numeric or match
    any unique subset of the dimnames, e.g. if the dimnames are
    male/female then "m" and "f" suffice.
    \item dimid: user variables are matched to dimensions by name.
    There can either be a \code{dimid} attribute containing those names,
    or the dimnames can be a named list.
    \item type: a vector with one element per dimension and values of
    1=numeric, 2= categorical, 3= date, and 4= date margin of a US
    ratetable.  In the \code{survexp.us} object the first dimension of age
    is numeric, the second of sex is categorical, and calendar year is of
    type 4.
    \item cutpoints: a list with one element per dimension.  Elements
    should be NULL for categorical dimensions, and otherwise they contain
    the starting value for each column.
    \item summary: a optional function that is called with the user's data
    and creates a summary of it, e.g.,
  "age ranges from 10.3 to 19.5 years".
}
   All numeric dimensions of a rate table must be in the same units.
   The \code{survexp.us} rate table contains daily hazard rates, the age
   cutpoints are in days, and the calendar year cutpoints are a Date.
}
\examples{
survexp.uswhite <- survexp.usr[,,"white",]
}
\keyword{survival}
\keyword{datasets}
