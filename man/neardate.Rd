\name{neardate}
\alias{neardate}
\title{
Find the index of the closest date in list 2, for given dates in list 1.
}
\description{
  Given a set of observations with dates and an optional identifier, find
  the index in a second set of data, for each of the first, that is
  closest in date.
}
\usage{
neardate(date1, date2, id1, id2, best = c("after", "prior"),
nomatch = NA_integer_)
}

\arguments{
  \item{date1}{the dates to be matched.  This need not be a date: any
  data type with an ordered numeric representation is allowed.}
  
  \item{date2}{the dates to match against.  Must be of the same data
    type as the \code{date1} argument.}

  \item{id1}{optional vector of identifiers.  If present, matches are
    only between observations with the same identifier.}

  \item{id2}{optional vector of identifiers for the target.  Must be
    present if and only if \code{id1} is present.}

  \item{best}{find the first date2 value prior to or equal to the target
    date (date1), for first one equal to or after the target date?}
  \item{nomatch}{the value to return for items without a match}
}
\details{
  Closest date matching is often needed in clinical studies.  For
  example data set 1 might contain the subject identifier and the date
  of some procedure and data set set 2 has the dates and values for
  subsequent laboratory tests, and the query is to find the first
  cholesterol value after the intervention but no closer than 7 days.

  If \code{date1} and \code{date2} are not of the same class the user is
  on their own.  Since there are many case pairs where the result will
  be sensible the routine will proceed under an "if in doubt assume that
  the user knows what they are doing" assumption.  In particular expression
  \code{c(date1, date2)} must give a sensible result.
  Be careful about matching Date and DateTime values and the impact of
  time zones; see \code{\link{as.POSIXct}}. 
}
\value{the index of the matching observations in the second data set, or
  the \code{nomatch} value for no successful match}
\author{Terry Therneau}

\examples{
data1 <- data.frame(id=1:10,
                    date=as.Date(paste("2011", 1:10, "5", sep='-')))
temp1 <- c(1,4,5,1,3,6,9, 2,7,8,12,4,6,7,10,12,3)
data2 <- data.frame(id=c(1,1,1,2,2,4,4,5,5,5,6,8,8,9,10,10,12),
                    date=as.Date(paste("2011", temp1, "1", sep='-')),
                    chol = round(runif(17, 130, 280)))

#first cholesterol on or after enrollment
indx1 <- neardate(data1$date, data2$date, data1$id, data2$id)
data2[indx1, "chol"]

# closest one before or after
indx2 <- neardate(data1$date, data2$date, data1$id, data2$id,
                   best="prior")
ifelse(is.na(indx1), indx2,
       ifelse(is.na(indx2), indx1,
       ifelse(abs(data2$date[indx2]- data1$date) <
              abs(data2$date[indx1]-data1$date), indx2, indx1)))

# closest date before or after, but no more than 7 days prior to index
indx2 <- neardate(data1$date, data2$date-7, data1$id, data2$id,
                            best="after")
ifelse(is.na(indx1), indx2,
       ifelse(is.na(indx2), indx1,
       ifelse(abs(data2$date[indx2]- data1$date) <
              abs(data2$date[indx1]-data1$date), indx2, indx1)))

# find matching segments from a copy number reference
#  (example is a dummy reference)
set.seed(135)
temp <- sort(runif(10, 1, 5000))
reference <- data.frame(segment= 1:10,
                        start =  round(temp),
                        end = round(temp + runif(10, 100, 200)))
target <- c(200, 1500, 2300, 2200)
index <- neardate(target, reference$start)
ifelse(reference$end[index] >= target, index, NA)
}
\keyword{ manip }
\keyword{ utilities }
