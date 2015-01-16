\name{neardate}
\alias{neardate}
\title{
  Find the index of the closest value in data set 2, for each entry in
  data set one.
}
\description{
  A common task in medical work is to find the closest lab value to some
  index date, for each subject.
  This routine accomplishes that task. 
}
\usage{
neardate(id1, id2, y1, y2, best = c("after", "prior"),
nomatch = NA_integer_)
}

\arguments{
  \item{id1}{vector of subject identifiers for the index group}
  \item{id2}{vector of identifiers for the reference group}
  \item{y1}{normally a vector of dates for the index group,
    but any orderable data type is allowed}
  \item{y2}{reference set of dates}

  \item{best}{find the index of the first y2 value prior to or equal to
    the target y1 value, for each subject, or the one after?}

  \item{nomatch}{the value to return for items without a match}
}
\details{
  Closest date matching is often needed in clinical studies.  For
  example data set 1 might contain the subject identifier and the date
  of some procedure and data set set 2 has the dates and values for
  laboratory tests, and the query is to find the first
  test value after the intervention but no closer than 7 days.

  If \code{y1} and \code{y2} are not of the same class the user is
  on their own.  Since there are unmatched pairs where the result would
  be sensible the routine will proceed under an "if in doubt, assume that
  the user knows what they are doing" assumption.
  The expression
  \code{c(date1, date2)} must give a sensible, sortable result.
  Be careful about matching Date and DateTime values and the impact of
  time zones; see \code{\link{as.POSIXct}}. 
}
\value{the index of the matching observations in the second data set, or
  the \code{nomatch} value for no successful match}
\author{Terry Therneau}

\examples{
data1 <- data.frame(id = 1:10,
                    entry.dt = as.Date(paste("2011", 1:10, "5", sep='-')))
temp1 <- c(1,4,5,1,3,6,9, 2,7,8,12,4,6,7,10,12,3)
data2 <- data.frame(id = c(1,1,1,2,2,4,4,5,5,5,6,8,8,9,10,10,12),
                    lab.dt = as.Date(paste("2011", temp1, "1", sep='-')),
                    chol = round(runif(17, 130, 280)))

#first cholesterol on or after enrollment
indx1 <- neardate(data1$id, data2$id, data1$entry.dt, data2$lab.dt)
data2[indx1, "chol"]

# Closest one, either before or after. 
# 
indx2 <- neardate(data1$id, data2$id, data1$entry.dt, data2$lab.dt, 
                   best="prior")
ifelse(is.na(indx1), indx2, # none after, take before
       ifelse(is.na(indx2), indx1, #none before
       ifelse(abs(data2$date[indx2]- data1$date) <
              abs(data2$date[indx1]-data1$date), indx2, indx1)))

# closest date before or after, but no more than 7 days prior to index
indx2 <- neardate(data1$id, data2$id, data2$date, data2$date-7,
                            best="after")
ifelse(is.na(indx1), indx2,
       ifelse(is.na(indx2), indx1,
       ifelse(abs(data2$date[indx2]- data1$date) <
              abs(data2$date[indx1]-data1$date), indx2, indx1)))
}
\seealso{\code{\link{match}}, \code{\link{findInterval}}}
\keyword{ manip }
\keyword{ utilities }
