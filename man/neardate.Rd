\name{neardate}
\alias{neardate}
\title{
  Find the index of the closest value in data set 2, for each entry in
  data set one.
}
\description{
  A common task in medical work is to find the closest lab value to some
  index date, for each subject.
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
  This routine is closely related to \code{match} and to
  \code{findInterval}, the first of which finds exact matches and
  the second closest matches.  This finds the closest matching date
  within sets of exactly matching identifiers.  
  Closest date matching is often needed in clinical studies.  For
  example data set 1 might contain the subject identifier and the date
  of some procedure and data set set 2 has the dates and values for
  laboratory tests, and the query is to find the first
  test value after the intervention but no closer than 7 days.

  The \code{id1} and \code{id2} arguments are similar to \code{match} in
  that we are searching for instances of \code{id1} that will be found
  in \code{id2}, and the result is the same length as \code{id1}.
  However, instead of returning the first match with \code{id2} this
  routine returns the one that best matches with respect to \code{y1}.

  The \code{y1} and \code{y2} arguments need not be dates, the function
  works for any data type such that the expression
  \code{c(y1, y2)} gives a sensible, sortable result.
  Be careful about matching Date and DateTime values and the impact of
  time zones, however, see \code{\link{as.POSIXct}}. 
  If \code{y1} and \code{y2} are not of the same class the user is
  on their own.
  Since there exist pairs of unmatched data types where the result could
  be sensible, the routine will in this case proceed under the assumption
  that "the user knows what they are doing".  Caveat emptor.
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
       ifelse(abs(data2$lab.dt[indx2]- data1$entry.dt) <
              abs(data2$lab.dt[indx1]- data1$entry.dt), indx2, indx1)))

# closest date before or after, but no more than 21 days prior to index
indx2 <- ifelse((data1$entry.dt - data2$lab.dt[indx2]) >21, NA, indx2)
ifelse(is.na(indx1), indx2, # none after, take before
       ifelse(is.na(indx2), indx1, #none before
       ifelse(abs(data2$lab.dt[indx2]- data1$entry.dt) <
              abs(data2$lab.dt[indx1]- data1$entry.dt), indx2, indx1)))
}
\seealso{\code{\link{match}}, \code{\link{findInterval}}}
\keyword{ manip }
\keyword{ utilities }
