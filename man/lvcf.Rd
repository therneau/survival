\name{lvcf}
\alias{lvcf}
\title{
  Last Value Carried Forward
}
\description{
  Replace missing values in a covariate with the last non-missing value,
  separately for each id.
}
\usage{
lvcf(id, x, time)
}
\arguments{
  \item{id}{subject identifier}
  \item{x}{the covariate which will be modified}
  \item{time}{optional; used to sort observations within subject}
}
\details{
If an analysis includes a covariate that was measured repeatedly,
cholesterol say, the data may often include visits where this was not
measured. It is convenient in that case to use, for each subject,
the most recently measured value.
When a value is extended beyond the last measurement, rather that just
replacement of intermediate values, there is the potential for bias,
e.g., if the reason for discontinuation is predictive of the future
trajectory of values.

If \code{time} is missing, then the order of the values within subject
is assumed to be the calendar time order.  It is not necessary for the
the data to be sorted by \code{id}.
}
\value{
an updated copy of \code{x}
}
\author{Terry Therneau}
\note{
Other R packages also implement an locf function, e.g., \code{na.locf}
in the \code{zoo} package or \code{LOCF} in \code{DescTools}.
However, they do not include the \code{id} argument, so per subjct application
requires further processing.
The tdc operation in \code{tmerge} also performs lvcf substitution.

We have used the phrase "last value" rather than "last observation"
since this function only addresses a single variable with each call, not
the entire row of a data frame (observation) as some others do.
}

\seealso{
  \code{\link{tmerge}}
}
\examples{
newplat <- with(pbcseq, lvcf(id, platelet))
table(is.na(pbcseq$platelet), is.na(newplat))
}
\keyword{ survival }
