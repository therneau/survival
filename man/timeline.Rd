\name{timeline}
\alias{timeline}
\alias{totimeline}
\alias{fromtimeline}
\title{
Convert to/from a timeline data set format}
\description{
Convert from a 'timeline' data set format for survival data to the
counting process form, and vice versa.
}
\usage{
totimeline(formula, data, id, istate)
fromtimeline(formula, data, id, istate="istate")
}
\arguments{
  \item{formula}{a model formula with a \code{Surv} object on the left}
  \item{data}{data set in which to evaluate the formula}
  \item{id}{the name of the identifier variable, which will be searched
    first in the \code{data}.  Multiple rows for the same subject will
    share the same \code{id} value.}
  \item{istate}{for \code{totimeline} the name of the variable in the
    counting process data set that contains the initial state.  For
    \code{fromtimeline} the variable name to use for the initial state.}
}
\details{
Counting process style data sets are heavily used in the \code{survival}
package for both time-dependent covariates and multistate data.  Each
row of the data will contain a time interval (t1, t2), status or state
at the end of the interval, covariate values that apply over the
interval, and an id variable.
A timeline data set will have a single time covariate, an \code{id}
variable, along with other covariate and outcome values that were
observed at that time point. If some covariates are observed at a
particular time point but others were not, these other values would be
missing for that row.  (The exception are covariates that are constant,
like birthdate or a genetic marker, which will normally appear across
all rows).

A disadvantage of the counting process form is that it requires special
tools for manipulation, e.g., \code{tmerge}; timeline data sets are much
simpler in structure and thus can benefit from a much wider variety of
tools in their creation.  They are also more direct wrt ensuring
validity: each row should encode what was \emph{actually observed} at
that time point.  Another potential advantage is for variables such as
diabetes, which might be used as an outcome in one model and a predictor
in another.  This requires two separate variables in a counting process
data set, since covariates change at the beginning of a time interval
and outcomes happen at the end of it.

The conversion from timeline to counting process form uses the same
rules with respect to missing values as \code{tmerge}, it is in fact
what is used behind the scenes to do the conversion.
}
\value{a data set of the proper form}
\note{
  This is at present an experimental feature.
}
\seealso{\code{\link{tmerge}}, \code{\link{survSplit}}
}
\keyword{ survival }

