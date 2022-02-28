\name{Surv-methods}
\alias{Math.Surv}
\alias{Ops.Surv}
\alias{Summary.Surv}
\alias{anyDuplicated.Surv}
\alias{as.character.Surv}
\alias{as.data.frame.Surv}
\alias{as.integer.Surv}
\alias{as.matrix.Surv}
\alias{as.numeric.Surv}
\alias{c.Surv}
\alias{duplicated.Surv}
\alias{format.Surv}
\alias{head.Surv}
\alias{is.na.Surv}
\alias{length.Surv}
\alias{mean.Surv}
\alias{median.Surv}
\alias{names.Surv}
\alias{names<-.Surv}
\alias{quantile.Surv}
\alias{plot.Surv}
\alias{rep.Surv}
\alias{rep.int.Surv}
\alias{rep_len.Surv}
\alias{rev.Surv}
\alias{t.Surv}
\alias{tail.Surv}
\alias{unique.Surv}
\title{Methods for Surv objects}
\description{The list of methods that apply to \code{Surv} objects}
\usage{
    \method{anyDuplicated}{Surv}(x, ...)
    \method{as.character}{Surv}(x, ...)
    \method{as.data.frame}{Surv}(x, ...)
    \method{as.integer}{Surv}(x, ...)
    \method{as.matrix}{Surv}(x, ...)
    \method{as.numeric}{Surv}(x, ...)
    \method{c}{Surv}(...)
    \method{duplicated}{Surv}(x, ...)
    \method{format}{Surv}(x, ...)
    \method{head}{Surv}(x, ...)
    \method{is.na}{Surv}(x)
    \method{length}{Surv}(x)
    \method{mean}{Surv}(x, ...)
    \method{median}{Surv}(x, na.rm=FALSE, ...)
    \method{names}{Surv}(x)
    \method{names}{Surv}(x) <- value
    \method{quantile}{Surv}(x, probs, na.rm=FALSE, ...)
    \method{plot}{Surv}(x, ...)
    \method{rep}{Surv}(x, ...)
    \method{rep.int}{Surv}(x, ...)
    \method{rep_len}{Surv}(x, ...)
    \method{rev}{Surv}(x)
    \method{t}{Surv}(x)
    \method{tail}{Surv}(x, ...)
    \method{unique}{Surv}(x, ...)
}
\arguments{
  \item{x}{a \code{Surv} object}
  \item{probs}{a vector of probabilities}
  \item{na.rm}{remove missing values from the calculation}
  \item{value}{a character vector of up to the same length as \code{x}, or
    \code{NULL}}
  \item{\ldots}{other arguments to the method}
}
\details{
   These functions extend the standard methods to \code{Surv} objects.
   The arguments and results from these are mostly as expected, with the
   following further details:
   \itemize{
     \item The \code{as.character} function uses "5+" for right censored
     at time 5, "5-" for left censored at time 5, "[2,7]" for an
     observation that was interval censored between 2 and 7,
     "(1,6]" for a counting process data denoting an observation which
     was at risk from time 1 to 6, with an event at time 6, and
     "(1,6+]" for an observation over the same interval but not ending
     with and event.
     For a multi-state survival object the type of event is appended to
     the event time using ":type".

     \item The \code{print} and \code{format} methods make use of
     \code{as.character}.

     \item The \code{as.numeric} and \code{as.integer} methods perform
     these actions on the survival times, but do not affect the
     censoring indicator.

     \item The \code{as.matrix} and \code{t} methods return a matrix

     \item The \code{length} of a \code{Surv} object is the number of
     survival times it contains, not the number of items required to
     encode it, e.g., \code{x <- Surv(1:4, 5:9, c(1,0,1,0)); length(x)}
     has a value of 4.
     Likewise \code{names(x)} will be NULL or a vector of length 4.
     (For technical reasons, any names are actually stored in the
     \code{rownames} attribute of the object.)

     \item For a multi-state survival object \code{levels} returns the
     names of the endpoints, otherwise it is NULL.

     \item The \code{median}, \code{quantile} and \code{plot} methods
     first construct a survival curve using \code{survfit}, then apply
     the appropriate method to that curve.

    \item The concatonation method \code{c()} is asymmetric, its first
     argument determines the exection path.  For instance
     \code{c(Surv(1:4), Surv(5:6))} will concatonate the two objects,
     \code{c(Surv(1:4), 5:6)} will give an error, and
     \code{c(5:6, Surv(1:4))} is equivalent to
     \code{c(5:6, as.vector(Surv(1:4)))}. 
   }
 }
 \seealso{ \code{\link{Surv}}}
 \keyword{survival}
 