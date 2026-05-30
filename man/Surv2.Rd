\name{Surv2}
\alias{Surv2}
\title{Create a survival object}
\description{
  Create a survival object from a timeline style data set. This will
  almost always be the response variable in a formula.
}
\usage{
Surv2(time, event, repeated=FALSE)
}
\arguments{
  \item{time}{a timeline variable, such as age, time from enrollment,
    date, etc.}
  \item{event}{the outcome at that time.  This can be a 0/1 variable,
    TRUE/FALSE, or a factor.
    If the latter, the first level of the factor corresponds to
    `no event was observed at this time'.}
  \item{repeated}{if the same level of the outcome repeats, without an
    intervening event of another type, should this be treated as a new event?}
}

\value{
  An object of class \code{Surv2}.  There are methods for \code{print},
  \code{is.na} and subscripting.
}

\details{
  When used in a \code{coxph} or \code{survfit} model,
  Surv2 acts as a trigger to internally convert a timeline style data
  set into counting process style data, which is then acted on by the
  routine.

  The \code{repeated} argument controls how repeated instances of the same event
  code are treated.  If TRUE, they are treated as new events, an example
  where this might be desired is repeated infections in a subject.
  If FALSE, then repeats are not a new
  event.  An example would be a data set where we wanted to use
  diabetes, say, as an endpoint, but this is diagnosis code is repeated
  in the data at each medical visit.

  The third option is repeated="first", in which case any repeated
  instance of any level is treated as censored.  In an data set with
  states of cognitively unimpaired (CU), mild cognitive impairment
  (MCI), dementia, and death, say, we might want an analysis with
  "ever MCI" and "ever dementia" as the states. A subject with observed
  states of (CU, CU, MCI, CU, MCI, MCI, dementia)
  would become (CU, censor, MCI, CU, MCI, censor, dementia)
  with \code{repeated =FALSE} and (CU, censor, MCI, censor, censor,
  censor, death) with \code{repeated = "first"}.
}
\seealso{
\code{\link{Surv2data}}, \code{\link{coxph}},  
\code{\link{survfit}}  
}
\keyword{survival}
