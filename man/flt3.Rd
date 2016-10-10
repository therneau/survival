\name{flt3}
\alias{flt3}
\docType{data}
\title{FLT3-mutated acute myeloid leukemia trial}
\description{
This simulated data set is based on Cancer and Luekemia Group B (CALGB)
trial 10603, midostaurin plus chemotherapy in FLT3-mutated acute
myeloid leukemia.
Information about the actual study is found in the reference.
}
\usage{data("flt3")}
\format{
  A data frame with 646 observations on the following 9 variables.
  \describe{
    \item{\code{id}}{subject identifier, 1-646}
    \item{\code{trt}}{treatment arm A or B}
    \item{\code{sex}}{a factor with levels of \code{f} \code{m}}
    \item{\code{group}}{stratification factor for the randomization}
    \item{\code{futime}}{time to death or last follow-up}
    \item{\code{death}}{1 if \code{futime} is a death, 0 for censoring}
    \item{\code{txtime}}{time to hematropetic stem cell transplant}
    \item{\code{crtime}}{time to complete response}
    \item{\code{rltime}}{time to relapse of disease}
  }
}
\details{
  Mutation in the gene encoding the trans membrane tyrosine kinase FLT3
  occurs in about 30 percent of adult newly diagnosed with 
  acute myeloid leukemia (AML); these patients have a poor outcome.
  This was a trial to see if an FLT inhibitor would be a beneficial
  addition to standard treatment.

  The subjects are stratified by the type of mutation: an internal
  tandem duplication (ITD), thought to confer poor prognosis, and tyrosine
  kinase point mutations, which are of uncertain effect.  To be eligible
  for the study at least one of these needed to be present.  The groups
  are TKD without ITD, ITD with an alleic ratio of less than 0.7 (with
  or without TKD),
  and ITD with an alleic ratio greater than or equal to 0.7 (with or
  without TKD).
}
\references{
  Richard M Stone, Sumithra J  Mandrekar, et al (29 authors).
  Midostaurin plus chemotharapy in FLT3-mutated acute myeloid leukemia.
  submitted, 2016.
  }
\examples{
coxph(Surv(futime, death) ~ trt, data=flt3)
# See the competing risk vignette for a more complete analysis
}
\keyword{datasets}
