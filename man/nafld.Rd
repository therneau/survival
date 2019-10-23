\name{nafld}
\alias{nafld1}
\alias{nafld2}
\alias{nafld3}
\docType{data}
\title{Non-alcohol fatty liver disease}
\description{
  Data sets containing the data from a population study of non-alcoholic
  fatty liver disease (NAFLD).  Subjects with the condition and a set of
  matched control subjects were followed forward for metabolic
  conditions, cardiac endpoints, and death.
}
\usage{nafld1
       nafld2
       nafld3}
\format{
  \code{nafld1} is a data frame with 17549 observations on the following 10 variables.
  \describe{
    \item{\code{id}}{subject identifier}
    \item{\code{age}}{age at entry to the study}
    \item{\code{male}}{0=female, 1=male}
    \item{\code{weight}}{weight in kg}
    \item{\code{height}}{height in cm}
    \item{\code{bmi}}{body mass index}
    \item{\code{case.id}}{the id of the NAFLD case to whom this subject
      is matched}
    \item{\code{futime}}{time to death or last follow-up}
    \item{\code{status}}{0= alive at last follow-up, 1=dead}
  }
  
  \code{nafld2} is a data frame with 400123 observations and 4 variables
  containing laboratory data
  \describe{
    \item{\code{id}}{subject identifier}
    \item{\code{days}}{days since index date}
    \item{\code{test}}{the type of value recorded}
    \item{\code{value}}{the numeric value}
  }
  
  \code{nafld3} is a data frame with 34340 observations and 3 variables
  containing outcomes
  \describe{
    \item{\code{id}}{subject identifier}
    \item{\code{days}}{days since index date}
    \item{\code{event}}{the endpoint that occurred}
  }
}
\details{
The primary reference for the NAFLD study is Allen (2018).
The incidence of non-alcoholic fatty liver disease (NAFLD) has been rising
rapidly in the last decade and it is now one of the main drivers of
hepatology practice \cite{Tapper2018}.
It is essentially the presence of excess fat in the
liver, and parallels the ongoing obesity epidemic.
Approximately 20-25\% of NAFLD patients will develop the inflammatory state
of non-alcoholic steatohepatitis (NASH), leading to fibrosis 
and eventual end-stage liver disease.
NAFLD can be accurately diagnosed by MRI methods, 
but NASH diagnosis currently requires a biopsy.

The current study constructed a population cohort of
all adult NAFLD subjects from 1997 to 2014  along with 4 potential
controls for each case.
To protect patient confidentiality all time intervals are in days since
the index date; none of the dates from the original data were retained.
Subject age is their integer age at the index date, and the subject
identifier is an arbitrary integer. 
As a final protection, we include only a 90\% random sample of the data.
As a consequence analyses results will not exactly match the
original paper.

There are 3 data sets: \code{nafld1} contains baseline data and has
one observation per subject, \code{nafld2} has one observation for
each (time dependent) continuous measurement,
and \code{nafld3} has one observation for
each yes/no outcome that occured.
}

\source{Data obtained from the author.}
\references{
AM Allen, TM Therneau, JJ Larson, A Coward, VK Somers and PS Kamath,
   Nonalcoholic Fatty Liver Disease Incidence and Impact on Metabolic
            Burden and Death: A 20 Year Community Study,
   Hepatology 67:1726-1736, 2018.
}

\keyword{datasets}
