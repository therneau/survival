\name{mgus2}
\alias{mgus2}
\docType{data}
\title{Monoclonal gammapothy data}

\description{Natural history of 1341 sequential patients with monoclonal
  gammapothy of undetermined significance (MGUS).
}
\usage{data("mgus2")}
\format{
  A data frame with 1384 observations on the following 10 variables.
  \describe{
    \item{\code{id}}{subject identifier}
    \item{\code{age}}{age at diagnosis, in years}
    \item{\code{sex}}{a factor with levels \code{F} \code{M}}
    \item{\code{hgb}}{hemoglobin}
    \item{\code{creat}}{creatinine}
    \item{\code{mspike}}{size of the monoclonal serum splike}
    \item{\code{ptime}}{time until progression to a plasma cell
      malignancy (PCM) or last contact, in months}
    \item{\code{pstat}}{occurence of PCM: 0=no, 1=yes }
    \item{\code{futime}}{time until death or last contact, in months}
    \item{\code{death}}{occurence of death: 0=no, 1=yes}
  }
}
\details{
  This is a larger follow-on study of the condition also found in data
  set \code{mgus}.
}
\source{Mayo Clinic data courtesy of Dr. Robert Kyle.  All patient
  identifiers have been removed, age rounded to the nearest year, and
  follow-up times rounded to the nearest month.}
\references{
  R. Kyle, T. Therneau, V. Rajkumar, J. Offord, D. Larson, M. Plevak,
  and L. J. Melton III, A long-terms study of prognosis in monoclonal
  gammopathy of undertermined significance. New Engl J Med, 346:564-569 (2002).
}
\keyword{datasets}
