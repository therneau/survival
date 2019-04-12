\name{mgus2}
\alias{mgus2}
\docType{data}
\title{Monoclonal gammopathy data}

\description{Natural history of 1341 sequential patients with monoclonal
  gammopathy of undetermined significance (MGUS).  This is a superset of
  the \code{mgus} data, at a later point in the accrual process
}
\usage{data("mgus2")}
\format{
  A data frame with 1384 observations on the following 10 variables.
  \describe{
    \item{\code{id}}{subject identifier}
    \item{\code{age}}{age at diagnosis, in years}
    \item{\code{sex}}{a factor with levels \code{F} \code{M}}
    \item{\code{dxyr}}{year of diagnosis}
    \item{\code{hgb}}{hemoglobin}
    \item{\code{creat}}{creatinine}
    \item{\code{mspike}}{size of the monoclonal serum splike}
    \item{\code{ptime}}{time until progression to a plasma cell
      malignancy (PCM) or last contact, in months}
    \item{\code{pstat}}{occurrence of PCM: 0=no, 1=yes }
    \item{\code{futime}}{time until death or last contact, in months}
    \item{\code{death}}{occurrence of death: 0=no, 1=yes}
  }
}
\details{
  This is an extension of the study found in the \code{mgus} data set,
  containing enrollment through 1994 and follow-up through 1999.
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
