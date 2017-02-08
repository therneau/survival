\name{summary.pyears}
\alias{summary.pyears}
\title{Summary function for pyears objecs}
\description{Create a printable table of a person-years result.}

\usage{
\method{summary}{pyears}(object, header = TRUE, call = header, n = TRUE,
event = TRUE, pyears = TRUE, expected = TRUE, rate = FALSE, rr =expected,
ci.r = FALSE, ci.rr = FALSE, totals=FALSE, legend = TRUE, vline = FALSE,
vertical= TRUE, nastring=".", conf.level = 0.95,
scale = 1, ...)
}

\arguments{
  \item{object}{a pyears object}
  \item{header}{print out a header giving the total number of
    observations, events, person-years, and total time (if any)
    omitted from the table}
  \item{call}{print out a copy of the call}
  \item{n, event, pyears, expected}{logical arguments: should these
    elements be printed in the table?}
  \item{rate, ci.r}{logical arguments: should the incidence rate and/or
    its confidence interval be given in the table?}
  \item{rr, ci.rr}{logical arguments: should the hazard ratio and/or
    its confidence interval be given in the table?}
  \item{totals}{should row and column totals be added?}
  \item{legend}{should a legend be included in the printout?}
  \item{vline}{should vertical lines be included in the printed tables?}
  \item{vertical}{when there is only a single predictor, should the
    table be printed with the predictor on the left (vertical=TRUE) or
    across the top (vertical=FALSE)?}
  \item{nastring}{what to use for missing values in the table.  Some of
  these are structural, e.g., risk ratios for a cell with no follow-up time.}
  \item{conf.level}{confidence level for any confidence intervals}
  \item{scale}{a scaling factor for printed rates}
  \item{\dots}{optional arguments which will be passed to the
    \code{format} function; common choices would be digits=2 or nsmall=1.}
  
}

\details{
  The \code{pyears} function is often used to create initial
  descriptions of a survival or time-to-event variable; the type of
  material that is often found in ``table 1'' of a paper.  The summary
  routine prints this information out using one of pandoc table styles.
  A primary reason for choosing this style is that Rstudio is then able
  to automatically render the results in multiple formats: html, rtf,
  latex, etc.

  If the \code{pyears} call has only a single covariate then the table
  will have that covariate as one margin and the statistics of interest
  as the other.
  If the \code{pyears} call has two predictors then those two predictors
  are used as margins of the table, while each cell of the table
  contains the statistics of interest as multiple rows within the cell.
  If there are more than two predictors then multiple tables are
  produced, in the same order as the standard R printout for an array.

  The "N" entry of a pyears object is the number of observations which
  contributed to a particular cell.  When the original call includes
  \code{tcut} objects then a single observation may contribute to
  multiple cells.
}
\section{Notes}{
  The pandoc system has four table types: with or without vertical bars,
  and with single or multiple rows of data in each cell.
  This routine produces all 4 styles depending on options, but currently
  not all of them are recognized by the Rstudio-pandoc pipeline.
  (And we don't yet see why.)
  }
\value{a copy of the object}
\author{Terry Therneau and Elizabeth Atkinson}
\seealso{\code{\link{cipoisson}}, \code{\link{pyears}}, \code{\link{format}}}
\keyword{ survival }

