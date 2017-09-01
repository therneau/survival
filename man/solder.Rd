\name{solder}
\alias{solder}
\docType{data}
\title{Data from a soldering experiment}
\description{
In 1988 an experiment was designed and implemented at one of AT&T's
factories to investigate alternatives in the "wave soldering" procedure
for mounting electronic componentes to printed circuit boards.
The experiment varied a number of factors relevant to the process.
The response, measured by eye, is the number of visible solder skips.
}
\usage{data("solder")}
\format{
  A data frame with 900 observations on the following 6 variables.
  \describe{
    \item{\code{Opening}}{the amount of clearance around the mounting
      pad (3 levels)}
    \item{\code{Solder}}{the amount of solder (2 levels)}
    \item{\code{Mask}}{type and thickness of the material used for the
      solder mask (5 levels)}
    \item{\code{PadType}}{the geometry and size of the mounting pad (10 levels)}
    \item{\code{Panel}}{each board was divided into 3 panels}
    \item{\code{skips}}{the number of skips}
  }
}
\details{
This data set is used as a detailed example in chapter 1 of Chambers and
Hastie, though they largely focus on a balanced subset of the data.
The data set was included in Splus.
In the survival package the data is used as part of the tests for the
\code{yates} command.
}
\references{
J Chambers and T Hastie, Statistical models in S.  Chapman and Hall, 1993.
}
\examples{
data(solder)
# The balanced subset contains the first 180 of each mask, and deletes
#  mask A6. 
index <- 1 + (1:nrow(solder)) - match(solder$Mask, solder$Mask)
solder.balance <- droplevels(subset(solder, Mask != "A6" & index <= 180))
}
\keyword{datasets}
