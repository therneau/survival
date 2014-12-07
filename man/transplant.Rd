\name{transplant}
\alias{transplant}
\docType{data}
\title{Liver transplant waiting list}
\description{
  Subjects on a liver transplant waiting list from 1990-1999, and their
  disposition: received a transplant, died while waiting, withdrew from
  the list, or censored. 
}
\usage{data("transplant")}
\format{
  A data frame with 815 observations on the following 6 variables.
  \describe{
    \item{\code{age}}{age at addition to the waiting list}
    \item{\code{sex}}{\code{m} or \code{f}}
    \item{\code{abo}}{blood type: \code{A}, \code{B}, \code{AB}  or \code{O}}
    \item{\code{year}}{year in which they entered the waiting list}
    \item{\code{futime}}{time from entry to final disposition}
    \item{\code{event}}{final disposition: \code{censored},
      \code{death},
      \code{ltx} or \code{withdraw}}
  }
}
\details{
This represents the transplant experience in a particular region,
over a time period in which liver transplant became much more widely
recognized as a viable treatment modality.
The number of liver transplants rises over the period, but the number of
subjects added to the liver transplant waiting list grew much faster.
Important questions addressed by the data are the change in waiting
time, who waits, and whether there was an consequent increase in deaths
while on the list.

Blood type is an important consideration.  Donor livers from subjects
with blood type O can be used by patients with A, B, AB or 0 blood
types, whereas a recipient of type B cannot accept an A or AB liver
for instance.
Thus type O subjects on the waiting list are at a disadvantage, since
the pool of competitors is larger for type O donor livers.

This data is of historical interest but has little relevance to current
practice.  Liver allocation policies have evolved and now depend
directly on individual patient's risk and need, assessments of which are
regularly updated while a patient is on the waiting list.
The overall organ shortage remains acute, however.
}
\examples{
period <- cut(transplant$year, c(1989, 1992, 1995, 1997, 2000),
              labels=c('90-92', '93-95', '96-97', '98-99'))
pfit <- survfit(Surv(futime, event) ~ period, transplant)
pfit[2,]  #time to liver transplant
plot(pfit[,2], mark.time=FALSE, col=1:4, lwd=2, xmax=735,
       xscale=30.5, xlab="Months", ylab="Fraction transplanted",
       xaxt = 'n')
temp <- c(0, 6, 12, 18, 24)
axis(1, temp, temp)

legend(15, .35, levels(period), lty=1, col=1:4, lwd=2, bty='n')
}
\keyword{datasets}
