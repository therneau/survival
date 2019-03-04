\name{transplant}
\alias{transplant}
\alias{transplant2}
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
\usage{data("transplant2")}
\format{
  A data frame with 861 observations on the following 6 variables.
  \describe{
    \item{\code{age}}{age at addition to the waiting list}
    \item{\code{sex}}{\code{m} or \code{f}}
    \item{\code{abo}}{blood type: \code{A}, \code{B}, \code{AB}  or \code{O}}
    \item{\code{year}}{year in which they entered the waiting list}
    \item{\code{futime}}{time from entry to final disposition}
    \item{\code{event}}{final disposition: \code{censored},
      \code{death},
      \code{ltx} or \code{withdraw}}
    \item{\code{creat}}{serum creatinine}
    \item{\code{bili}}{serum bilirubin}
    \item{\code{inr}}{International Normalized Ratio, a measure of the
      blood's clotting ability}
    \item{\code meld}{calculated MELD score}
    \item{diag}{primary diagnosis: alcoholic liver disease,
      cholestatic liver disease, hepatitis B, hepatitis C, or other}
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
types, whereas an AB liver can only be used by an AB recipient.
Thus type O subjects on the waiting list are at a disadvantage, since
the pool of competitors is larger for type O donor livers.

This data is of historical interest and provides a useful example of
competing risks, but it has little relevance to current
practice.  Liver allocation policies have evolved and now depend
directly on each individual patient's risk and need, assessments of which are
regularly updated while a patient is on the waiting list.
The overall organ shortage remains acute, however.

The \code{transplant} data set was a version early in the analysis,
\code{transplant2} was the final data set and matches the paper.
}
\examples{
#since event is a factor, survfit creates competing risk curves
pfit <- survfit(Surv(futime, event) ~ abo, transplant)
pfit[,2]  #time to liver transplant, by blood type
plot(pfit[,2], mark.time=FALSE, col=1:4, lwd=2, xmax=735,
       xscale=30.5, xlab="Months", ylab="Fraction transplanted",
       xaxt = 'n')
temp <- c(0, 6, 12, 18, 24)
axis(1, temp*30.5, temp)
legend(450, .35, levels(transplant$abo), lty=1, col=1:4, lwd=2, bty='n')

# competing risks for type O
plot(pfit[4,], xscale=30.5, xmax=735, col=1:3, lwd=2)
legend(450, .4, c("Death", "Transpant", "Withdrawal"), col=1:3, lwd=2)
}
\references{
 Kim WR, Therneau TM, Benson JT, Kremers WK, Rosen CB, Gores GJ, Dickson
 ER. 
 Deaths on the liver transplant waiting list: An analysis of competing risks. 
 Hepatology 2006 Feb; 43(2):345-51.
 }
\keyword{datasets}
