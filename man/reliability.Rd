\name{reliability}
\docType{data}
\alias{reliability}
\alias{capacitor}
\alias{cracks}
\alias{genfan}
\alias{ifluid}
\alias{imotor}
\alias{turbine}
\alias{valveSeat}
\title{Reliability data sets}
\description{
  A set of data for simple reliablility analyses, taken from the book by
Meeker and Escobar.
}
\usage{data(reliability, package="survival")
}
\details{
  \itemize{
    \item \code{capacitor}:
    Data from a factorial experiment on the life of glass capacitors as a
    function of voltage and operating temperature.  There were 8 capacitors
    at each combination of temperature and voltage.
    Testing at each combination was terminated after the fourth failure.
    \itemize{
      \item \code{temperature}: temperature in degrees celcius
      \item \code{voltage}: applied voltage
      \item \code{time}: time to failure
      \item \code{status}: 1=failed, 0=censored
    }

    \item \code{cracks}: Data on the time until the development of cracks
    in a set of 167 identical turbine parts.
    The parts were inspected at 8 selected times.
    \itemize{
      \item day: time of inspection
      \item fail: number of fans found to have cracks, at this inspection
    }

    \item Data set \code{genfan}: Time to failure of 70 diesel engine fans.
    \itemize{
      \item \code{hours}: hours of service
      \item \code{status}: 1=failure, 0=censored
    }

    Data set \code{ifluid}:
    A data frame with two variables describing the time to electrical
    breakdown of an insulating fluid.
    \itemize{
      \item \code{time}: hours to breakdown
      \item \code{voltage}: test voltage in kV
    }

    \item Data set \code{imotor}: Breakdown of motor insulation as a function of
    temperature.
    \itemize{
      \item temp: temperature of the test
      \item time: time to failure or censoring
      \item status: 0=censored, 1=failed
    }
  
    \item Data set \code{turbine}:
    Each of 432 turbine wheels was inspected
    once to determine whether a crack had developed in the wheel or not.
    \itemize{
      \item hours: time of inspection (100s of hours)
      \item inspected: number that were inspected
      \item failed: number that failed
    }

    \item Data set \code{valveSeat}:
    Time to replacement of valve seats for 41 diesel engines.  More than
    one seat may be replaced at a particular service, leading to duplicate
    times in the data set.  The final inspection time for each engine will
    have status=0.
    \itemize{
      \item id: engine identifier
      \item time: time of the inspection, in days
      \item status: 1=replacement occured, 0= not
    }
  }
}
\references{
  Meeker and Escobar, Statistical Methods for Reliability Data, 1998.
}

\examples{
survreg(Surv(time, status) ~ temperature + voltage, capacitor)

# Replacement of valve seats.  In this case the cumulative hazard is the 
#  natural target, an estimate of the number of replacements by a given time
#  (known as the cumulative mean function = CMF in relability).
# When two valve seats failed at the same inspection, we need to jitter one
#  of the times, to avoid a (time1, time2) interval of length 0
ties <- which(with(valveSeat, diff(id)==0 & diff(time)==0))  #first of a tie
temp <- valveSeat$time
temp[ties] <- temp[ties] - .1 # jittered time
vdata <- valveSeat
vdata$time1 <- ifelse(!duplicated(vdata$id), 0, c(0, temp[-length(temp)]))
vdata$time2 <- temp
fit2 <- survfit(Surv(time1, time2, status) ~1, vdata, id=id)
\dontrun{
plot(fit2, cumhaz= TRUE, xscale= 365.25, 
      xlab="Years in service", ylab = "Expected number of repairs")
}
}
\keyword{datasets}
