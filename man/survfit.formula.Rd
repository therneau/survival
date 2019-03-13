\name{survfit.formula}
\alias{survfit.formula}
\alias{[.survfit}
\title{
  Compute a Survival Curve for Censored Data 
}
\description{
  Computes an estimate of a survival curve for censored data using
  the Aalen-Johansen estimator.  For ordinary (single event) survival
  this reduces to the Kaplan-Meier estimate.
}
\usage{
\method{survfit}{formula}(formula, data, weights, subset, na.action,  
        type=1, id, cluster, istate, timefix=TRUE,
        etype, error,  ...)
}
\arguments{
  \item{formula}{
    a formula object, which must have a 
    \code{Surv} object as the  
    response on the left of the \code{~} operator and, if desired, terms  
    separated by + operators on the right. 
    One of the terms may be a \code{strata} object.
    For a single survival curve the right hand side should be \code{~ 1}.
  }
  \item{data}{
    a data frame in which to interpret the variables named in the formula, 
    \code{subset} and \code{weights} arguments. 
  }
  \item{weights}{
    The weights must be nonnegative and it is strongly recommended that  
    they be strictly positive, since zero weights are ambiguous, compared 
    to use of the \code{subset} argument.
  }
  \item{subset}{
    expression saying that only a subset of the rows of the data 
    should be used in the fit. 
  }
  \item{na.action}{
    a missing-data filter function, applied to the model frame, after any 
    \code{subset} argument has been used. 
    Default is \code{options()$na.action}. 
  }
  \item{type}{the type of estimates to produce; see details below}

  \item{id}{
    identifies individual subjects, used when a given person can have multiple
    lines of data or to signal a robust error estimate
  }

  \item{id}{
    identifies individual subjects, when a given person can have multiple
    lines of data.
  }
  \item{cluster}{used to group observations for the infinitesmal
    jackknife variance estimate, defaults to the value of id}
  
  \item{istate}{for multi-state models, identifies the initial state of
    each subject or observation}

  \item{timefix}{process times through the \code{aeqSurv} function to
  eliminate potential roundoff issues.}

  \item{etype}{
    a variable giving the type of event.  This has been superseded by
    multi-state Surv objects and is depricated; see example below.
  }
  \item{error}{this argument is no longer used}
  \item{\dots}{
    The following additional arguments are passed to internal functions
    called by \code{survfit}.
    \describe{
      \item{se.fit}{logical value, default is TRUE.  If FALSE then
	standard error computations are omitted.
	}
      \item{conf.type}{
	One of \code{"none"}, \code{"plain"}, \code{"log"} (the default),
	\code{"log-log"} or \code{"logit"}.  Only
	enough of the string to uniquely identify it is necessary.
	The first option causes confidence intervals not to be
	generated.  The second causes the standard intervals
	\code{curve +- k *se(curve)}, where k is determined from
	\code{conf.int}.  The log option calculates intervals based on the
	cumulative hazard or log(survival). The log-log option bases the
	intervals on the log hazard or log(-log(survival)), and the
	logit option on log(survival/(1-survival)). 
      }
      \item{conf.lower}{
	a character string to specify modified lower limits to the curve, the 
	upper limit remains unchanged.  
	Possible values are \code{"usual"} (unmodified), 
	\code{"peto"}, 
	and \code{"modified"}.  The modified lower limit 
	is based on an "effective n" argument.  The confidence 
	bands will agree with the usual calculation at each death time, but unlike 
	the usual bands the confidence interval becomes wider at each censored 
	observation.  The extra width is obtained by multiplying the usual 
	variance by a factor m/n, where n is the number currently at risk and 
	m is the number at risk at the last death time.  (The bands thus agree 
	with the un-modified bands at each death time.) 
	This is especially useful for survival curves with a long flat tail. 
	
	The Peto lower limit is based on the same "effective n" argument as the 
	modified limit, but also replaces the usual Greenwood variance term with 
	a simple approximation.  It is known to be conservative. 
      }
      \item{start.time}{
	numeric value specifying a time to start calculating survival
	information.
	The resulting curve is the survival conditional on surviving to
	\code{start.time}.
      }
      \item{conf.int}{
	the level for a two-sided confidence interval on the survival curve(s). 
	Default is 0.95. 
      }
      \item{se.fit}{
	a logical value indicating whether standard errors should be 
	computed.  Default is \code{TRUE}.
      }
      \item{influence}{a logical value indicating whether to return the
	infinitesimal jackknife (influence) values for each subject.
	These contain the values of the derivative of each value
	with respect to the case weights of each subject i:
	\eqn{\partial p/\partial w_i}{dp/dw[i]},
	evaluated at the vector of weights.
      }

      \item{p0}{this applies only to multi-state curves. 
	An optional vector giving the initial probability across
	the states. If this is missing, then p0 is estimated using the
	frequency of the starting states of all observations at risk
	at \code{start.time}, or if that is not specified, at
	the time of the first event.}
    }
  }
}
\value{
an object of class \code{"survfit"}.  
See \code{survfit.object} for 
details. Methods defined for survfit objects are  
\code{print}, \code{plot}, 
\code{lines}, and \code{points}. 
}
\details{
  If there is a \code{data} argument, then variables in the \code{formula},
  code{weights}, \code{subset}, \code{id}, \code{cluster} and
  \code{istate} arguments will be searched for in that data set.

  The routine returns both an estimated probability in state and an
  estimated cumulative hazard estimate.
  The cumulative hazard estimate is the Nelson-Aalen (NA) estimate or the
  Fleming-Harrington (FH) estimate, and the probability in state will be
  either the Aalen-Johansen (AJ) estimate or the exponential of the
  cumulative hazard.
  For \code{type=1} gives the AJ + NA, \code{type=2} the AJ + FH,
  \code{type=3} exp(NA) + NA and \code{type=4} exp(FH) + FH.
  For single state data the AJ estimate reduces to the Kaplan-Meier and
  the probability in state to the survival curve; 
  for competing risks data the AJ reduces to the cumulative incidence (CI)
  estimator.
  
  For backward compatability the \code{type} argument can also be a
  character string of "kaplan-meier", "fleming-harrington" or "fh2",
  which map to type = 1,3, and 4, respectively.  Any unique abbreviation
  of the string is also allowed.

  If a \code{id} or \code{cluster} argument is present, or for multi-state
  curves, then the standard
  errors of the results will be based on an infinitesimal jackknife (IJ)
  estimate, otherwise the standard model based estimate will be used.
  With the IJ estimate, the leverage values themselves can be returned
  as arrays with dimensions: number of subjects, number of unique times,
  and for a multi-state model, the number of unique states.
  Be forwarned that these arrays can be huge.  If there is a
  \code{cluster} argument this first dimension will be the number of
  clusters and the variance will be a grouped IJ estimate; this can be
  an important tool for reducing the size.
  A numeric value for the \code{influence} argument allows finer
  control: 0= return neither (same as FALSE), 1= return the influence
  array for probability in state, 2= return the influence array for the
  cumulative hazard, 3= both (same as TRUE).

  When the data set includes left censored or interval censored data (or both),
  then the EM approach of Turnbull is used to compute the overall curve.
  Currently this algorithm is very slow, only a survival curve is
  produced, and it does not support the IJ variance.
  }
\section{References}{
Dorey, F. J. and Korn, E. L. (1987).  Effective sample sizes for confidence 
intervals for survival probabilities.  \emph{Statistics in Medicine} 
\bold{6}, 679-87. 

Fleming, T. H. and Harrington, D. P. (1984).  Nonparametric estimation of the 
survival distribution in censored data.  \emph{Comm. in Statistics}  
\bold{13}, 2469-86. 

Kalbfleisch, J. D. and Prentice, R. L. (1980).
\emph{The Statistical Analysis of Failure Time Data.}
New York:Wiley. 

Kyle, R. A. (1997).
Moncolonal gammopathy of undetermined significance and solitary
      plasmacytoma. Implications for progression to overt multiple myeloma\},
\emph{Hematology/Oncology Clinics N. Amer.}
\bold{11}, 71-87.

Link, C. L. (1984). Confidence intervals for the survival 
function using Cox's proportional hazards model with  
covariates.  \emph{Biometrics}  
\bold{40}, 601-610.

Turnbull, B. W. (1974).  Nonparametric estimation of a survivorship
function with doubly censored data. \emph{J Am Stat Assoc},
\bold{69}, 169-173. 

}
\seealso{
  \code{\link{survfit.coxph}} for survival curves from Cox models,
  \code{\link{survfit.object}} for a description of the components of a
  survfit object,
\code{\link{print.survfit}},  
\code{\link{plot.survfit}},  
\code{\link{lines.survfit}},   
\code{\link{coxph}},  
\code{\link{Surv}}.  
}
\examples{
#fit a Kaplan-Meier and plot it 
fit <- survfit(Surv(time, status) ~ x, data = aml) 
plot(fit, lty = 2:3) 
legend(100, .8, c("Maintained", "Nonmaintained"), lty = 2:3) 

#fit a Cox proportional hazards model and plot the  
#predicted survival for a 60 year old 
fit <- coxph(Surv(futime, fustat) ~ age, data = ovarian) 
plot(survfit(fit, newdata=data.frame(age=60)),
     xscale=365.25, xlab = "Years", ylab="Survival") 

# Here is the data set from Turnbull
#  There are no interval censored subjects, only left-censored (status=3),
#  right-censored (status 0) and observed events (status 1)
#
#                             Time
#                         1    2   3   4
# Type of observation
#           death        12    6   2   3
#          losses         3    2   0   3
#      late entry         2    4   2   5
#
tdata <- data.frame(time  =c(1,1,1,2,2,2,3,3,3,4,4,4),
                    status=rep(c(1,0,2),4),
                    n     =c(12,3,2,6,2,4,2,0,2,3,3,5))
fit  <- survfit(Surv(time, time, status, type='interval') ~1, 
              data=tdata, weight=n)

#
# Time to progression/death for patients with monoclonal gammopathy
#  Competing risk curves (cumulative incidence)
#
fitKM <- survfit(Surv(stop, event=='pcm') ~1, data=mgus1,
                    subset=(start==0))

fitCI <- survfit(Surv(stop, event) ~1,
                    data=mgus1, subset=(start==0))
\dontrun{
# CI curves show the probability in state
plot(fitCI, xscale=365.25, xmax=7300, mark.time=FALSE,
            col=2:3, xlab="Years post diagnosis of MGUS",
            ylab="P(state)")
lines(fitKM, fun='event', xmax=7300, mark.time=FALSE,
            conf.int=FALSE)
text(3652, .4, "Competing risk: death", col=3)
text(5840, .15,"Competing risk: progression", col=2)
text(5480, .30,"KM:prog")
}
}
\keyword{survival}
