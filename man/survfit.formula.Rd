\name{survfit.formula}
\alias{survfit.formula}
\alias{[.survfit}
\title{
  Compute a Survival Curve for Censored Data 
}
\description{
Computes an estimate of a survival curve for censored data 
using either the Kaplan-Meier or the Fleming-Harrington method. 
For competing risks data it computes the
cumulative incidence curve. 
}
\usage{
\method{survfit}{formula}(formula, data, weights, subset, na.action,  
        etype, id, istate, ...)
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
  \item{etype}{
    a variable giving the type of event.  This has been superseded by
    multi-state Surv objects; see example below.
  }
  \item{id}{
    identifies individual subjects, when a given person can have multiple
    lines of data.
  }
  \item{istate}{for multi-state models, identifies the initial state of
  each subject} 

  \item{\dots}{
    The following additional arguments are passed to internal functions
    called by \code{survfit}.
    \describe{
      \item{type}{
	a character string specifying the type of survival curve. 
	Possible values are \code{"kaplan-meier"}, 
	\code{"fleming-harrington"} or \code{"fh2"} 
	if a formula is given.  This is ignored for competing risks or when
	the Turnbull estimator is used.
      }
      \item{error}{
	a character string specifying the error.  Possible values are 
	\code{"greenwood"} for the Greenwood formula or 
	\code{"tsiatis"} for the Tsiatis formula, 
	(only the first character is 
	necessary).  
      }
      \item{conf.type}{
	One of \code{"none"}, \code{"plain"}, \code{"log"} (the default),
	or \code{"log-log"}.  Only
	enough of the string to uniquely identify it is necessary.
	The first option causes confidence intervals not to be
	generated.  The second causes the standard intervals
	\code{curve +- k *se(curve)}, where k is determined from
	\code{conf.int}.  The log option calculates intervals based on the
	cumulative hazard or log(survival). The last option bases
	intervals on the log hazard or log(-log(survival)). 
      }
      \item{conf.lower}{
	a character string to specify modified lower limits to the curve, the 
	upper limit remains unchanged.  
	Possible values are \code{"usual"} (unmodified), 
	\code{"peto"}, 
	and \code{"modified"}.  T
	he modified lower limit 
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
The estimates used are the Kalbfleisch-Prentice 
(Kalbfleisch and Prentice, 1980, p.86) and the Tsiatis/Link/Breslow, 
which reduce to the Kaplan-Meier and Fleming-Harrington estimates, 
respectively, when the weights are unity.  

The Greenwood formula for the variance is a sum of terms 
d/(n*(n-m)), where d is the number of deaths at a given time point, n 
is the sum of weights for all individuals still at risk at that time, and 
m is the sum of weights for the deaths at that time.  The 
justification is based on a binomial argument when weights are all 
equal to one; extension to the weighted case is ad hoc.  Tsiatis 
(1981) proposes a sum of terms d/(n*n), based on a counting process 
argument which includes the weighted case. 

The two variants of the F-H estimate have to do with how ties are handled.
If there were 3 deaths out of 10 at risk, then the first
increments the hazard by 3/10 and the second
by 1/10 + 1/9 + 1/8.
For the first method S(t) = exp(H), where H is 
the Nelson-Aalen cumulative hazard estimate,
whereas the \code{fh2} method will 
give results S(t) results closer to the Kaplan-Meier. 

When the data set includes left censored or interval censored data (or both),
then the EM approach of Turnbull is used to compute the overall curve.
When the baseline method is the Kaplan-Meier, this is known to converge to
the maximum likelihood estimate.

The cumulative incidence curve is an alternative to the Kaplan-Meier
for competing risks data.
For instance, in patients with MGUS, conversion to an overt plasma
cell malignancy occurs at a nearly constant rate among those still alive.    
A Kaplan-Meier estimate, treating death due to other causes as censored,
gives a 20 year cumulate rate of 33\% for the 241 early patients of Kyle.  
This estimates the incidence of conversion if all other causes of death were
removed, which is an unrealistic assumption given the mean starting age
of 63 and a median follow up of over 21 years.
  
The CI estimate, on the other hand, estimates the total number of
conversions that will actually occur.  Because the population is older,
this is much smaller than the KM, 22\% at 20 years for Kyle's data.
If there were no censoring, then CI(t) could very simply be computed as
total number
of patients with progression by time t divided by the sample size n.

}
\section{References}{
Dorey, F. J. and Korn, E. L. (1987).  Effective sample sizes for confidence 
intervals for survival probabilities.  \emph{Statistics in Medicine} 
\bold{6}, 679-87. 

Fleming, T. H. and Harrington, D. P. (1984).  Nonparametric estimation of the 
survival distribution in censored data.  \emph{Comm. in Statistics}  
\bold{13}, 2469-86. 

Kablfleisch, J. D. and Prentice, R. L. (1980).  
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
\code{\link{survfit.coxph}} for survival curves from Cox models.
  
\code{\link{print}},  
\code{\link{plot}},  
\code{\link{lines}},   
\code{\link{coxph}},  
\code{\link{Surv}},  
\code{\link{strata}}.   
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
fitKM <- survfit(Surv(stop, event=='progression') ~1, data=mgus1,
                    subset=(start==0))

fitCI <- survfit(Surv(stop, status*as.numeric(event), type="mstate") ~1,
                    data=mgus1, subset=(start==0))

# CI curves are always plotted from 0 upwards, rather than 1 down
plot(fitCI, xscale=365.25, xmax=7300, mark.time=FALSE,
            col=2:3, xlab="Years post diagnosis of MGUS")
lines(fitKM, fun='event', xscale=365.25, xmax=7300, mark.time=FALSE,
            conf.int=FALSE)
text(10, .4, "Competing risk: death", col=3)
text(16, .15,"Competing risk: progression", col=2)
text(15, .30,"KM:prog")
}
\keyword{survival}
