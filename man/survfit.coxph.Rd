\name{survfit.coxph}
\alias{survfit.coxph}
\title{
  Compute a Survival Curve from a Cox model
}
\description{
Computes the predicted survivor function for a Cox proportional 
hazards model. 
}
\usage{
\method{survfit}{coxph}(formula, newdata, 
        se.fit=TRUE, conf.int=.95,
        individual=FALSE, 
        type,vartype,
        conf.type=c("log","log-log","plain","none"), censor=TRUE, id,
        newstrata, na.action=na.pass, ...)
}
\arguments{
  \item{formula}{
    A \code{coxph} object. 
  }
  \item{newdata}{
    a data frame with the same variable names as those that appear 
    in the \code{coxph} formula.  It is also valid to use a vector, if
    the data frame would consist of a single row.
    
    The curve(s) produced will be representative of a cohort whose 
    covariates correspond to the values in \code{newdata}. 
    Default is the mean of the covariates used in the 
    \code{coxph} fit. 
  }
  \item{individual}{
    This argument has been superseded by the \code{id} argument and is
    present only for backwards compatability.
    A logical value indicating whether each row of \code{newdata} represents
    a distinct individual (FALSE, the default), 
    or if each row of the data frame represents different 
    time epochs for only one individual (TRUE).  
    In the former case the result will have one curve for each row in
    \code{newdata}, in the latter only a single curve will be produced.
  }
  \item{conf.int}{
    the level for a two-sided confidence interval on the survival curve(s). 
    Default is 0.95. 
  }
  \item{se.fit}{
    a logical value indicating whether standard errors should be 
    computed.  Default is \code{TRUE}. 
  }
  \item{type,vartype}{
    a character string specifying the type of survival curve. 
    Possible values are  \code{"aalen"}, \code{"efron"}, or
    \code{"kalbfleish-prentice"} 
    (only the first two characters are necessary). 
    The default is to match the computation used in the Cox model.
    The Nelson-Aalen-Breslow estimate for \code{ties='breslow'}, the Efron
    estimate for \code{ties='efron'} and the Kalbfleisch-Prentice
    estimate for a discrete time model \code{ties='exact'}.
    Variance estimates are the Aalen-Link-Tsiatis, Efron, and
    Greenwood.
    The default will be the Efron
    estimate for \code{ties='efron'} and the Aalen estimate otherwise.    
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
  \item{censor}{if FALSE time points at which there are no events (only
    censoring) are not included in the result.}
  \item{id}{optional variable name of subject identifiers.  If this is
    present, then each group of rows with the same subject id represents
    the covariate path through time of a single subject, and the result
    will contain one curve per subject.  If the \code{coxph} fit had
    strata then that must also be specified in \code{newdata}.
    If missing, then each
    individual row of \code{newdata} is presumed to represent a distinct
    subject and there will be \code{nrow(newdata)} times the number of
    strata curves in the result (one for each strata/subject combination).
    result.}

  \item{newstrata}{if the original coxph model had strata, should the
    predictions be for all strata, or only for those in newdata?
    The default for this is TRUE if either the \code{id} or
    \code{individual} argument is present, as these require strata for
    the resulting curves to make sense.
    Otherwise the default is FALSE, which means to ignore any strata
    variable in the newdata data set, and produce predicted survivals
    for the entire set of strata in the original model. In this case
    some components of the output will be matrices with one column for
    each row of newdata.}
  \item{na.action}{the na.action to be used on the newdata argument}
  
   \item{\dots}{for future methods}
}
\value{
an object of class \code{"survfit"}.  
See \code{survfit.object} for 
details. Methods defined for survfit objects are  
\code{print}, \code{plot}, 
\code{lines}, and \code{points}. 
}
\details{
    Serious thought has been given to removing the `default' for
    \code{newdata}, which is to use a single "psuedo" subject with
    covariate values equal to the means of the data set.
    The resulting curve(s) almost never make sense.
    It remains due to the unwarranted attachment to the option shown by
    some users and by other packages.  Two particularly egregious examples
    are factor variables and interactions.  Suppose one were studying
    interspecies transmission of a virus, and the data set has a factor
    variable with levels ("pig", "chicken") and about equal numbers of
    observations for each.  The ``mean'' covariate level will be 1/2 --
    is this a flying pig?  As to interactions assume data with sex coded as 0/1,
    ages ranging from 50 to 80, and a model with age*sex.  The ``mean''
    value for the age:sex interaction term will be about 30, a value
    that does not occur in the data.

    Users are strongly advised to use the newdata argument.  Note that
    this data set needs to contain values for the main effects but not
    for any interaction terms.

When the original model contains time-dependent covariates, then the
path of that covariate through time needs to be specified in order to
obtain a predicted curve. This requires \code{newdata} to contain
multiple lines for each hypothetical subject which gives the covariate
values, time interval, and strata for each line (a subject can change
strata), along with an \code{id} variable
 which demarks which rows belong to each subject.
The time interval must have the same (start, stop, status)
variables as the original model: although the status variable is not
used and thus can be set to a dummy value of 0 or 1, it is necessary for
the variables to be recognized as a \code{Surv} object.
Last, although predictions with a time-dependent covariate path can be
useful, it is very easy to create a prediction that is senseless.  Users
are encouraged to seek out a text that discusses the issue in detail.

When all the coefficients are zero,
the Kalbfleisch-Prentice estimator reduces to the Kaplan-Meier,
the Aalen estimate to the exponential of Nelson's
cumulative hazard estimate, and the Efron estimate to the
Fleming-Harrington estimate of survival.
The variances of the curves from a Cox model are larger than these
non-parametrec estimates, however, due
to the variance of the coefficients.

See \code{\link{survfit}} for more details about the counts (number of
events, number at risk, etc.)

The censor argument was fixed at FALSE in earlier versions of the code
and not made  available to the user.
The default argument is sensible in most instances --- and causes the
familiar + sign to appear on plots --- it is not sensible for time
dependent covariates since it may lead to a large number of spurious marks.
}
\section{References}{
Fleming, T. H. and Harrington, D. P. (1984).  Nonparametric estimation of the 
survival distribution in censored data.  \emph{Comm. in Statistics}  
\bold{13}, 2469-86. 

Kablfleisch, J. D. and Prentice, R. L. (1980).  
\emph{The Statistical Analysis of Failure Time Data.}
New York:Wiley. 

Link, C. L. (1984). Confidence intervals for the survival 
function using Cox's proportional hazards model with  
covariates.  \emph{Biometrics}  
\bold{40}, 601-610.

Therneau T and Grambsch P (2000), Modeling Survival Data: Extending the
Cox Model, Springer-Verlag.

Tsiatis, A. (1981). A large sample study of the estimate 
for the integrated hazard function in Cox's regression 
model for survival data. \emph{Annals of Statistics}  
\bold{9}, 93-108. 
}
\seealso{
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
fit1 <- survfit(Surv(stop, event=='progression') ~1, data=mgus1,
                    subset=(start==0))
fit2 <- survfit(Surv(stop, status) ~1, data=mgus1,
                    subset=(start==0), etype=event) #competing risks
# CI curves are always plotted from 0 upwards, rather than 1 down
plot(fit2, fun='event', xscale=365.25, xmax=7300, mark.time=FALSE,
            col=2:3, xlab="Years post diagnosis of MGUS")
lines(fit1, fun='event', xscale=365.25, xmax=7300, mark.time=FALSE,
            conf.int=FALSE)
text(10, .4, "Competing Risk: death", col=3)
text(16, .15,"Competing Risk: progression", col=2)
text(15, .30,"KM:prog")
}
\keyword{survival}
