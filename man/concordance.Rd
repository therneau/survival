\name{concordance}
\alias{concordance.fit}
\title{Compute the concordance statistic for data or a model}
\description{
The concordance statistic compute the agreement between an observed
response and a predictor.  It is closely related to Kendall's tau-a and
tau-b, Goodman's gamma, and Somers' d, all of which can also be
calculated from the results of this function.
}
\usage{
concordance(x, \ldots)
\method{concordance}{formula}(x, data, weights, subset, na.action,
group, \ldots)
\method{concordance}{lm}(x, \ldots)
\method{concordance}{coxph}(x, reverse=TRUE, \ldots)
\method{concordance}{survreg}(x, ldots)
concordance.fit(y, x, strata, weight, ymin, ymax,
   timewt= c("n", "S", "S/G", "n/G", "n/G^2"), group,
   influence = 0, residuals = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a fitted model or formula.  The formula should be of
  the form \code{y ~x}  or \code{Surv(y1, y2) ~ x} with a single
  numeric or survival response and a single covariate.}

 \item{data}{
    a data.frame in which to interpret the variables named in 
    the \code{formula}, or in the \code{subset} and the \code{weights}
    argument. 
  }
  \item{weights}{
    vector of case weights, see the note below.
    For a thorough discussion of these see the
    book by Therneau and Grambsch.
  }
  \item{subset}{
    expression indicating which subset of the rows of data should be used in 
    the fit.    All observations are included by default. 
  }
  \item{na.action}{
    a missing-data filter function.  This is applied to the model.frame
    after any 
   subset argument has been used.  Default is \code{options()\$na.action}. 
  }

  \item{reverse}{if true the concordance is reversed: high values of
    \code{x} are assumed to predict low values of \code{y}}
  \item{group}{optional grouping vector, for calculating the robust
    variance}
  
  \item{y}{the response, which may be numeric or survival}

  \item{x}{the predictor variable}

  \item{strata}{optional stratification variable.  Counts of
    concordant, discordant and tied pairs 
    are computed separately per stratum, and then added.
    This will appear as \code{strata(z)} in a formula.
  }

   \item{ymin, ymax}{compute the concordance over the restricted range
     ymin <= y <= ymax.
}
    \item{timewt}{the weighting to be applied.  The overall value is a
      weighted mean over event times of the rank of the observation with
      the event, within all those at risk at that time
}
  \item{influence}{1= return the dfbeta vector, 2= return the full
    influence matrix.
}
  \item{residuals}{if TRUE, return a data frame containing the
     individual ranks that make up the overall score
}
}
\details{
  Users will rarely if ever call the concordance.fit routine directly,
  it is called by the formula and model methods to do the actual work.
}
\value{
  An object of class \code{concordance} containing the following
  components:
  \item{concordance}{the estimated concordance value}
  \item{var}{a vector containing the estimated variance of the
     concordance, first using the infinitesimal jackknife (IJ) method
     and the second using the variance of the score statistic from a
     proportional hazards model}
   \item{count}{a vector containing the number of concordant pairs,
     discordant, tied on x but not y, tied on y but not x, and tied on
     both x and y}
   \item{dfbeta}{optional, the vector of leverage estimates for the
     concordance}
   \item{influence}{optional, the matrix of leverage values for each of
     the counts, one row per observation}
   \item{residuals}{optional, a data frame containing the Somers' d rank
     at each event time, along with the time weight, case weight of the
     observation with an event, and variance (contribution to the
     proportional hazards model information matrix)}
 }
  
\author{Terry Therneau}
\seealso{\code{\link{coxph}}}
\examples{
fit1 <- coxph(Surv(time, status) ~ age + sex + ph.ecog, lung)
concordance(fit1, timewt="n")  # 

# logistic regression
fit2 <- glm(status ~ age + sex + ph.ecog, binomial, data=lung)
concordance(fit2)  # equal to the AUC
}

\keyword{ survival }
