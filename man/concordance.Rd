\name{concordance}
\alias{concordance}
\alias{concordance.coxph}
\alias{concordance.formula}
\alias{concordance.lm}
\alias{concordance.survreg}
\title{Compute the concordance statistic for data or a model}
\description{
The concordance statistic compute the agreement between an observed
response and a predictor.  It is closely related to Kendall's tau-a and
tau-b, Goodman's gamma, and Somers' d, all of which can also be
calculated from the results of this function.
}
\usage{
concordance(object, \ldots)
\method{concordance}{formula}(object, data, weights, subset, na.action,
  cluster, ymin, ymax, timewt= c("n", "S", "S/G", "n/G", "n/G2", "I"),
  influence=0, ranks = FALSE, reverse=FALSE, timefix=TRUE, keepstrata=10, \ldots)
\method{concordance}{lm}(object, \ldots, newdata, cluster, ymin, ymax,
  influence=0, ranks=FALSE, timefix=TRUE, keepstrata=10)
\method{concordance}{coxph}(object, \ldots, newdata, cluster, ymin, ymax,
  timewt= c("n", "S", "S/G", "n/G", "n/G2", "I"), influence=0,
  ranks=FALSE, timefix=TRUE, keepstrata=10)
\method{concordance}{survreg}(object, \ldots, newdata, cluster, ymin, ymax,
  timewt= c("n", "S", "S/G", "n/G", "n/G2", "I"), influence=0,
  ranks=FALSE, timefix=TRUE, keepstrata=10)
}

\arguments{
  \item{object}{a fitted model or a formula.  The formula should be of
  the form \code{y ~x}  or \code{y ~ x + strata(z)} with a single
  numeric or survival response and a single predictor.
  Counts of concordant, discordant and tied pairs 
  are computed separately per stratum, and then added.
}

 \item{data}{
    a data.frame in which to interpret the variables named in 
    the \code{formula}, or in the \code{subset} and the \code{weights}
    argument. Only applicable if \code{object} is a formula.
  }
  \item{weights}{
    optional vector of case weights.
    Only applicable if \code{object} is a formula.
  }
  \item{subset}{
    expression indicating which subset of the rows of data should be used in 
    the fit.   Only applicable if \code{object} is a formula.
  }
  \item{na.action}{
    a missing-data filter function.  This is applied to the model.frame
    after any subset argument has been used.  Default is
   \code{options()\$na.action}. Only applicable if \code{object} is a formula.
  }

  \item{\ldots}{multiple fitted models are allowed.  Only applicable if
    \code{object} is a model object.}
  
  \item{newdata}{optional, a new data frame in which to evaluate (but
    not refit) the models}
  
  \item{cluster}{optional grouping vector for calculating the robust
    variance}
  
  \item{ymin, ymax}{compute the concordance over the restricted range
     ymin <= y <= ymax.  (For survival data this is a time range.)
  }
  \item{timewt}{the weighting to be applied.  The overall statistic is a
      weighted mean over event times.
    }
  \item{influence}{1= return the dfbeta vector, 2= return the full
    influence matrix, 3 = return both
  }
  \item{ranks}{if TRUE, return a data frame containing the
    individual ranks that make up the overall score.  
  }
  \item{reverse}{if TRUE then assume that larger \code{x} values predict
    smaller response values \code{y}; a proportional hazards model is
    the common example of this. }

  \item{timefix}{correct for possible rounding error.  See the
    vignette on tied times for more explanation. Essentially, exact ties
    are an important part of the concordance computatation, but "exact"
    can be a subtle issue with floating point numbers.
  }
  \item{keepstrata}{either TRUE, FALSE, or an integer value.
    Computations are always done within stratum, then added. If the
    total number of strata greater than \code{keepstrata}, or
    \code{keepstrata=FALSE}, those subtotals are not kept in the output.
    }
}
\details{
  At each event time, compute the rank of the subject who had the
  event as compared to all others with a longer survival, where the
  rank is value between 0 and 1.  The concordance is a weighted mean
  of these values, determined by the \code{timewt} option.
  For uncensored data each unique response value is compared to all
  those which are larger.

  Using the default value for \code{timewt} gives the area
  under the receiver operating curve (AUC) for a binary response,
  and (d+1)/2 when y is continuous, where d is Somers' d.
  For a survival time, \code{timewt} of n gives Harrell's c-statistic,
  which is closely related to the Gehan-Wilcoxon test; here n is short
  for \eqn{n(t)}{n(t)}, the number at risk at each event time.
  S corresponds to the Peto-Wilcoxon, n/G^2 is the weight advocated
  by Uno, and S/G the weighting proposed by Schemper.
  Peto pointed out that \eqn{n(t) \approx n(0) S(t-) G(t-)}
  {n(t) approx= n(0)S(t-)G(t-)} where S and G are the survival and censoring
  distributions; this implies that the Schemper and Uno weightings
  should give very similar results. 
  
  When the number of strata is very large, such as in a conditional
  logistic regression for instance (\code{clogit} function), a much
  faster computation is available when the individual strata results
  are not retained; use \code{keepstrata=FALSE} or \code{keepstrata=0}
  to do so. In the general case the \code{keepstrata = 10}
  default simply keeps the printout managable: it retains and prints
  per-strata information if the number of strata is <= 10.
}
\value{
  An object of class \code{concordance} containing the following
  components:
  \item{concordance}{the estimated concordance value or values}
  \item{count}{a vector containing the number of concordant pairs,
     discordant, tied on x but not y, tied on y but not x, and tied on
     both x and y}
  \item{n}{the number of observations}
  \item{var}{a vector containing the estimated variance of the
    concordance based on the infinitesimal jackknife (IJ) method.
    If there are multiple models it contains the estimtated
    variance/covariance matrix.}
  \item{cvar}{a vector containing the estimated variance(s) of the
    concordance values, based on the variance formula for the associated
    score test from a proportional hazards model.  (This was the primary
    variance used in the \code{survConcordance} function.)}
   \item{dfbeta}{optional, the vector of leverage estimates for the
     concordance}
   \item{influence}{optional, the matrix of leverage values for each of
     the counts, one row per observation}
   \item{ranks}{optional, a data frame containing the Somers' d rank
     at each event time, along with the time weight, case weight of the
     observation with an event, and variance (contribution to the
     proportional hazards model information matrix).
     A weighted mean of the ranks equals Somer's d.}
 }

\note{A coxph model that has a numeric failure may have undefined
  predicted values, in which case the concordance will be NULL.

   Computation for an existing coxph model along with \code{newdata} has
  some subtleties with respect to extra arguments in the original call.
  These include
  \itemize{
    \item tt() terms in the model.  This is not supported with newdata.
    \item subset.  Any subset clause in the original call is ignored,
    i.e., not applied to the new data.
    \item strata() terms in the model.  The new data is expected to
    have the strata variable(s) found in the original data set,
    with concordance computed within strata.
    The levels of the strata variable need not be
    the same as in the original data.
    \item id or cluster directives.  This has not yet been sorted out.
  }
}

\author{Terry Therneau}
\seealso{\code{\link{coxph}}}
\examples{
fit1 <- coxph(Surv(ptime, pstat) ~ age + sex + mspike, mgus2)
concordance(fit1, timewt="n") 

# logistic regression
fit2 <- glm(pstat ~ age + sex + mspike, binomial, data= mgus2)
concordance(fit2)  # equal to the AUC
}
\references{
  F Harrell, R Califf, D Pryor, K Lee and R Rosati, 
  Evaluating the yield of medical tests, J Am Medical Assoc, 1982.

  R Peto and J Peto,
  Asymptotically efficient rank invariant test procedures (with
  discussion), J Royal Stat Soc A, 1972.
  
  M Schemper, Cox analysis of survival data with non-proportional
  hazard functions, The Statistician, 1992.

  H Uno, T Cai, M Pencina, R D'Agnostino and Lj Wei,
  On the C-statistics for evaluating overall adequacy of risk 
        prediction procedures with censored survival data,
	Statistics in Medicine, 2011.
}

\keyword{ survival }
