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
  cluster, ymin, ymax, timewt= c("n", "S", "S/G", "n/G2", "I"),
  influence=0, ranks = FALSE, reverse=FALSE, timefix=TRUE, keepstrata=10, \ldots)
\method{concordance}{lm}(object, \ldots, newdata, cluster, ymin, ymax,
  influence=0, ranks=FALSE, timefix=TRUE, keepstrata=10)
\method{concordance}{coxph}(object, \ldots, newdata, cluster, ymin, ymax,
  timewt= c("n", "S", "S/G", "n/G2", "I"), influence=0,
  ranks=FALSE, timefix=TRUE, keepstrata=10)
\method{concordance}{survreg}(object, \ldots, newdata, cluster, ymin, ymax,
  timewt= c("n", "S", "S/G", "n/G2", "I"), influence=0,
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
    scaled ranks that make up the overall score.  
  }
  \item{reverse}{if TRUE then assume that larger \code{x} values predict
    smaller response values \code{y}; a proportional hazards model is
    the common example of this, larger hazard = shorter survival.}

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
  The concordance is an estimate of
  \eqn{Pr(x_i < x_j | y_i < y_j)}{Pr(x[i] < x[j] | y[i] < y[j])},
  for a model fit replace \eqn{x}{x} with \eqn{\hat y}{yhat}, the
  predicted response from the model.
  For a survival outcome some pairs of values
  are not comparable, e.g., censored at time 5 and a death at time 6,
  as we do not know if the first observation will or will not outlive
  the second.  In this case the total number of evaluable pairs is smaller.

  A naive compuation requires adding up over all n(n-1)/2 comparisons,
  which can be quite slow for large data sets.
  This routine uses an O(n log(n)) algorithm.
  At each uncensored event time y, compute the rank of x for the subject
  who had the event as compared to the x values for all others with a longer
  survival, where the rank has value between 0 and 1.
  The concordance is a weighted mean of these ranks,
  determined by the \code{timewt} option. The rank vector can be
  efficiently updated as subjects are added to the risk set.
  For further details see the vignette.

  The \code{timewt} option is only applicable to censored data.  In this
  case the default corresponds to Harrell's C statistic, which is
  closely related to the Gehan-Wilcoxon test;
  \code{timewt="S"} corrsponds to the Peto-Wilcoxon,
  \code{timewt="S/G"} is suggested by Schemper, and
  \code{timewt="n/G2"} corresponds to Uno's C.
  It turns out that the Schemper and Uno weights are computationally
  identical, we have retained both option labels as a user convenience.
  The \code{timewt= "I"} option is related to the log-rank
  statistic.
  
  Relatations to other statistics:
  For continuous x and y, 2C- 1 is equal to Somers' d.
  If the response is binary, C is equal to the area under the reciever
  operating curve or AUC.
  For a survival response and binary predictor C is the numerator of
  the Genan-Wilcoxon test.
  
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
     at each event time, along with the time weight, and the case weight of
     the observation.  The time weighted sum of the ranks will equal
     concordant pairs - discordant pairs.}
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
concordance(fit1, timewt="n/G2")  # Uno's weighting

# logistic regression
fit2 <- glm(pstat ~ age + sex + mspike, binomial, data= mgus2)
concordance(fit2)  # equal to the AUC

# compare two non-nested models
options(na.action = na.exclude)   # predict all 1341 obs, including missing
fit3 <- coxph(Surv(ptime, pstat) ~ age + sex + mspike, mgus2)
fit4 <- coxph(Surv(ptime, pstat) ~ age + sex + hgb + creat, mgus2)
            family= poisson, data= mgus2)
p3 <- -predict(fit3)  # high risk scores predict shorter survival
p4 <- -predict(fit4)
options(na.action = na.omit)      # return to the R default

cfit <- concordance(Surv(ptime, pstat) ~ p3 + p4, mgus2)
cfit
round(coef(cfit), 3)
signif(vcov(cfit), 3)

test <- c(1, -1)  # contrast vector for model 1 - model 2 concordance
round(c(difference = test %*% coef(cfit),
        sd= sqrt(test %*% vcov(cfit) %*%test)), 3)
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
