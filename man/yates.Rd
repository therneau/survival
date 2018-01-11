\name{yates}
\alias{yates}
\title{Population prediction}
\description{Compute the mean predicted value from a model fit, for
  a chosen population and statistic.
}
\usage{
yates(fit, term, population = c("data", "factorial", "sas"),
levels, test = c("global", "trend", "pairwise"), predict = "linear",
options, nsim = 1000, method = c("direct", "sgtt"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{a model fit.  Examples using lm, glm, and coxph objects
    are given in the vignette.
}
\item{term}{the term from the model whic is to be evaluated.
  This can be written as a character string or as a formula.
}
  \item{population}{the population to be used for the adjusting
    variables.  User can supply their own data frame or select one
    of the built in choices.
    The argument also allows "emprical" and "yates" as aliases for
    data and factorial, respectively, and ignores case.
}
  \item{levels}{optional, what values for \code{term} should be used.
}
  \item{test}{the test for comparing the population predictions.
}
\item{predict}{what to predict.  For a glm model this might be the
  'link' or 'response'.  For a coxph model it can be linear, risk, or
  survival.  User written functions are allowed.
}
  \item{options}{optional arguments for the prediction method.
}
\item{nsim}{number of simulations used to compute a variance for the
  predictions.  This is not needed for the linear predictor.
}
\item{method}{the computational approach for testing equality of the
  population predictions.  Either the direct approach or the algorithm
  used by the SAS glim procedure for "type 3" tests.
}
}

\details{
  The many options and details of this function are best described in a
  vignette on population prediction.
}
\value{an object of class \code{yates} with components of
  \item{estimate}{a data frame with one row for each level of the term,
    and columns containing the level, the mean population predicted
    value (mppv) and its standard deviation.}
  \item{tests}{a matrix giving the test statistics}
  \item{mvar}{the full variance-covariance matrix of the mppv values}
  \item{summary}{optional: any further summary if the values provided by
    the prediction method.}
}
\author{Terry Therneau}
\examples{
fit1 <- lm(skips ~ Solder*Opening + Mask, data = solder)
yates(fit1, ~Opening, population = "factorial")

fit2 <- coxph(Surv(time, status) ~ factor(ph.ecog)*sex + age, lung)
yates(fit2, ~ ph.ecog, predict="risk")  # hazard ratio
}
\keyword{ models }
\keyword{ survival }
