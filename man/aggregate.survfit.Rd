\name{aggregate.survfit}
\alias{aggregate.survfit}
\title{Average survival curves}
\description{
For a survfit object containing multiple curves, create average curves
over a grouping.
}
\usage{
\method{aggregate}{survfit}(x, by = NULL, FUN = mean, ...)
}
\arguments{
  \item{x}{a \code{survfit} object which has a data dimension.}
  \item{by}{an optional list or vector of grouping elements, each as
    long as \code{dim(x)['data']}. }
  \item{FUN}{a function to compute the summary statistic of interest. }
  \item{\dots}{optional further arguments to FUN.}
}
\details{
  The primary use of this is to take an average over multiple survival
  curves that were created from a modeling function.  That is, a
  marginal estimate of the survival.
  It is primarily used to average over multiple predicted curves from a
  Cox model.
}
\value{a \code{survfit} object of lower dimension.}
\seealso{\code{\link{\survfit}}}
\examples{
cfit <- coxph(Surv(futime, death) ~ sex + age*hgb, data=mgus2)
# marginal effect of sex, after adjusting for the others
dummy <- rbind(mgus2, mgus2)
dummy$sex <- rep(c("F", "M"), each=nrow(mgus2)) # population data set
dummy <- na.omit(dummy)   # don't count missing hgb in our "population
csurv <- survfit(cfit, newdata=dummy)
dim(csurv)  # 2 * 1384 survival curves
csurv2 <- aggregate(csurv, dummy$sex)
\keyword{ survival }
