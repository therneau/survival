\name{frailty}
\alias{frailty}
\alias{frailty.gamma}
\alias{frailty.gaussian}
\alias{frailty.t}
\title{
  Random effects terms
}
\description{
The frailty function allows one to add a simple random effects term to
a Cox or survreg model.
}
\usage{
frailty(x, distribution="gamma", ...)
frailty.gamma(x, sparse = (nclass > 5), theta, df, eps = 1e-05, method = c("em","aic", "df", "fixed"), ...) 
frailty.gaussian(x, sparse = (nclass > 5), theta, df, method =c("reml","aic", "df", "fixed"), ...)
frailty.t(x, sparse = (nclass > 5), theta, df, eps = 1e-05, tdf = 5,method = c("aic", "df", "fixed"), ...)
}

\arguments{
  \item{x}{
    the variable to be entered as a random effect.  
    It is always treated as a factor.
  }
  \item{distribution}{
    either the \code{gamma}, 
    \code{gaussian} or \code{t}
    distribution may be specified.
    The routines \code{frailty.gamma},
    \code{frailty.gaussian} and 
    \code{frailty.t} do the actual work.
  }
  \item{\dots}{Arguments for specific distribution, including (but not
    limited to) }
  \item{sparse}{
    cutoff for using a sparse coding of the data matrix.  
    If the total number of levels of \code{x} is larger
    than this value, then a sparse matrix approximation is used.
    The correct cutoff is still a matter of exploration: if the number of
    levels is very large (thousands) then the non-sparse calculation may not be
    feasable in terms of both memory and compute time.  
    Likewise, the accuracy of the sparse approximation appears to be related to
    the maximum proportion of subjects in any one class, being best when no
    one class has a large membership.
  }
  \item{theta}{
    if specified, this fixes the variance of the random effect.
    If not, the variance is a parameter, and a best solution is sought.
    Specifying this implies \code{method='fixed'}.
  }
  \item{df}{
    if specified, this fixes the degrees of freedom for the random effect.
    Specifying this implies \code{method='df'}.
    Only one of \code{theta} or 
    \code{df} should be specified.
  }
  \item{method}{
    the method used to select a solution for theta, the variance of the
    random effect.  
    The \code{fixed} corresponds to a user-specified
    value, and no iteration is done.
    The \code{df} selects the variance such that the
    degrees of freedom for the random effect matches a user specified value.
    The \code{aic} method seeks to 
    maximize Akiake's information criteria 
    2*(partial likelihood - df).
    The \code{em} and \code{reml}
    methods are specific to Cox models with gamma and gaussian random effects,
    respectively.
    Please see further discussion below.
  }
  \item{tdf}{
    the degrees of freedom for the t-distribution.
  }
  \item{eps}{
    convergence critera for the iteration on theta.
  }
}
\value{
this function is used in the model statment of either
\code{coxph} or \code{survreg}.
It's results are used internally.
}
\details{
The \code{frailty} plugs into the general penalized
modeling framework provided by the \code{coxph} 
and \code{survreg} routines.  
This framework deals with likelihood, penalties, and degrees of freedom;
these aspects work well with either parent routine.

Therneau, Grambsch, and Pankratz show how maximum likelihood estimation for
the Cox model with a gamma frailty can be accomplished using a general
penalized routine, and Ripatti and Palmgren work through a similar argument
for the Cox model with a gaussian frailty.  Both of these are specific to
the Cox model.  
Use of gamma/ml or gaussian/reml with 
\code{survreg} does not lead to valid results. 

The extensible structure of the penalized methods is such that the penalty
function, such as \code{frailty} or
\code{pspine}, is completely separate from the modeling
routine.  The strength of this is that a user can plug in any penalization
routine they choose.  A weakness is that it is very difficult for the
modeling routine to know whether a sensible penalty routine has been
supplied.

Note that use of a frailty term implies a mixed effects model and use of
a cluster term implies a GEE approach; these cannot be mixed.

The \code{coxme} package has superseded
this method.  It is faster, more stable, and more flexible.
}
\section{References}{
S Ripatti and J Palmgren, Estimation of multivariate frailty 
models using penalized partial likelihood,
Biometrics, 56:1016-1022, 2000.

T Therneau, P Grambsch and VS Pankratz,
Penalized survival models and frailty,
J Computational and Graphical Statistics, 12:156-175, 2003.
}
\seealso{
\link{coxph}, \link{survreg}
}
\examples{
# Random institutional effect
coxph(Surv(time, status) ~ age + frailty(inst, df=4), lung)

# Litter effects for the rats data
rfit2a <- survreg(Surv(time, status) ~ rx +
                  frailty.gaussian(litter, df=13, sparse=FALSE), rats )
rfit2b <- survreg(Surv(time, status) ~ rx +
                  frailty.gaussian(litter, df=13, sparse=TRUE), rats )
}
\keyword{survival}

