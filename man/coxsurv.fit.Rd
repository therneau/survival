\name{coxsurv.fit}
\alias{coxsurv.fit}
\title{
A direct interface to the `computational engine' of survfit.coxph
}
\description{
  This program is mainly supplied to allow other packages to invoke the
  survfit.coxph function at a `data' level rather than a `user' level.
  It does no checks on the input data that is provided, which can lead
  to unexpected errors if that data is wrong.
}
\usage{
coxsurv.fit(ctype, stype, se.fit, varmat, cluster, 
            y, x, wt, risk, position, strata, oldid,
            y2, x2, risk2, strata2, id2, unlist=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stype}{survival curve computation: 1=direct, 2=exp(-cumulative
    hazard)}
  \item{ctype}{cumulative hazard computation: 1=Breslow, 2=Efron}
  \item{se.fit}{if TRUE, compute standard errors}
  
  \item{varmat}{the variance matrix of the coefficients
  }
  \item{cluster}{vector to control robust variance}
  
  \item{y}{the response variable used in the Cox model.  (Missing values
    removed of course.)
  }
  \item{x}{covariate matrix used in the Cox model
  }
  \item{wt}{weight vector for the Cox model. If the model was unweighted
    use a vector of 1s.
  }
  \item{risk}{the risk score exp(X beta + offset) from the fitted Cox model.}
  \item{position}{optional argument controlling what is counted as
    'censored'.  Due to time dependent covariates, for instance, a
    subject might have start, stop times of (1,5)(5,30)(30,100).  Times
    5 and 30 are not 'real' censorings.  Position is 1 for a real start,
    2 for an actual end, 3 for both, 0 for neither.}
  \item{strata}{strata variable used in the Cox model. This will be a
    factor.}
  \item{oldid}{identifier for subjects with multiple rows in the
    original data.}
  \item{y2, x2, risk2, strata2}{variables for the hypothetical subjects,
    for which prediction is desired}
  \item{id2}{optional; if present and not NULL this should be
    a vector of identifiers of length \code{nrow(x2)}.
    A non-null value signifies that \code{x2} contains time dependent
    covariates, in which case this identifies which rows of \code{x2} go
    with each subject.
  }
  \item{unlist}{if \code{FALSE} the result will be a list with one
    element for each strata.  Otherwise the strata are ``unpacked'' into
    the form found in a \code{survfit} object.}
}
\value{a list containing nearly all the components of a \code{survfit}
  object.  All that is missing is to add the confidence intervals, the
  type of the original model's response (as in a coxph object), and the
  class.
}
\note{The source code for for both this function and
  \code{survfit.coxph} is written using noweb.  For complete
  documentation see the \code{inst/sourcecode.pdf} file.
}  
\author{Terry Therneau}
\seealso{\code{\link{survfit.coxph}}
}
\keyword{survival}

