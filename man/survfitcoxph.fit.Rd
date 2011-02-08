\name{survfitcoxph.fit}
\alias{survfitcoxph.fit}
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
survfitcoxph.fit(y, x, wt, x2, risk, newrisk, strata, se.fit, survtype,
vartype, varmat, id, y2, strata2, unlist=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{the response variable used in the Cox model.  (Missing values
    removed of course.)
}
  \item{x}{covariate matrix used in the Cox model
}
\item{wt}{weight vector for the Cox model. If the model was unweighted
  use a vector of 1s.
}
  \item{x2}{matrix describing the hypothetical subjects for which a
    curve is desired.  Must have the same number of columns as \code{x}.
  }
  \item{risk}{the risk score exp(X beta) from the fitted Cox model.  If
 the model had an offset, include it in the argument to exp.
  }
  \item{newrisk}{risk scores for the hypothetical subjects
  }
  \item{strata}{strata variable used in the Cox model. This will be a factor.
  }
  \item{se.fit}{if \code{TRUE} the standard errors of the curve(s) are returned
  }
  \item{survtype}{1=Kalbfleish-Prentice, 2=Nelson-Aalen, 3=Efron.  It is
    usual to match this to the approximation for ties used in the
    \code{coxph} model: KP for `exact', N-A for `breslow' and Efron for `efron'.
  }
  \item{vartype}{1=Greenwood, 2=Aalen, 3=Efron
  }
  \item{varmat}{the variance matrix of the coefficients
  }
  \item{id}{optional; if present and not NULL this should be
    a vector of identifiers of length \code{nrow(x2)}.
    A mon-null value signifies that \code{x2} contains time dependent
    covariates, in which case this identifies which rows of \code{x2} go
    with each subject.
  }
  \item{y2}{survival times, for time dependent prediction.  It gives
    the time range (time1,time2] for each row of \code{x2}.  Note: this
    must be a Surv object and thus contains a status indicator, which is
    never used in the routine, however.  
  }
  \item{strata2}{vector of strata indicators for \code{x2}.  This must
    be a factor.
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

