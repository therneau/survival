\name{yates_setup}
\alias{yates_setup}
\title{Method for adding new models to the \code{yates} function.
}
\description{This is a method which is called by the \code{yates}
  function, in order to setup the code to handle a particular
  model type.  Methods for glm, coxph, and default are part of
  the survival package.
}
\usage{
yates_setup(fit, \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{a fitted model object}

  \item{\ldots}{optional arguments for some methods}
}
\details{
If the predicted value should be the linear predictor, the function
should return NULL.  The \code{yates} routine has particularly efficient
code for this case.
Otherwise it should return a prediction function or a list of two
elements containing the prediction function and a summary function.
The prediction function will be passed the linear predictor as a single
argument and should return a vector of predicted values.
}

\author{Terry Therneau}
\note{See the vignette on population prediction for more details.}
\seealso{
   \code{\link{yates}}
}
\keyword{ models }
\keyword{ survival }
