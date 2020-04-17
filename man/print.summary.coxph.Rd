\name{print.summary.coxph}
\alias{print.summary.coxph}
\title{
  Print method for summary.coxph objects
}
\description{
  Produces a printed summary of a fitted coxph model 
}
\usage{
\method{print}{summary.coxph}(x, digits=max(getOption("digits") - 3, 3),
  signif.stars = getOption("show.signif.stars"), expand=FALSE, ...)
}
\arguments{
  \item{x}{
    the result of a call to \code{summary.coxph} 
  }
  \item{digits}{significant digits to print}
  \item{signif.stars}{
    Show stars to highlight small p-values
  }
  \item{expand}{if the summary is for a multi-state coxph fit, print the
    results in an expanded format.}
  \item{\dots}{For future methods}
}
