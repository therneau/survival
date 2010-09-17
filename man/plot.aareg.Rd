\name{plot.aareg}
\alias{plot.aareg}
\title{
  Plot an aareg object.
}
\description{
Plot the estimated coefficient function(s) from a fit
of Aalen's additive regression model.
}
\usage{
\method{plot}{aareg}(x, se=TRUE, maxtime, type='s', ...)
}
\arguments{
  \item{x}{
    the result of a call to the \code{aareg} function
  }
  \item{se}{
    if TRUE, standard error bands are included on the plot
  }
  \item{maxtime}{
    upper limit for the x-axis.
  }
  \item{type}{
    graphical parameter for the type of line, default is "steps".
  }
  \item{\dots }{
    other graphical parameters such as line type, color, or axis labels.
  }
}
\section{Side Effects}{
A plot is produced on the current graphical device.
}
\section{References}{
Aalen, O.O. (1989). A linear regression model for the analysis of life times.
Statistics in Medicine, 8:907-925.
}
\seealso{
aareg
}

