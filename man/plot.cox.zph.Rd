\name{plot.cox.zph}
\alias{plot.cox.zph}
\title{
  Graphical Test of Proportional Hazards 
}
\description{
Displays a graph of the scaled Schoenfeld residuals, along with a smooth curve. 
}
\usage{
\method{plot}{cox.zph}(x, resid=TRUE, se=TRUE, df=4, nsmo=40, var,
        xlab="Time", ylab, lty=1:2, col=1, lwd=1, pch=1, cex=1,
    hr=FALSE, plot=TRUE, ...)
}
\arguments{
\item{x}{
result of the \code{cox.zph} function. 
}
\item{resid}{
a logical value, if \code{TRUE} the residuals are included on the plot, as well as the smooth fit. 
}
\item{se}{
a logical value, if \code{TRUE}, confidence bands at two standard errors 
will be added. 
}
\item{df}{
the degrees of freedom for the fitted natural spline, \code{df=2} leads 
to a linear fit. 
}

\item{nsmo}{number of points to use for the lines}

\item{var}{
the set of variables for which plots are desired.  By default, plots are 
produced in turn for each variable of a model.  Selection of a single variable 
allows other features to be added to the plot, e.g., a horizontal line at 
zero or a main title. 

This has been superseded by a subscripting method; see the example below. 
}
\item{hr}{if TRUE, label the y-axis using the estimated hazard ratio
  rather than the estimated coefficient.  (The plot does not change,
  only the axis label.)}

\item{xlab}{label for the x-axis of the plot}
\item{ylab}{optional label for the y-axis of the plot.  If missing a
  default label is provided.  This can be a vector of labels.}
\item{lty, col, lwd}{line type, color, and line width for the overlaid
  curve.  Each of these can be vector of length 2, in which case the
  second element is used for the confidence interval.}
\item{plot}{if FALSE, return a list containing the x and y values of
  the curve, instead of drawing a plot}
\item{pch}{used for points on the plot, see \code{points}}
\item{cex}{used for points on the plot, see \code{points}}
\item{\dots}{
additional graphical arguments passed to the \code{plot} function. 
}
}
\section{Side Effects}{
a plot is produced on the current graphics device. 
}
\seealso{
\code{\link{coxph}},  \code{\link{cox.zph}}.   
}
\examples{
vfit <- coxph(Surv(time,status) ~ trt + factor(celltype) + 
              karno + age, data=veteran, x=TRUE) 
temp <- cox.zph(vfit) 
plot(temp, var=3)      # Look at Karnofsy score, old way of doing plot 
plot(temp[3])     # New way with subscripting 
abline(0, 0, lty=3) 
# Add the linear fit as well  
abline(lm(temp$y[,3] ~ temp$x)$coefficients, lty=4, col=3)  
title(main="VA Lung Study") 
}
\keyword{survival}

