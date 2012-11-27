\name{lines.survfit}
\alias{lines.survfit}
\alias{points.survfit}
\title{
  Add Lines or Points to a Survival Plot 
}
\description{
Often used to add the expected survival curve(s) to a Kaplan-Meier plot 
generated with \code{plot.survfit}. 
}
\usage{
\method{lines}{survfit}(x, type="s", mark=3, col=1, lty=1,
lwd=1, mark.time=TRUE, 
	xscale=1,  firstx=0, firsty=1, xmax, fun, conf.int=FALSE,  ...)
\method{points}{survfit}(x, xscale, xmax, fun, ...)
}

\arguments{
\item{x}{
a survival object, generated from the \code{survfit} or \code{survexp} functions. 
}
\item{type}{
the line type, as described in \code{lines}.  The default is a step function 
for \code{survfit} objects, and a connected line for \code{survexp} objects. 
}
\item{mark, col, lty, lwd}{
vectors giving the mark symbol, color, line type and line width for the 
added curves. 
}
\item{...}{other graphical parameters}
\item{mark.time}{
controls the labeling of the curves.   
If \code{FALSE}, no labeling is done.   
If \code{TRUE}, then curves are marked at each censoring time.   
If \code{mark.time} is a numeric vector, then curves are marked at  
the specified time points. 
}
\item{xscale}{
a number used to divide the x values.  If time was originally in days, a 
value of 365.25 would give a plotted scale in years. 
}
\item{firstx, firsty}{
the starting point for the survival curves.  If either of these is set to 
\code{NA} or < blank > the plot will start at the first time point of the curve. 
}
\item{xmax}{
the maximum horizontal plot coordinate.     
This shortens the curve before plotting it, so unlike using the 
\code{xlim} graphical parameter, warning messages about out of bounds points are 
not generated. 
}
\item{fun}{
an arbitrary function defining a transformation of the survival curve. 
For example \code{fun=log} is an alternative way to draw a log-survival curve 
(but with the axis labeled with log(S) values). 
Four often used transformations can be specified with a character 
argument instead: "log" is the same as using the \code{log=T} option, 
"event" plots cumulative events (f(y) = 1-y), 
"cumhaz" plots the cumulative hazard function (f(y) = -log(y)) 
and "cloglog" creates a complimentary log-log survival plot  
(f(y) = log(-log(y) along with log scale for the x-axis). 
}
\item{conf.int}{
if \code{TRUE}, confidence bands for the curves are also plotted. 
If set to \code{"only"}, then only the CI bands are plotted, and the curve 
itself is left off.   
This can be useful for fine control over the colors or line types of a plot. 
}
}
\value{
a list with components \code{x} and \code{y}, containing the coordinates of the 
last point on each of the curves (but not of the confidence limits). 
This may be useful for labeling. 
}
\section{Side Effects}{
one or more curves are added to the current plot. 
}
\seealso{
\code{\link{lines}},  \code{\link{par}},  \code{\link{plot.survfit}},  \code{\link{survfit}},  \code{\link{survexp}}.   
}
\details{
  When the \code{survfit} function creates a multi-state survival curve
  the resulting object has class `survfitms'.  The only difference in
  the plots is that that it defaults to a curve that goes from lower
  left to upper right (starting at 0), where survival curves default
  to starting at 1 and going down.  All other options are identical.
}    
\examples{
fit <- survfit(Surv(time, status==2) ~ sex, pbc,subset=1:312)
plot(fit, mark.time=FALSE, xscale=365.25,
        xlab='Years', ylab='Survival')
lines(fit[1], lwd=2, xscale=365.24)    #darken the first curve and add marks


# Add expected survival curves for the two groups,
#   based on the US census data
# The data set does not have entry date, use the midpoint of the study
efit <- survexp(~ ratetable(sex=sex,age=age*365.35,year=as.Date('1979/1/1')) +
            sex, data=pbc,  times=(0:24)*182)
temp <- lines(efit, lty=2, xscale=365.24, lwd=2:1)
text(temp, c("Male", "Female"), adj= -.1) #labels just past the ends
title(main="Primary Biliary Cirrhosis, Observed and Expected")

}
\keyword{survival}
