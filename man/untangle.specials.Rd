\name{untangle.specials}
\alias{untangle.specials}
\title{
Help Process the `specials' Argument of the `terms' Function.
}
\description{
Given a \code{terms} structure and a desired special name, this returns an
index appropriate for subscripting the \code{terms} structure and another
appropriate for the data frame.
}
\usage{
untangle.specials(tt, special, order=1)
}
\arguments{
\item{tt}{
a \code{terms} object.
}
\item{special}{
the name of a special function, presumably used in the terms object.
}
\item{order}{
the order of the desired terms.  If set to 2, interactions with the special
function will be included.
}}
\value{
a list with two components:
\item{vars}{
a vector of variable names, as would be found in the data frame, of the
specials.
}
\item{terms}{
a numeric vector, suitable for subscripting the terms structure, that indexes
the terms in the expanded model formula which involve the special.
}}
\examples{
formula<-Surv(tt,ss)~x+z*strata(id)
tms<-terms(formula,specials="strata")
## the specials attribute
attr(tms,"specials")
## main effects 
untangle.specials(tms,"strata")
## and interactions
untangle.specials(tms,"strata",order=1:2)
}
\keyword{survival}
% Converted by Sd2Rd version 0.3-2.
