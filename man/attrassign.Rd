\name{attrassign}
\alias{attrassign.default}
\alias{attrassign}
\alias{attrassign.lm}
\title{Create new-style "assign" attribute}
\description{
The \code{"assign"} attribute on model matrices describes which columns
come from which terms in the model formula. It has two versions. R uses
the original version, but the alternate version found in
S-plus is sometimes useful.
}
\usage{
\method{attrassign}{default}(object, tt,...)
\method{attrassign}{lm}(object,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{model matrix or linear model object}
  \item{tt}{terms object}
  \item{...}{ignored}
}

\value{
 A list with names corresponding to the term names and elements
 that are vectors indicating which columns come from which terms
}

\details{
  For instance consider the following
  \preformatted{
    survreg(Surv(time, status) ~ age + sex + factor(ph.ecog), lung)
  }
  R gives the compact for for assign, a vector (0, 1, 2, 3, 3, 3);
  which can be
  read as ``the first column of the X matrix (intercept) goes with none
  of the terms, the second column of X goes with term 1 of the model
  equation, the third column of X with term 2, and columns 4-6 with
  term 3''.
 
  The alternate (S-Plus default) form is a list
  \preformatted{
       $(Intercept)     1
       $age             2
       $sex             3
       $factor(ph.ecog) 4 5 6
     }
   }
   
\seealso{\code{\link{terms}},\code{\link{model.matrix}}}

\examples{
formula <- Surv(time,status)~factor(ph.ecog)
tt <- terms(formula)
mf <- model.frame(tt,data=lung)
mm <- model.matrix(tt,mf)
## a few rows of data
mm[1:3,]
## old-style assign attribute
attr(mm,"assign")
## alternate style assign attribute
attrassign(mm,tt)
}
\keyword{models}
