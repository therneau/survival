\name{survregDtest}
\alias{survregDtest}
\title{Verify a survreg distribution}
\description{
  This routine is called by \code{survreg} to verify that a distribution
object is valid.  
}
\usage{
survregDtest(dlist, verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dlist}{the list describing a survival distribution}
  \item{verbose}{return a simple TRUE/FALSE from the test for validity
(the default), or a verbose description of any flaws.}
}
\details{
  If the \code{survreg} function rejects your user-supplied distribution
as invalid, this routine will tell you why it did so.
}
\value{
  TRUE if the distribution object passes the tests, and either FALSE or a
vector of character strings if not.
}
\author{Terry Therneau}
\seealso{\code{\link{survreg.distributions}}, \code{\link{survreg}}}
\examples{
# An invalid distribution (it should have "init =" on line 2)
#  surveg would give an error message
mycauchy <- list(name='Cauchy',
                 init<- function(x, weights, ...) 
                      c(median(x), mad(x)),
                 density= function(x, parms) {
                      temp <- 1/(1 + x^2)
                      cbind(.5 + atan(temp)/pi, .5+ atan(-temp)/pi,
                            temp/pi, -2 *x*temp, 2*temp^2*(4*x^2*temp -1))
                      },
                 quantile= function(p, parms) tan((p-.5)*pi),
                 deviance= function(...) stop('deviance residuals not defined')
                 )

survregDtest(mycauchy, TRUE)
}
\keyword{survival}
