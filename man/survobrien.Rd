\name{survobrien}
\alias{survobrien}
\title{
 O'Brien's Test for Association of a Single Variable with Survival
}
\description{
Peter O'Brien's test for association of a single variable with survival 
This test is proposed in Biometrics, June 1978. 
}
\usage{
survobrien(formula, data,
    transform = function(r) (r - 0.5)/(0.5 + length(r) - r))
}
\arguments{
\item{formula}{
a valid formula for a cox model, without time dependent covariates. 
}
\item{data}{
a data frame. 
    }
\item{transform}{the transformation function to be applied at each
   time point. The default is O'Brien's suggestion.}
}
\value{
a new data frame.  The original time and status variables are removed, 
and have been replaced with \code{start}, \code{stop}, and \code{event}.  If a 
predictor variable is a factor or is protected with \code{I()}, it is 
retained as is.  Other predictor variables have been replaced with 
time-dependent logit scores. 

Because of the time dependent variables, the new data frame will have many 
more rows that the original data, approximately \#rows * \#deaths /2. 
}
\section{Method}{
A time-dependent cox model can now be fit to the new data. 
The univariate statistic, as originally proposed, is equivalent to 
single variable score tests from the time-dependent model. 
This equivalence is the rationale for using the time dependent model as a 
multivariate extension of the original paper. 

In O'Brien's method, the x variables are re-ranked at each death time.  A 
simpler method, proposed by Prentice, ranks the data only once at the 
start. The results are usually similar. 
}
\references{
O'Brien, Peter, "A Nonparametric Test for Association with Censored Data", 
\emph{Biometrics} 34: 243-250, 1978. 
}
\seealso{
\code{\link{survdiff}} 
}
\keyword{survival}
\examples{
xx <- survobrien(Surv(futime, fustat) ~ age + factor(rx) + I(ecog.ps), 
			       data=ovarian) 
coxph(Surv(start, stop, event) ~ age, data=xx) 
coxph(Surv(start, stop, event) ~ age + rx + ecog.ps, data=xx) 
}
