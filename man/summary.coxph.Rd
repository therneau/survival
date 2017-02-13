\name{summary.coxph}
\alias{summary.coxph}
\title{
  Summary method for Cox models 
}
\description{
Produces a summary of a fitted coxph model 
}
\usage{
\method{summary}{coxph}(object, conf.int=0.95, scale=1,...)
}
\arguments{
  \item{object}{
    the result of a coxph fit 
  }
  \item{conf.int}{
    level for computation of the confidence intervals. 
    If set to FALSE no confidence intervals are printed 
  }
  \item{scale}{
    vector of scale factors for the coefficients, defaults to 1. 
    The confidence limits are for the risk change associated with one scale unit. 
  }
  \item{\dots}{for future methods}
}
\value{
An object of class \code{summary.coxph}. 
}
\seealso{
  coxph, print.coxph  
}
\examples{
fit <- coxph(Surv(time, status) ~ age + sex, lung) 
summary(fit) 
\dontrun{
Call: 
coxph(formula = Surv(time, status) ~ age + sex, data = lung) 
  
  n= 228  
  
      coef exp(coef) se(coef)     z      p  
age  0.017     1.017  0.00922  1.85 0.0650 
sex -0.513     0.599  0.16745 -3.06 0.0022 
  
    exp(coef) exp(-coef) lower .95 upper .95  
age     1.017      0.983     0.999     1.036 
sex     0.599      1.670     0.431     0.831 
  
Rsquare= 0.06   (max possible= 0.999 ) 
Likelihood ratio test= 14.1  on 2 df,   p=0.000857 
Wald test            = 13.5  on 2 df,   p=0.00119 
Score (logrank) test = 13.7  on 2 df,   p=0.00105 
}}
\keyword{survival}
% docclass is function
% Converted by Sd2Rd version 37351.
