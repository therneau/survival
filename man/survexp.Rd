\name{survexp}
\alias{survexp}
\alias{print.survexp}
\title{
  Compute Expected Survival 
}
\description{
Returns either the expected survival of a cohort of subjects, or the 
individual expected survival for each subject. 
}
\usage{
survexp(formula, data, weights, subset, na.action, rmap, times, cohort=TRUE, 
        conditional=FALSE, ratetable=survexp.us, scale=1, npoints, 
        se.fit, model=FALSE, x=FALSE, y=FALSE)
}
\arguments{
\item{formula}{
formula object.  The response variable is a vector of follow-up times 
and is optional.  The predictors consist of optional grouping variables 
separated by the \code{+} operator (as in \code{survfit}), and is often
\code{~1}, i.e., expected survival for the entire group.
}
\item{data}{
data frame in which to interpret the variables named in 
the \code{formula}, \code{subset} and \code{weights} arguments. 
}
\item{weights}{
case weights. 
}
\item{subset}{
expression indicating a subset of the rows of \code{data} to be used in the fit. 
}
\item{na.action}{
function to filter missing data. This is applied to the model frame after  
\code{subset} has been applied.  Default is \code{options()$na.action}. A possible 
value for \code{na.action} is \code{na.omit}, which deletes observations that contain 
one or more missing values. 
}
\item{rmap}{
 an optional list that maps data set names to the ratetable names.  See
 the details section below.
 }
\item{times}{
vector of follow-up times at which the resulting survival curve is  
evaluated.  If absent, the result will be reported for each unique  
value of the vector of follow-up times supplied in \code{formula}. 
}
\item{cohort}{
logical value: if \code{FALSE}, each subject is treated as a subgroup of size 1. 
The default is \code{TRUE}. 
}
\item{conditional}{
logical value: if \code{TRUE}, the follow-up times supplied in \code{formula} 
are death times and conditional expected survival is computed. 
If \code{FALSE}, the follow-up times are potential censoring times.  
If follow-up times are missing in \code{formula}, this argument is ignored.   
}
\item{ratetable}{
a table of event rates, such as \code{survexp.uswhite}, or a fitted Cox model. 
}
\item{scale}{
numeric value to scale the results.  If \code{ratetable} is in units/day, 
\code{scale = 365.25} causes the output to be reported in years. 
}
\item{npoints}{
number of points at which to calculate intermediate results, evenly spaced  
over the range of the follow-up times.  The usual (exact) calculation is done  
at each unique follow-up time. For very large data sets specifying \code{npoints}  
can reduce the amount of memory and computation required. 
For a prediction from a Cox model \code{npoints} is ignored. 
}
\item{se.fit}{
compute the standard error of the predicted survival.  
The default is to compute standard errors whenever  
possible, which at this time is only for the Ederer method and a Cox  
model as the rate table. 
}
\item{model,x,y}{
flags to control what is returned.  If any of these is true, then the
model frame, the model matrix, and/or the vector of response times will be
returned as components of the final result, with the same names as the
flag arguments.
}}

\value{
if \code{cohort=TRUE} an object of class \code{survexp}, 
otherwise a vector of per-subject expected survival values. 
The former contains the number of subjects at risk 
and the expected survival for the cohort at each requested time. 
}
\details{
Individual expected survival is usually used in models or testing, to 
`correct' for the age and sex composition of a group of subjects. 
For instance, assume that birth date, entry date into the study, 
sex and actual survival time are all known for a group of subjects. 
The \code{survexp.us} population tables contain expected death rates 
based on calendar year, sex and age. 
Then 
\preformatted{
  haz <- -log(survexp(fu.time ~ 1, data=mydata, 
                      rmap = list(year=entry.dt, age=(birth.dt-entry.dt)),
                      cohort=FALSE)) 
}
gives for each subject the total hazard experienced up to their observed 
death time or last follow-up time (variable fu.time) 
This probability can be used as a rescaled time value in models: 
\preformatted{
glm(status ~ 1 + offset(log(haz)), family=poisson) 
glm(status ~ x + offset(log(haz)), family=poisson) 
}
In the first model, a test for intercept=0 is the one sample log-rank 
test of whether the observed group of subjects has equivalent survival to 
the baseline population.  The second model tests for an effect of variable 
\code{x} after adjustment for age and sex. 

The ratetable being used may have different variable names than the user's
data set, this is dealt with by the \code{rmap} argument.  
The rate table for the above calculation was \code{survexp.us}, a call to
\code{summary{survexp.us}} reveals that it expects to have variables 
\code{age} = age in days, \code{sex}, and \code{year} = the date of study
entry, we create them in the \code{rmap} line.  The sex variable is not
mapped, therefore the code assumes that it exists in \code{mydata} in the
correct format.  (Note: for factors such as sex, the program will match on
any unique abbreviation, ignoring case.)

Cohort survival is used to produce an overall survival curve.  This is then 
added to the Kaplan-Meier plot of the study group for visual comparison 
between these subjects and the population at large.  There are three common 
methods of computing cohort survival. 
In the "exact method" of Ederer the cohort is not censored; this corresponds 
to having no response variable in the formula.  Hakulinen recommends censoring 
the cohort at the anticipated censoring time of each patient, and Verheul 
recommends censoring the cohort at the actual observation time of each 
patient. 
The last of these is the conditional method. 
These are obtained by using the respective time values as the 
follow-up time or response in the formula. 
}
\references{
Berry, G. (1983). The analysis of mortality by the subject-years method. 
\emph{Biometrics}, 39:173-84.

Ederer, F., Axtell, L. and Cutler, S. (1961). 
The relative survival rate: a statistical methodology. 
\emph{Natl Cancer Inst Monogr}, 6:101-21.

Hakulinen, T. (1982). 
Cancer survival corrected for heterogeneity in patient withdrawal. 
\emph{Biometrics}, 38:933-942.

Verheul, H., Dekker, E., Bossuyt, P., Moulijn, A. and Dunning, A. (1993). 
Background mortality in clinical survival studies. 
\emph{Lancet}, 341: 872-875.
}
\seealso{
\code{\link{survfit}}, \code{\link{pyears}},  \code{\link{survexp.us}},  \code{\link{survexp.fit}}.   
}
\examples{
# 
# Stanford heart transplant data
#  We don't have sex in the data set, but know it to be nearly all males.
# Estimate of conditional survival  
survexp(futime ~ 1, rmap=list(sex="male", year=accept.dt,   
          age=(accept.dt-birth.dt)), conditional=TRUE, data=jasa) 
# Estimate of expected  survival stratified by prior surgery 
survexp(futime ~ surgery, rmap= list(sex="male", year=accept.dt,  
	age=(accept.dt-birth.dt)), conditional=TRUE, data=jasa) 

## Compare the survival curves for the Mayo PBC data to Cox model fit
## 
pfit <-coxph(Surv(time,status>0) ~ trt + log(bili) + log(protime) + age +
                platelet, data=pbc)
plot(survfit(Surv(time, status>0) ~ trt, data=pbc))
lines(survexp( ~ trt, ratetable=pfit, data=pbc), col='purple')
}
\keyword{survival}

