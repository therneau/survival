\name{pyears}
\alias{pyears}
\title{
  Person Years 
}
\description{
This function computes the person-years of follow-up time contributed by a 
cohort of subjects, stratified into subgroups. 
It also computes the number of subjects who contribute to each cell of the 
output table, and optionally the number of events and/or expected number of 
events in each cell. 
}
\usage{
pyears(formula, data, weights, subset, na.action,  rmap,
       ratetable, scale=365.25, expect=c('event', 'pyears'),
       model=FALSE, x=FALSE, y=FALSE, data.frame=FALSE)
}
\arguments{
  \item{formula}{
    a formula object. 
    The response variable will be a vector of follow-up times 
    for each subject, or a \code{Surv} object containing the survival
    time and an event indicator. 
    The predictors consist of optional grouping variables 
    separated by + operators (exactly as in \code{survfit}),
    time-dependent grouping 
    variables such as age (specified with \code{tcut}), and optionally a 
    \code{ratetable} term. 
    This latter matches each subject to his/her expected cohort. 
  }
  \item{data}{
    a data frame in which to interpret the variables named in 
    the \code{formula}, or in the \code{subset} and the \code{weights}
    argument. 
  }
  \item{weights}{
    case weights. 
  }
  \item{subset}{
    expression saying that only a subset of the rows of the data 
    should be used in the fit. 
  }
  \item{na.action}{
    a missing-data filter function, applied to the model.frame, after any 
    \code{subset} argument has been used.
    Default is \code{options()$na.action}. 
  }
  \item{rmap}{
    an optional list that maps data set names to the ratetable names.  See
    the details section below.
  }
  \item{ratetable}{
    a table of event rates, such as \code{survexp.uswhite}. 
  }
  \item{scale}{
    a scaling for the results.  As most rate tables are in units/day, the 
    default value of 365.25 causes the output to be reported in years. 
  }
  \item{expect}{
    should the output table include the expected number of events, or the
    expected number of person-years of observation.  This is only valid with
    a rate table.
  }
  \item{data.frame}{
    return a data frame rather than a set of arrays.}
  \item{model, x, y}{
    If any of these is true, then the
    model frame, the model matrix, and/or the vector of response times will be
    returned as components of the final result.
  }
}
\value{
a list with components: 

\item{pyears}{
an array containing the person-years of exposure. (Or other units, depending 
on the rate table and the scale).
The dimension and dimmanes of the array correspond to the variables on
the right hand side of the model equation.
}
\item{n}{
an array containing the number of subjects who contribute time to each cell 
of the \code{pyears} array. 
}
\item{event}{
an array containing the observed number of events.  This will be present only 
if the response variable is a \code{Surv} object. 
}
\item{expected}{
an array containing the expected number of events (or person years if 
\code{expect ="pyears"}).  This will be present only if there was a 
\code{ratetable} term. 
}
\item{data}{
if the \code{data.frame} option was set, a data frame containing the
variables \code{n}, \code{event}, \code{pyears} and \code{event} that
supplants the four arrays listed above,
along with variables corresponding to each dimension.
There will be one row for each cell in the arrays.}
\item{offtable}{
the number of person-years of exposure in the cohort that was not part of 
any cell in the \code{pyears} array.  This is often useful as an error check;
if 
there is a mismatch of units between two variables, nearly all the person 
years may be off table. 
}
\item{summary}{
a summary of the rate-table matching. This is also useful as an error 
check. 
}
\item{call}{
an image of the call to the function. 
}
\item{na.action}{
the \code{na.action} attribute contributed by an \code{na.action}
routine, if any. 
}
}
\details{
Because \code{pyears} may have several time variables, it is necessary that all 
of them be in the same units.  For instance, in the call 
\preformatted{
  py <- pyears(futime ~ rx, rmap=list(age=age, sex=sex, year=entry.dt),
                    ratetable=survexp.us) 
}
the natural unit of the ratetable is hazard per day, it is important that
\code{futime}, 
\code{age} and \code{entry.dt} all be in days. 
Given the wide range of possible inputs, 
it is difficult for the routine to do sanity checks of this aspect. 

The ratetable being used may have different variable names than the user's
data set, this is dealt with by the \code{rmap} argument.  
The rate table for the above calculation was \code{survexp.us}, a call to
\code{summary{survexp.us}} reveals that it expects to have variables 
\code{age} = age in days, \code{sex}, and \code{year} = the date of study
entry, we create them in the \code{rmap} line.  The sex variable is not
mapped, therefore the code assumes that it exists in \code{mydata} in the
correct format.  (Note: for factors such as sex, the program will match on
any unique abbreviation, ignoring case.)

A special function \code{tcut} is needed to specify time-dependent cutpoints. 
For instance, assume that age is in years, and that the desired final 
arrays have as one of their margins the age groups 0-2, 2-10, 10-25, and 25+. 
A subject who enters the study at age 4 and remains under observation for 
10 years will contribute follow-up time to both the 2-10 and 10-25 
subsets.  If \code{cut(age, c(0,2,10,25,100))} were used in the formula, the 
subject would be classified according to his starting age only. 
The \code{tcut} function has the same arguments as \code{cut}, 
but produces a different output object which allows the \code{pyears} function 
to correctly track the subject. 

The results of \code{pyears} are normally used as input to further calculations. 
The \code{print} routine, therefore, is designed to give only a summary of the 
table. 
}
\seealso{
\code{\link{ratetable}},  \code{\link{survexp}},  \code{\link{Surv}}.   
}
\examples{
# Look at progression rates jointly by calendar date and age
# 
temp.yr  <- tcut(mgus$dxyr, 55:92, labels=as.character(55:91)) 
temp.age <- tcut(mgus$age, 34:101, labels=as.character(34:100))
ptime <- ifelse(is.na(mgus$pctime), mgus$futime, mgus$pctime)
pstat <- ifelse(is.na(mgus$pctime), 0, 1)
pfit <- pyears(Surv(ptime/365.25, pstat) ~ temp.yr + temp.age + sex,  mgus,
     data.frame=TRUE) 
# Turn the factor back into numerics for regression
tdata <- pfit$data
tdata$age <- as.numeric(as.character(tdata$temp.age))
tdata$year<- as.numeric(as.character(tdata$temp.yr))
fit1 <- glm(event ~ year + age+ sex +offset(log(pyears)),
             data=tdata, family=poisson)
\dontrun{
# fit a gam model 
gfit.m <- gam(y ~ s(age) + s(year) + offset(log(time)),  
                        family = poisson, data = tdata) 
}

# Example #2  Create the hearta data frame: 
hearta <- by(heart, heart$id,  
             function(x)x[x$stop == max(x$stop),]) 
hearta <- do.call("rbind", hearta) 
# Produce pyears table of death rates on the surgical arm
#  The first is by age at randomization, the second by current age
fit1 <- pyears(Surv(stop/365.25, event) ~ cut(age + 48, c(0,50,60,70,100)) + 
       surgery, data = hearta, scale = 1)
fit2 <- pyears(Surv(stop/365.25, event) ~ tcut(age + 48, c(0,50,60,70,100)) + 
       surgery, data = hearta, scale = 1)
fit1$event/fit1$pyears  #death rates on the surgery and non-surg arm

fit2$event/fit2$pyears  #death rates on the surgery and non-surg arm
}
\keyword{survival}
