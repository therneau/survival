\name{pyears2html}
\alias{pyears2html}
\title{
  Create HTML output from pyears object
}
\description{
Create an html file with tables based on the results stored in a pyears object.
}
\usage{
pyears2html(object, filename=paste(deparse(substitute(object)), ".html", sep = ""), 
            title=NULL, exp=TRUE, r=False, rr=TRUE,
            ci.rr= False, conf.level= 0.95, digits=2)
}
\arguments{
\item{object}{
result from a call to the function pyears
}
\item{filename}{
character string naming the file to write to.  The default uses the name of
the object.
}
\item{title}{
character string used as table caption.
}
\item{exp}{
Logical: Default=TRUE, indicates whether or not the expected number of events
should be included in the table.
}
\item{r}{
Logical: Default=False, indicates whether or not the observed/person-years ratio
should be included in the table.
}
\item{rr}{
Logical: Default=TRUE, indicates whether or not the observed/expected ratio
should be included in the table.
}
\item{ci.rr }{
Logical: Default=False, indicates whether or not the confidence interval for
the observed/expected ratio
should be included in the table.
}
\item{conf.level }{
If ci.rr=TRUE, the confidence level used to construct the RR confidence interval.
}
\item{digits}{
number of significant digits to use when rounding the main variables.
}
}
\value{
HTML is written to the specified file and the table content is returned invisibly.
}
\section{Side Effects}{
an HTML file is created that includes multiple tables with the output from
the pyears function.
}
\details{
Currently at most 4 varibles can be included in the pyears function.
}
\seealso{
pyears, html.table
}
\examples{
# Relative death rates for MGUS patients, by time from diagnosis
# (MGUS is often diagnosed incidentally when a patient is being seen
# for something more serious.)
# Skip the 2 subjects with MGUS found at autopsy (futime=0)
temp.time <- tcut(rep(0, nrow(mgus)), c(0, 1, 5, 10, 100)*365.25,
       labels=c('0-1 year', '1-5 years', '5-10 years', '10+ years'))
temp.start <- as.Date(paste('19', mgus$dxyr, '/07/01',  sep=''))
pfit <- pyears(Surv(futime, death) ~ temp.time + sex +
                     ratetable(year=temp.start, age=age*365.25, sex=sex),
                  ratetable=survexp.mn, data=mgus, subset=futime>0 )
pyears2html(pfit)
}
\keyword{pyears, table}
% docclass is function
% Converted by Sd2Rd version 47442.
