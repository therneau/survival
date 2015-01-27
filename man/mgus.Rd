\name{mgus}
\alias{mgus}
\alias{mgus1}
\docType{data}
\title{Monoclonal gammapothy data}
\description{
Natural history of 241 subjects with monoclonal gammapothy of
undetermined significance (MGUS).
}
\usage{
mgus
mgus1
}
\format{
  mgus: A data frame with 241 observations on the following 12 variables.
  \tabular{ll}{
    id:\tab subject id \cr
    age:\tab age in years at the detection of MGUS \cr
    sex:\tab \code{male} or \code{female} \cr
    dxyr:\tab year of diagnosis \cr
    pcdx:\tab for subjects who progress to a plasma cell malignancy \cr
       \tab the subtype of malignancy: multiple myeloma (MM) is the \cr
       \tab most common, followed by amyloidosis (AM), macroglobulinemia (MA),\cr
      \tab and other lymphprolifative disorders (LP) \cr
    pctime:\tab days from MGUS until diagnosis of a plasma cell malignancy \cr
    futime:\tab days from diagnosis to last follow-up \cr
    death:\tab 1= follow-up is until death \cr
    alb:\tab albumin level at MGUS diagnosis \cr
    creat:\tab creatinine at MGUS diagnosis \cr
    hgb:\tab hemoglobin at MGUS diagnosis \cr
    mspike:\tab size of the monoclonal protien spike at diagnosis \cr
    }

   mgus1: The same data set in start,stop format. Contains the id, age, sex, and
  laboratory variable described above along with
  \tabular{ll}{ 
    start, stop:\tab sequential intervals of time for each
      subject \cr
    status:\tab =1 if the interval ends in an event \cr
    event:\tab a factor containing the event type: censor, death, or plasma cell malignancy \cr
    enum: \tab event number for each subject: 1 or 2
  }
} 
\details{
 Plasma cells are responsible for manufacturing immunoglobulins, an
important part of the immune defense. At any given time there are
estimated to be about \eqn{10^6} different immunoglobulins in the circulation
at any one time.  When a patient has a plasma cell malignancy the
distribuion will become dominated by a single isotype, the product of
the malignant clone, visible as a spike on a serum protein
electrophoresis. Monoclonal gammapothy of undertermined significance
(MGUS) is the presence of such a spike, but in a patient with no
evidence of overt malignancy.  This data set of 241 sequential subjects
at Mayo Clinic was the groundbreaking study defining the natural history
of such subjects.
Due to the diligence of the principle investigator 0 subjects have
been lost to follow-up.

Three subjects had MGUS detected on the day of death.  In data set
\code{mgus1} these subjects have the time to MGUS coded as .5 day before
the death in order to avoid tied times.

These data sets were updated in Jan 2015 to correct some small errors.
}
\source{
Mayo Clinic data courtesy of Dr. Robert Kyle.
}
\examples{
# Create the competing risk curves for time to first of death or PCM
sfit <- survfit(Surv(start, stop, event) ~ sex, mgus1, subset=(enum==1))
print(sfit)  # the order of printout is the order in which they plot

plot(sfit, xscale=365.25, lty=c(2,1,2,1), col=c(1,1,2,2),
     xlab="Years after MGUS detection", ylab="Proportion")
legend(0, .8, c("Death/male", "Death/female", "PCM/male", "PCM/female"),
       lty=c(1,1,2,2), col=c(2,1,2,1), bty='n')

title("Curves for the first of plasma cell malignancy or death")
# The plot shows that males have a higher death rate than females (no
# surprise) but their rates of conversion to PCM are essentially the same.
}
\references{
 R Kyle, Benign monoclonal gammopathy -- after 20 to 35 years of
 follow-up,
 Mayo Clinic Proc 1993; 68:26-36. 
}
\keyword{datasets}
\keyword{survival}
