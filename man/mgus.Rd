\name{mgus}
\alias{mgus}
\alias{mgus1}
\alias{mgus2}
\docType{data}
\title{Monoclonal gammapothy data}
\description{
Natural history of 241 subjects with monoclonal gammapothy of
undetermined significance (MGUS).
}
\usage{
mgus
mgus1
mgus2
}
\format{
  mgus: A data frame with 241 observations on the following 12 variables.
  \tabular{ll}{
    id:\tab subject id \cr
    age:\tab age in years \cr
    sex:\tab \code{male} or \code{female} \cr
    dxyr:\tab year of diagnosis \cr
    pcdx:\tab for subjects who progress to a plasma cell malignancy \cr
       \tab the subtype of malignancy: multiple myeloma (MM) is the \cr
       \tab most common, followed by amyloidosis (AM), macroglobulinemia (MA),\cr
      \tab and other lymphprolifative (LP) \cr
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
    event:\tab the event type \cr
  }

  mgus2: The mgus data, but formatted in the competing risks style.  Each
  subject has three observations, one for time to death, one for time to
  MM, and one for time to a PC malignancy other than MM.
  Contains the id, age, sex, and
  laboratory variable described above along with 
    \tabular{ll}{
    time:\tab time to event or censoring \cr
    status:\tab 1 if the event occured, 0 otherwise \cr
    event:\tab death, myeloma, or other \cr
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
}
\source{
Mayo Clinic data courtesy of Dr. Robert Kyle.
}
\references{
 R Kyle, Benign monoclonal gammopathy -- after 20 to 35 years of
 follow-up,
 Mayo Clinic Proc 1993; 68:26-36. 
}
\keyword{datasets}
\keyword{survival}
