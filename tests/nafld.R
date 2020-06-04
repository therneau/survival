#
# check of Surv2 on a more complex data set
#
# First, baseline and death data, with age and study time (days) as the time scales
library(survival)

ndata <- data.frame(nafld1[,1:7],
                    days= 0,  death=0, iage=nafld1$age)
ndata <- merge(ndata, with(nafld1, data.frame(id=id, days=futime, death=status)),
               all=TRUE)
ndata$nafld0 <- 1*(ndata$id == ndata$case.id)

# Add in lab tests.  Some subjects have multiple tests before enrollment, we
#  want only the most recent for "value at baseline"
thin <- function(df) {  
    df$days <- pmax(0, df$days)
    keep <- !duplicated(df[,c("id", "days")], fromLast=TRUE)
    df[keep,]
}
    
ndata <- merge(ndata, with(thin(subset(nafld2, test=="chol")),
                           data.frame(id=id, days=days, chol=value)), 
               all=TRUE)
ndata <- merge(ndata, with(thin(subset(nafld2, test=="hdl")),
                           data.frame(id=id, days= days, hdl=value)),
               all=TRUE, by=c("id", "days"))

# Add in the comorbidities of interest.  None of these 4 happen to have
#  duplicates (MI does, for instance).  
# 
ndata <- merge(ndata, with(subset(nafld3, event=="diabetes"),
                           data.frame(id=id, days=pmax(0,days), diabetes=1)),
                   all=TRUE, by=c("id", "days"))
           
ndata <- merge(ndata, with(subset(nafld3, event=="htn"),
                           data.frame(id=id, days=pmax(0,days), htn=1)),
                     all=TRUE, by=c("id", "days"))
         
ndata <- merge(ndata, with(subset(nafld3, event=="dyslipidemia"),
                           data.frame(id=id, days= pmax(0, days), dyslipid=1)),
                       all=TRUE, by=c("id", "days"))
             
ndata <- merge(ndata, with(subset(nafld3, event=="nafld"),
                           data.frame(id=id, days= pmax(0,days), nafld=1)),
                     all=TRUE, by=c("id", "days"))
# age as a time scale
ndata$age <- nafld1$age[match(ndata$id, nafld1$id)] + round(ndata$days,2)

event <- function(id, time, status, istate, cumulative=FALSE, firstonly=FALSE) {
    # do all the work on ordered data
    ord <- order(id, time)
    id2 <- id[ord]
    time2 <- time[ord]
    stat2 <- ifelse(is.na(status[ord]), 0, status[ord])
    firstid <- !duplicated(id)
    if (cumulative) {
        csum <- cumsum(stat2)
        indx <- match(id2, id2)
        cstat<- csum + stat2[indx] - csum[indx]  
        cstat[stat2==0] <- 0
    }
    else cstat <- stat2
          
    if (!missing(istate)) cstat[firstid] <- istate

    keep <- (firstid | (!is.na(stat2)& stat2 !=0))
    newdata <- data.frame(id=id2[keep], time=time2[keep], status=cstat[keep])
    if (firstonly) { # only keep the first endpoint
        dummy <- seq(along= newdata$id)
        count <- dummy - dummy[match(newdata$id, newdata$id)] 
        newdata <- newdata[count<2,]
    } 

    newdata
}
             
# create the 4-way endpoint, people can have 2 comorbidities show up on the
#  same day.
temp1 <- rowSums(ndata[,c('diabetes', 'htn', 'dyslipid')], na.rm=TRUE)

temp2 <- with(ndata, event(id, days, pmax(temp1, 4*death, na.rm=TRUE),
                           cumulative=TRUE))
temp3 <- with(temp2, data.frame(id=id, days=time,
                     state=factor(pmin(status, 4), -1:4,
                               c("censor", paste0(0:3, "comorbid"), "death"))))
ndata <- merge(ndata, temp3, all=TRUE,  by=c("id", "days"))
check1 <- survcheck(Surv2(days, state) ~ 1, id=id, ndata)


# Repeat the process with tmerge
tdata <- tmerge(nafld1[,1:7], nafld1, id=id, death=event(futime, status))
# add the lab tests
tdata <- tmerge(tdata, subset(nafld2, test="chol"), id=id,
                chol = tdc(days, value))
tdata <- tmerge(tdata, subset(nafld2, test="hdl"), id=id,
                hdl = tdc(days, value))
#
# add the number of comorbidities
tdata <- tmerge(tdata, subset(nafld3, event %in% c('diabetes', 'dyslipidemia',
                                                   'htn')),
                id=id, ecount = cumevent(days), pcount= cumtdc(days))
tdata <- tmerge(tdata, subset(nafld3, event == "nafld"), id=id, 
                nafld=tdc(days))

#After the last tmerge call, we can create a combined endpoint
temp <- with(tdata, ifelse(death==1, 4, ecount))
tdata$state <- factor(temp, 0:4, c("censor", paste0(1:3, "comorbid"), "death"))

# Add the current state for each
tdata$cstate <- factor(tdata$pcount, 0:3, paste0(0:3, "comorbid"))

# Check it out
check2 <- survcheck(Surv(tstart, tstop, state) ~ 1, id=id, istate=cstate, tdata)

all.equal(check1$transitions, check2$transitions)
all.equal(check1$events, check2$events)
