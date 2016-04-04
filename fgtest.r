#
# Try out several models
#
library(survival)
library(mstate)
library(cmprsk)

tdata <- mgus2

tdata$etime <- with(mgus2, ifelse(pstat==0, futime, ptime))
tdata$event <- with(mgus2, ifelse(pstat==0, 2*death, 1))
tdata$event <- factor(tdata$event, 0:2, labels=c("censor", "pcm", "death"))
tdata$id <- 1:nrow(tdata)
tdata$age <- tdata$age/10   # use decades = nicer size coefficient


# Fine-Gray
fstat <- as.numeric(tdata$event) -1
xmat <-  with(tdata, cbind(sex= ifelse(sex=='F', 0, 1), 
                           age = age, mspike=mspike))
                           
fit1 <- crr(tdata$etime, fstat, xmat, failcode =1, cencode=0)


# mstate

msdata1 <- crprep("etime", fstat, trans=1, cens=0, data=tdata,
                 keep= c("sex", "age", "mspike", "id"))
fit2a <- coxph(Surv(Tstart, Tstop, status==1) ~ sex + age + mspike,
              data=msdata1, weight= weight.cens, ties='breslow')

msdata2 <- crprep("etime", fstat, trans=1, cens=0, data=tdata,
                 keep= c("sex", "age", "mspike", "id"), prec.factor=1e12)
fit2b <-  coxph(Surv(Tstart, Tstop, status==1) ~ sex + age + mspike,
              data=msdata2, weight= weight.cens, ties='breslow')

# survival
fgdata <- finegray(Surv(etime, event) ~ sex + age + mspike + id, tdata,
                   na.action=na.pass)
fit3 <- coxph(Surv(fgtime1, fgtime2, fgstatus) ~ sex + age + mspike, 
                   data=fgdata, weight=fgwt, ties='breslow')

rbind(cmprsk = fit1$coef, mstate1 = fit2a$coef, 
      mstate2= fit2b$coef, surv = fit3$coef)


# Now look at a data set with left truncation
#
fdata <- data.frame(time  =      c(1,2,3,4,4,4,5,5,6,8,8, 9,10,12),
                    init =       c(0,0,0,3,2,0,0,1,0,7,5, 0, 0, 0),
                   status=factor(c(1,2,0,1,0,0,2,1,0,0,2, 0, 1 ,0), 0:2,
                             c("cen", "type1", "type2")),      
                    x     =c(5,4,3,1,2,1,1,2,2,4,6,1,2, 0),
                    id = 1:14)
fstat <- as.numeric(fdata$status) -1
msdata3 <- crprep('time', fstat, fdata, Tstart='init', keep=c("x", "id"),
                  prec.factor=1e13) 

fg3 <- finegray(Surv(init, time, status) ~. + cluster(id), fdata)

source("~/temp/crtest.R")
source("~/Rsrc/mstate/R/create.wData.omega.R")
mtest <- crtest('time', fstat, fdata, Tstart='init', keep=c("x", "id"),
                  prec.factor=1e13) 
