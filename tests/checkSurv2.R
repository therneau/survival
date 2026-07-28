library(survival)
#
# check of the Surv2 function
#
# Build a flat form of the mgus2 data set.  Mix up the data set order, to test
#  out that part of the underlying code.
# I had a trial run at autodetecting timeline data sets, i.e., one could
#  use Surv instead of Surv2. It worked well for these examples but we 
#  decided not to go that route: the sfit3, cfit3, cfit6 tests have been left
#  in for possible future reference, but are commented out. 
set.seed(1953)
m2 <- mgus2[sample(1:nrow(mgus2), nrow(mgus2),replace=FALSE),]

temp1 <- data.frame(m2[,1:7], ftime=0)
temp2 <- with(subset(m2, pstat==1), 
              data.frame(id=id, ftime=ptime, event="progression"))

# competing risks: use only the first of death and progression
temp3 <- with(subset(m2, pstat==0),
              data.frame(id=id, ftime=futime, 
                                event=ifelse(death==0, "censor", "death")))
mflat <- merge(temp1, rbind(temp2, temp3), all=TRUE)
mflat$event <- factor(mflat$event, c("censor", "progression", "death"))
mflatc <- fromtimeline(Surv2(ftime, event) ~., mflat, id=id)

# now compare it to the usual way
etime <- with(mgus2, ifelse(pstat==1, ptime, futime))
estat <- with(mgus2, ifelse(pstat==1, 1, 2*death))
estat <- factor(estat, 0:2, c("censor", "progression", "death"))

# The sfit3 lines are from the brief time that I thought I could dispense with
#  Surv2, (time, status) data with multiple rows per subject implied it. But
#  Beth/Cindy pointed out the diabetic retinopathy data
sfit1 <- survfit(Surv(etime, estat) ~ sex, mgus2)  # original way
sfit2 <- survfit(Surv2(ftime, event) ~ sex, mflat, id=id) # timeline with Surv2
#sfit3 <- survfit(Surv(ftime, event) ~ sex, mflat, id=id)  # timeline with Surv
sfit4 <- survfit(Surv(ftime, event) ~ sex, mflatc, id=id)# converted timeline

all.equal(sfit1$pstate, sfit2$pstate)
#all.equal(sfit1$pstate, sfit3$pstate)
all.equal(sfit1$pstate, sfit4$pstate)

# Cox model
cfit1 <- coxph(Surv(etime, estat)  ~ sex + age, data=mgus2, id=id)
cfit2 <- coxph(Surv2(ftime, event) ~ sex + age, data=mflat, id=id)
#cfit3 <- coxph(Surv(ftime, event)  ~ sex + age, data=mflat, id=id)
cfit4 <- coxph(Surv(ftime, event)  ~ sex + age, data=mflatc, id=id)
all.equal(cfit1[c("coefficients", "var", "loglik", "score")],
          cfit2[c("coefficients", "var", "loglik", "score")])
#all.equal(cfit1[c("coefficients", "var", "loglik", "score")],
#          cfit3[c("coefficients", "var", "loglik", "score")])
all.equal(cfit1[c("coefficients", "var", "loglik", "score")],
          cfit4[c("coefficients", "var", "loglik", "score")])

# Full 3 state model.  We need to make progressions that are tied with
#  deaths be just a bit sooner.
temp2b <- with(subset(m2, pstat==1),
               data.frame(id=id, ftime= ifelse(ptime==futime & death==1, 
                                               ptime-.1, ptime),
               event="progression"))
temp3b <- with(m2,
               data.frame(id=id, ftime=futime, event=ifelse(death==0,   
                                                            "censor", "death")))
mflat3 <- merge(temp1, rbind(temp2b, temp3b), all=TRUE)
mflat3$event <- factor(mflat3$event, c("censor", "progression", "death"))

# For a standard start-stop data set use tmerge
m3 <- tmerge(m2[,1:7], subset(m2,,c(id, futime, death)), id=id, 
                                  event= event(futime, 2*death))
m3 <- tmerge(m3, temp2b, id=id, event= event(ftime))
m3$event <- factor(m3$event, 0:2, c("censor", "progression", "death"))

cfit4 <- coxph(Surv(tstart, tstop, event) ~ sex + age + mspike, m3, id=id)
cfit5 <- coxph(Surv2(ftime, event) ~ sex + age + mspike, mflat3, id=id)
#cfit6 <- coxph(Surv(ftime, event) ~ sex + age + mspike, mflat3, id=id)
mflat3c <- fromtimeline(Surv2(ftime, event) ~ ., mflat3, id=id)
cfit7 <- coxph(Surv(ftime1, ftime2, event) ~ sex + age + mspike, mflat3c, id=id)

all.equal(cfit4[c("coefficients", "var", "loglik", "score")],
          cfit5[c("coefficients", "var", "loglik", "score")])
#all.equal(cfit4[c("coefficients", "var", "loglik", "score")],
#          cfit6[c("coefficients", "var", "loglik", "score")])
all.equal(cfit4[c("coefficients", "var", "loglik", "score")],
          cfit7[c("coefficients", "var", "loglik", "score")])

