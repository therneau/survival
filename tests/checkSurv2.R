library(survival)
#
# check of the Surv2 function
#
# Build a flat form of the mgus2 data set.  Mix up the data set order, to test
#  out that part of the underlying code.
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

sfit1 <- survfit(Surv2(ftime, event) ~ sex, mflat, id=id)

# now compare it to the usual way
etime <- with(mgus2, ifelse(pstat==1, ptime, futime))
estat <- with(mgus2, ifelse(pstat==1, 1, 2*death))
estat <- factor(estat, 0:2, c("censor", "progression", "death"))
sfit2 <- survfit(Surv(etime, estat) ~ sex, mgus2)

all.equal(sfit1$pstate, sfit2$pstate)

# Cox model
cfit1 <- coxph(Surv2(ftime, event) ~ sex + age, data=mflat, id=id)
cfit2 <- coxph(Surv(etime, estat)  ~ sex + age, data=mgus2, id=id)
all.equal(cfit1[c("coefficients", "var", "loglik", "score")],
          cfit2[c("coefficients", "var", "loglik", "score")])

# And using the explicit call to build a data set
sdata <- Surv2data(Surv2(ftime, event) ~ ., data=mflat, id=id)
cfit3 <- coxph(Surv2.y ~ sex + age, data=sdata, id=id)
all.equal(cfit1[c("coefficients", "var", "loglik", "score")],
          cfit3[c("coefficients", "var", "loglik", "score")])

# Create a data set with error = two events on the same day
#  A model with this data will generate an error.
temp4 <- with(m2, 
              data.frame(id=id, ftime=futime, 
                                event=ifelse(death==0, "censor", "death")))
mflat2 <- merge(temp1, rbind(temp2, temp4), all=TRUE)
mflat2$event <- factor(mflat2$event, c("censor", "prog", "death"))
stemp <- survcheck(Surv2(ftime, event) ~ sex, data=mflat2, id=id)
all.equal(stemp$duplicate$row, which(duplicated(mflat2[,c("id", "ftime")])))

# Full 3 state model.  We need to make progressions that are tied with
#  deaths be just a bit sooner.
temp2b <- with(subset(m2, pstat==1),
               data.frame(id=id, ftime= ifelse(ptime==futime & death==1, ptime-.1, ptime),
               event="progression"))
temp3b <- with(m2,
               data.frame(id=id, ftime=futime, event=ifelse(death==0, "censor", "death")))

mflat3 <- merge(temp1, rbind(temp2b, temp3b), all=TRUE)
mflat3$event <- factor(mflat3$event, c("censor", "progression", "death"))

cfit4 <- coxph(Surv2(ftime, event) ~ sex + age + mspike, mflat3, id=id)

# For a standard start-stop data set use tmerge
m3 <- tmerge(m2[,1:7], subset(m2,,c(id, futime, death)), id=id, 
                                  event= event(futime, 2*death))
m3 <- tmerge(m3, temp2b, id=id, event= event(ftime))
m3$event <- factor(m3$event, 0:2, c("censor", "progression", "death"))
cfit5 <- coxph(Surv(tstart, tstop, event) ~ sex + age + mspike, m3, id=id)

all.equal(cfit4[c("coefficients", "var", "loglik", "score")],
          cfit5[c("coefficients", "var", "loglik", "score")])

