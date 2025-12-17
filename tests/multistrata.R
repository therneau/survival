library(survival)

# work out the consequences of complex strata in a multistate model
# one subject has CR and SCT on the same day, make CR a day earlier
tdata <- myeloid  # temporary working copy
tied <- with(tdata, (!is.na(crtime) & !is.na(txtime) & crtime==txtime))
tdata$crtime[tied] <- tdata$crtime[tied] -1
mdata <- tmerge(tdata[,1:4], tdata,  id=id,  death= event(futime, death),
                sct = event(txtime), cr = event(crtime), 
                relapse = event(rltime), priorrel = tdc(rltime),
                priorcr = tdc(crtime), priortx = tdc(txtime))
temp <- with(mdata, cr + 2*sct  + 4*relapse + 8*death)
table(temp)
mdata$event <- factor(temp, c(0,1,2,4,8),
                       c("none", "CR", "SCT", "relapse", "death"))
check1 <- survcheck(Surv(tstart, tstop, event) ~1, mdata, id=id)

fail <- with(mdata, c(1,1,1,2,3)[event] + 3*priorrel)
mdata$fail <- factor(c(0,1,2,0,0,3)[fail], 0:3, c("censor", "relapse", "death", 
                                  "death after relapse"))
check2 <- survcheck(Surv(tstart, tstop, fail) ~1, mdata, id=id)

temp3     <- with(mdata, c(0,2,2,0)[fail] + sct + 3*priortx)
# failure is death or relapse
mdata$tx2 <- factor(c(0,1,2,0,0,3)[1+temp3], 0:3, 
                    c("censor", "SCT", "fail w/o SCT", "fail after SCT"))
check3 <- survcheck(Surv(tstart, tstop, tx2) ~1, mdata, id=id,
                    subset= (priorrel==0))
mdata2 <- subset(mdata, priorrel==0)  # no fu after failure
hfit1 <- coxph(Surv(tstart, tstop, tx2) ~ trt + flt3 + strata(sex), 
               mdata2, id=id)
