# An example from the epifit package that blew up coxph

id <- 1:8
group <- c(0, 0, 0, 0, 1, 1, 1, 1)
time <- c(4, 5, 7, 9, 6, 10, 11, 12)
event <- c(1, 1, 0, 1, 1, 1, 1, 0)
dat5 <- data.frame(id, group, time, event)

cov0 <- data.frame(id=id, time=0, value=0*group)
cov4 <- data.frame(id=id, time=3.9, value=4*group)
cov5 <- data.frame(id=id, time=4.9, value=5*group)
cov6 <- data.frame(id=id, time=5.9, value=6*group)
cov9 <- data.frame(id=id, time=8.9, value=9*group)
cov10 <- data.frame(id=id, time=9.9, value=10*group)
cov11 <- data.frame(id=id, time=10.9, value=11*group)
tdata <- data.frame(id=id, group=group)
tdata <- tmerge(tdata, dat5, id, status=event(time, event))
tdata <- tmerge(tdata, cov0, id, t_g=tdc(time, value))
tdata <- tmerge(tdata, cov4, id, t_g=tdc(time, value))
tdata <- tmerge(tdata, cov5, id, t_g=tdc(time, value))
tdata <- tmerge(tdata, cov6, id, t_g=tdc(time, value))
tdata <- tmerge(tdata, cov9, id, t_g=tdc(time, value))
tdata <- tmerge(tdata, cov10, id, t_g=tdc(time, value))
tdata <- tmerge(tdata, cov11, id, t_g=tdc(time, value))

fail <- coxph(Surv(tstart, tstop, status) ~ group + t_g, data=tdata)
