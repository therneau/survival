library(survival)
# Make sure that the old-style and new-style calls both work

# new style
vet2 <- survSplit(Surv(time, status) ~ ., data= veteran, cut=c(90, 180), 
                  episode= "tgroup", id="id")
vet2[1:7, c("id", "tstart", "time", "status", "tgroup", "age", "karno")]

# old style
vet3 <- survSplit(veteran, end='time', event='status', cut=c(90,180),
                  episode="tgroup", id="id")
all.equal(vet2, vet3)

all.equal(nrow(vet2), nrow(veteran) + sum(veteran$time >90) + 
                      sum(veteran$time > 180))


# Do a parallel computation using survSplit, and pyears/tcut.  We should get
#  the same answer.
# Break subjects up by year of entry and current time on study.  Most of
#  the deaths are within 3 months
# pyears complains (justifiably) about the obs with 0 days of fu, so we add 1
#
data1 <- jasa
data1$ayear <- as.numeric(substring(as.character(jasa$accept.dt), 1,4))
temp <- round(c(0, .25, 1,2,5)*365.25)   # years
ftime <- tcut(rep(0, nrow(jasa)), temp, 
                  labels=paste(c(0, .25, 1:2), c(.25, 1,2,5), sep='-'))
pfit <- pyears(Surv(futime +1,fustat) ~ ayear + ftime, data1,
               scale=1)

data2 <- survSplit(Surv(futime+1, fustat) ~ ., cut=temp, data=data1,
                   episode = "tgroup")

tab1 <- with(data2, tapply(fustat, list(ayear, tgroup), sum))
tab1 <- ifelse(is.na(tab1), 0, tab1)

all.equal(as.vector(tab1), as.vector(pfit$event))  # ignore dimnames

tab2 <- with(data2, tapply(tstop-tstart, list(ayear, tgroup), sum))
tab2 <- ifelse(is.na(tab2), 0, tab2)

all.equal(as.vector(tab2), as.vector(pfit$pyears))

# double check that the "data" option gives the same values
pfit2 <- pyears(Surv(futime +1,fustat) ~ ayear + ftime, data1,
               scale=1, data.frame=TRUE)$data
all.equal(pfit2$pyears, pfit$pyears[pfit$pyears >0])
all.equal(pfit2$event,  pfit$event[pfit$pyears >0])

# and that the rows of data2 have the right labels
keep <- which(pfit$pyears >0)  # these are not in the data
rname <- rownames(pfit$pyears)[row(pfit$pyears)[keep]]
all.equal(rname, as.character(pfit2$ayear))

cname <- colnames(pfit$pyears)[col(pfit$pyears)[keep]]
all.equal(cname, as.character(pfit2$ftime))

                   
