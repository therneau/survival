library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

# Using id as a character string is no longer supported (all mention of it has
#  be removed from documentation), but we do allow it, silently, for backward
#  compatability
vet2 <- survSplit(Surv(time, status) ~ ., data= veteran, cut=c(90, 180), 
                  episode= "tgroup", id="id")
vet2[1:7, c("id", "tstart", "time", "status", "tgroup", "age", "karno")]
all(vet2$id[1:7] == c(1,2,2,2,3,3,3))
all(vet2$tstart[1:7] == c(0,0, 90, 180, 0, 90, 180))
all(vet2$tstop[1:7] ==  c(72, 90, 180, 411, 90, 180, 288))

if (FALSE) { # we no longer support this style
# old style call
vet3 <- survSplit(veteran, end='time', event='status', cut=c(90,180),
                  episode="tgroup", id="id")
all.equal(vet2, vet3)

all.equal(nrow(vet2), nrow(veteran) + sum(veteran$time >90) + 
                      sum(veteran$time > 180))
}

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
aeq(tab1, pfit$event)

tab2 <- with(data2, tapply(tstop-tstart, list(ayear, tgroup), sum))
tab2 <- ifelse(is.na(tab2), 0, tab2)
aeq(tab2, pfit$pyears)

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

# make sure that ties between the data and the cutpoints are okay
# a roundoff issue caused duplicate times, at one time
tdata <- data.frame(id=c(1,1,1,1),
                    t1=c(0, 1.2, 2.8, 4.0),
                    t2=c(1.2, 2.8, 4.0, 5.9),
                    status=c(0,0,0,0))
new <- survSplit(Surv(t1,t2, status) ~., data=tdata, cut=c(1,4))
all(!duplicated(round(new$t1, 5)))

#
# Test using timeline data
tdata <- data.frame(id=c(1,2,2,2,3,3,3,4,4), time=c(5, 3, 8, 11, 0, 5,9, 2,7),
                  state=factor(c(2,2,1,3,3,1,2, 2,3), 1:3, c("none", 'A', 'B')),
                  x=1:9, y= 10:2)

target <- data.frame(id=   c(1, 2,2,2,2, 2, 3,3,3,3,3, 4,4,4), 
                     time =c(5, 3,6,8,9,11, 0,2,5,6,9, 2,6,7),
                     state=c(2, 2,1,1,1, 3, 3,1,1,1,2, 2,1,3),
                     x    =c(1, 2,2,3,3, 4, 5,5,6,6,7, 8,8,9),
                     y    =c(10,9,9,8,8, 7, 6,6,5,5,4, 3,3,2))
target$state <- factor(target$state, 1:3, c("censor", 'A', 'B'))

s1 <- survSplit(Surv2(time, state)~., tdata, id=id, cut= c(2,6,9))
all.equal(s1, target)
