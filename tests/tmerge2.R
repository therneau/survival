library(survival)

# This test is based on a user report that a 0/1 variable would not reset
#  to zero.  It turned out to be a bug when data2 was not sorted

baseline <- data.frame(idd=1:5,  futime=c(20, 30, 40, 30, 20),
                       status= c(0, 1, 0, 1, 0))
tests <- data.frame(idd = c(2,3,3,3,4,4,5),
                    date = c(25, -1, 15, 23, 17, 19, 14),
                    onoff= c( 1, 1, 0, 1, 1, 0, 1))
tests <- tests[c(7,2,6,3,4,1,5),]  #scramble data2

mydata <- tmerge(baseline, baseline, id=idd, death=event(futime, status))
mydata <- tmerge(mydata, tests, id=idd, ondrug=tdc(date, onoff))

all.equal(mydata$ondrug, c(NA, NA,1, 1,0,1, NA, 1,0, NA, 1))


# Check out addition of a factor, character, and logical
tests$ff <- factor(tests$onoff, 0:1, letters[4:5])
tests$fchar <- as.character(tests$ff)
tests$logic <- as.logical(tests$onoff)

mydata <- tmerge(mydata, tests, id=idd, fgrp= tdc(date, ff),
                 chgrp = tdc(date, fchar), 
                 options=list(tdcstart="new"))
all.equal(mydata$fgrp, 
          factor(c(3,3,2,2,1,2,3,2,1,3,2), labels=c("d", "e", "new")))
all.equal(mydata$chgrp, 
          c("d", "e", "new")[c(3,3,2,2,1,2,3,2,1,3,2)])

mydat2  <-  tmerge(mydata, tests, id=idd, 
                 logic1 = tdc(date, logic), logic2= event(date, logic))
all.equal(mydat2$logic1, c(FALSE, TRUE, NA)[as.numeric(mydat2$fgrp)])
all.equal(mydat2$logic2, c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE,
                           FALSE, FALSE, TRUE, FALSE))

# Multiple chained calls.  
newcgd <- tmerge(data1=cgd0[, 1:13], data2=cgd0, id=id, tstop=futime)
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime1)) 
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime2)) 
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime3)) 
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime4)) 
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime5)) 
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime6)) 
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime7)) 
newcgd <- tmerge(newcgd, newcgd, id, enum=cumtdc(tstart))
all.equal(dim(newcgd), c(203,17))
all.equal(as.vector(table(newcgd$infect)), c(127, 76))

tcount <- attr(newcgd, "tcount")
all(tcount[,1:3] ==0)  # no early, late, or gap

# table with number of subjects who have etime1 < futime (row 1)
#  and  etime1==futime (row 2)
# the table command ignores the missings
temp <- subset(cgd0, select=etime1:etime7)
counts <- sapply(temp, function(x) 
    as.vector(table(factor(x>= cgd0$futime, c(FALSE, TRUE)))))

all(tcount[1:7, c("within", "trailing")] == t(counts))


#
# Merging with a date as the time variable.  In this case tstart/tstop are required
#  A default start of 0 has no meaning
#
base2 <- baseline
base2$date1 <- as.Date("1953-03-10")   # everyone enrolled that day
base2$date2 <- as.Date("1953-03-10") + base2$futime
base2$futime <- NULL
test2 <- tests
test2$date <- as.Date("1953-03-10") + test2$date

mydata2 <- tmerge(base2, base2, id=idd, death=event(date2, status), 
                 tstart = date1, tstop= date2,
                 options=list(tstartname="date1", tstopname="date2"))
mydata2 <- tmerge(mydata2, test2, id=idd, ondrug=tdc(date, onoff))
all.equal(mydata$ondrug, c(NA, NA,1, 1,0,1, NA, 1,0, NA, 1))
