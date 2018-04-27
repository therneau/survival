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


# Check out addition of a factor
tests$ff <- factor(tests$onoff, 0:1, letters[4:5])
mydata <- tmerge(mydata, tests, id=idd, fgrp= tdc(date, ff),
                 options=list(tdcstart="new"))

all.equal(mydata$fgrp, 
          factor(c(3,3,2,2,1,2,3,2,1,3,2), labels=c("d", "e", "new")))


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

