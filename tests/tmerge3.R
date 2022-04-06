# Understanding edge cases
library(survival)

#
# this is from a user report of a problem with cumevents.  When there is
#  a row merged in that is a censor, don't mark it as a cumevent.
#
base <- data.frame(
  id = 1:2, tstart = c(0, 0), tstop  = c(10, 10), got_flu = c(0, 0), 
  has_flu = factor(c("no", "no"), levels = c("no", "yes")))
base <- tmerge(base, base, id = id, got_flu = event(tstop, got_flu))

# add time-varying covariates
vars <- data.frame(id = c(1, rep(2, 5)), time = c(0, (0:4) * 2), x = rnorm(6))
base <- tmerge(base, vars, id = id, x = tdc(time, x))

# add cumevents, using a covariate
events <- data.frame(
  id = c(2, 2, 2), 
  # notice the zero -- the second row should not add an event
  got_flu = c(1,0,2), 
  has_flu = c("yes", "no", "yes"),
  time = c(3, 5, 8))
b2 <- tmerge(base, events, id = id, got_flu = cumevent(time, got_flu),
               has_flu = tdc(time, has_flu))

all.equal(b2$got_flu, c(0,0,1,0,0,0,3,0))


# Tied times in the merger data set
# for all of them missings are essentially ignored
# last obs wins for tdc and event
tiedat <- data.frame(id=c(1, 1, 1, 2,2,2), time=c(3,4, 4, 3, 5, 5), 
                     x=c(1, NA,0, 2,3,4))
b3 <- tmerge(base, tiedat, id=id, x1= tdc(time, x), x2=cumtdc(time, x),
             x3= event(time, x), x4 = cumevent(time, x))
all.equal(b3$x1, c(NA, 1, 0, NA, NA, 2,2, 4,4,4))
all.equal(b3$x2, c(NA, 1, 1, NA, NA, 2,2, 9,9,9))
all.equal(b3$x3, c(1,0,0,0,2,0,4,0,0,0))
all.equal(b3$x4, c(1,0,0,0,2,0,9,0,0,0))

# Multiple overlapping time windows in the first step.
#  Should generate an error message
test <- tryCatch(
            {tmerge(pbcseq[, c("id", "trt", "age", "sex")], pbcseq, id,
               death = event(futime, status==2))},
            error= function(cond) {
                if (grepl("duplicate identifiers", cond)) 
                    cat("successful tmerge error test\n")
            }
)

 
