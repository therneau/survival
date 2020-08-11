library(survival)
#
# Tests of the residuals.survfit function
#
# The influence argument of survfit returns all the residuals at every time
#  point, but for large data sets the result will be huge.  This function uses
#  a different algorithm which should be faster when the number of time
#  points being reported out is small. 
# For testing, we can compare it to the results from the influence arg, which
#  is tested in the mstate.R file.
#

aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))
tdata <- data.frame(id= c(1, 1, 1,  2,  3,  4, 4, 4,  5,  5),
                    t1= c(0, 4, 9,  1,  2,  0, 2, 8,  1,  3),
                    t2= c(4, 9, 10, 5,  9,  2, 8, 9,  3, 11),
                    st= c(1, 2,  1, 2,  3,  1, 3, 0,  3,  0),
                    i0= c(1, 2,  3, 2,  1,  1, 2, 4,  3,  4),
                    wt= 1:10)

tdata$st <- factor(tdata$st, c(0:3),
                    labels=c("censor", "a", "b", "c"))
tdata$i0 <- factor(tdata$i0, 1:4,
                    labels=c("entry", "a", "b", "c"))  

tfun <- function(data=tdata) {
    reorder <- c(10, 9, 1, 2, 5, 4, 3, 7, 8, 6)
    new <- data[reorder,]
    new
}
mtest2 <- tfun(tdata)  # scrambled version


mfit1 <- survfit(Surv(t1, t2, st) ~ 1, tdata, id=id, istate=i0,
                 influence=1)

test1 <- resid(mfit1, time=c(3, 7))
aeq(aperm(test1, c(1,3,2)), mfit1$influence.pstate[,c(3,5),])


