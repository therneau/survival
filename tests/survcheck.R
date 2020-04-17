library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

# dummy data 1, simple survival
data1 <- data.frame(id=c(1,2,2,2,3,4), time=c(10, 20, 25, 30, 30, 10),
                    status=factor(c(0,0,1,1, 2,3)))
fit1 <- survcheck(Surv(time, status) ~ 1, data1, id=id)
aeq(fit1$flag, c(2,0,0,0,0))
aeq(fit1$overlap$row, 3:4)
aeq(fit1$overlap$id, 2)
aeq(fit1$transitions[1,], c(2,1,1,1))


# dummy data 2, no initial values, start stop data
# A: (0, 10, 0), (10, 20, 1), (20, 30, 2)   # no issues
# B: (0, 20, 1), (15, 24, 2), (25,26, 0)
# C: (10,13, 1), (15, 18, 0), (18,25,3)

data2 <- data.frame(id=rep(LETTERS[1:3], each=3),
                    t1 = c(0, 10, 20, 0, 15,25, 10, 15, 18),
                    t2 = c(10, 20,30, 20, 24, 26, 13, 18, 25),
                    status= factor(c(0, 1, 2, 1,2, 0, 1, 0, 3)),
                    x = c(1:5, NA, 7:9),
                    stringsAsFactors = FALSE)
fit2 <- survcheck(Surv(t1, t2, status) ~ 1, data2, id=id)

aeq(fit2$flag , c(1,2,0,0,0))
aeq(fit2$transition, rbind(c(3,0,0,0), c(0,2,1,0), c(0,0,0,1), 0))
(fit2$overlap$id == 'B')
(fit2$overlap$row ==5)
all(fit2$gap$id == c("B", "C"))
aeq(fit2$gap$row, c(6,8))

# scramble
reord <- c(9,2,1,4,3,5,6,8,7) 
tfit <- survcheck(Surv(t1, t2, status) ~ 1, data2[reord,], id=id)
all.equal(fit2[1:4], tfit[1:4])   

# let a missing value in
fit2b <- survcheck(Surv(t1, t2, status) ~ x, data2, id=id)
aeq(fit2b$flag , c(1,1,0,0,0))
aeq(fit2b$transition, rbind(c(3,0,0), c(0,2,1), 0,0))
(fit2b$overlap$id == 'B')
(fit2b$overlap$row ==5)
all(fit2b$gap$id == "C")
aeq(fit2b$gap$row, 8)

# designed to trigger all 4 error types
data3 <- data2
levels(data3$status) <- c("cens", "mgus", "recur", "death")
data3$istate <- c("entry", "entry", "recur",  "entry", "recur", "recur",
                  "entry", "recur", "recur")
fit3 <- survcheck(Surv(t1, t2, status) ~ 1, data3, id=id, istate=istate)

aeq(fit3$flag, c(1, 1, 1, 2, 0))
aeq(fit3$transitions, rbind(c(3,0,0,0), 0, c(0,2,1,1), 0))
all.equal(fit3$overlap, fit2$overlap)
all(fit3$teleport$id == c("A", "C"))
all(fit3$teleport$row == c(3,9))
all(fit3$jump$id == "C")
all(fit3$jump$row == 8)
all.equal(fit3$gap, list(row=6L, id= "B"))


