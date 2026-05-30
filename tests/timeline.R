library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

#
# Follow on to the checkSurv2.R tests
#  Multi-state data set with repeated states and covariates
#
data1 <- subset(pbcseq,, c(id, day, age, bili, albumin, edema, ast))
data1$bstat <- as.numeric( cut(data1$bili, c(0, 1, 4, 100)))
temp <- subset(pbcseq, !duplicated(id), c(id, futime, status))
names(temp) <- c("id", "day", "bstat")
temp$bstat <- 4*(temp$bstat==2)
pdata <- merge(data1, temp, all=TRUE)
pdata$bstat <- factor(pdata$bstat, 0:4, 
                      c("censor", "normal", "1-4", "4+", "death"))

# Throw in some missings: id 1 is missing edema, id 2 missing even visit
#  albumins, id 3 missing the first ast, bstat2 is missing on row 2.
# Models with edema should throw out subject 1, those with albumin will have
#  albumin filled in, and those with ast will have delayed entry of subject 3.
# A model using bstat2 should fail survcheck since id 1 will have a "hole"
pdata$edema[pdata$id==1] <- NA
indx <- which(pdata$id ==2)
indx <- indx[indx%%2== 1]
pdata$albumin[indx] <- NA
pdata$ast[which(pdata$id==3)[1]] <- NA
pdata$bstat2 <- pdata$bstat
pdata$bstat2[2] <- NA

subset(pdata, id<4)  # the printout looks good
cdata <- fromtimeline(Surv(day, bstat) ~ ., pdata, id=id)
subset(cdata, id<4)
# in the above, bstat2 also gets the lvcf treatment -- it's just a covariate

fit1a <- coxph(list(Surv2(day, bstat) ~ 1,
                    (1:3):4 ~ edema + albumin + ast), pdata, id=id)

fit1b <- coxph(list(Surv(day1, day2, bstat) ~ 1,
                    (1:3):4 ~ edema + albumin + ast), 
               cdata, id=id, istate=istate)
all.equal(fit1a[c("coefficients", "var", "loglik", "score")],
          fit1b[c("coefficients", "var", "loglik", "score")])

# reconstruct the model frame
fit1c <- coxph(list(Surv2(day, bstat) ~ 1,
                    (1:3):4 ~ edema + albumin + ast), pdata, id=id,
               model=TRUE)
m1 <- model.frame(fit1a)
m2 <- model.frame(fit1b)
m3 <- model.frame(fit1c)
all.equal(m1,m3)
# m2 will have a different name for the first column, and that first col
#  has a different "inputAttributes" attribute
all.equal(m2[,-1], m3[,-1])
aeq(m2[,1], m3[,1])

dummy <- expand.grid(edema=0:1, albumin=3, ast=c(75, 150))
s1 <- survfit(fit1a, newdata=dummy)
s2 <- survfit(fit1b, newdata=dummy)
all.equal(unclass(s1)[c("time", "pstate", "n.risk", "cumhaz")],
          unclass(s2)[c("time", "pstate", "n.risk", "cumhaz")])
