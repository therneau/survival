# Some checks of the internal survexpm routine
# It is not exported so I need the ::: style call
#
library(survival)
library(Matrix)

q1 <- matrix(0, 5, 5)  # the simple model of the NAFLD data
q1[1,2] <- q1[2,3] <- q1[3,4] <- 1
q1[1:4,5] <- 1

set.seed(1960)
rmat <- q1
rmat[rmat>0] <- exp(runif(7, -1, 1))
diag(rmat) <- -rowSums(rmat)

s1 <- survival:::survexpmsetup(rmat)
e1 <- survival:::survexpm(rmat, 2, s1)  # use the decomposition method
e2 <- survival:::survexpm(rmat, 2)      # use my Pade
e3 <- as.matrix(expm(2*rmat)) # use Matrix

all.equal(e1, e2)
all.equal(e1, e3)

#
# Compute derivatives
#
dR <- array(0, dim=c(5,5,7))  # the setup is a bit of a nuisance.
indx1 <- row(q1)[q1>0]
indx2 <- col(q1)[q1>0]
for (i in 1:7) dR[indx1[i], indx2[i],i] <- 1  # what derivatives do I want

d1 <- survival:::expmderiv(rmat, 2, dR, s1)
all.equal(d1$P, e1)

eps <- 1e-1
# This portion not yet correct
for (i in 1:7) {
    rtemp <- rmat
    rtemp[indx1[i], indx2[i]] <- rtemp[indx1[i], indx2[i]] + eps
    rtemp[indx1[i], indx1[i]] <- rtemp[indx1[i], indx1[i]] - eps
    ptemp <- as.matrix(expm(2*rtemp))
    delta <- (ptemp- e1)/eps
}
