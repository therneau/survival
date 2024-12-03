#
# Test the survexpm code, which is not exported so we need :::
#
library(survival)
survexpm <- survival:::survexp

nstate <- 10
set.seed(1953)
rmat <- matrix(runif(nstate^2,0, .3), nstate, nstate)
rmat[col(rmat) <= row(rmat)] <- 0
diag(rmat) <- -rowSums(rmat)

setup <- survival:::survexpmsetup(rmat) #should= 0, the matrix is in normal form
setup ==0

indx <- sample(1:10, 10, replace=FALSE)
rmat2 <- rmat[indx, indx]
setup2 <- survival:::survexpmsetup(rmat2)
all.equal(rmat2[setup2, setup2], rmat)  # it found the right order

e1 <- as.matrix(Matrix::expm(rmat))
e2 <- survexpm(rmat, 1, 0)  # diagonal method
e3 <- survexpm(rmat, 1, -1) # Pade method
all.equal(e1, e2)
all.equal(e1, e3)

e4 <- survexpm(rmat2, 1, setup2)
all.equal(e4,as.matrix(Matrix::expm(rmat2)))
all.equal(e4[setup2, setup2], e1)

if (FALSE) {  
    # only do this offline. It won't give the same values on two
    #  machines, so would trigger a difference in R CMD check
    n <- 1e5
    t1 <- system.time({for (i in 1:n) z <- survexpm(rmat, 1, 0)})
    t2 <- system.time({for (i in 1:n) z <- survexpm(rmat, 1, -1)})
    t3 <- system.time({for (i in 1:n) z <- Matrix::expm(rmat)})

    # For nstate=10 the cholesky approach just barely beats Pade, for 
    #   nstate =6 Pade won.  I was surprised.  
    # Expm is fastest, perhaps due to sparse matrix routines?
    #
}
# Now for derivatives, all nstate* (nstate-1)/2 of them
# The dr array is nstate x nstate x (number of derivatives I want).
nderiv <- nstate * (nstate-1)/2
dr <- array(0, dim=c(nstate,nstate, nderiv))
k <- 1
for (i in 1:(nstate-1)) {
    for (j in (i+1):nstate) {
        dr[i,j,k] <- 1
        k <- k+1
    }
}

test1 <- expmderiv(rmat, 1, dr, 0)
test2 <- expmderiv(rmat, 1, dr, -1)
all.equal(test1, test2)

if (FALSE) {
    n <- 2e4
    t3 <- system.time({for (i in 1:n) z <- expmderiv(rmat, 1, dr, 0)})
    t4 <- system.time({for (i in 1:n) z <- expmderiv(rmat, 1, dr, -1)})

    # The diagonal algorithm is substantially faster than Pade for 10 states (7x)
}
eps <- 1e-6
check <- logical(nderiv) # derivatives by hand
for (i in 1:nderiv) {
    d2 <- (as.matrix(expm(rmat+ eps*dr[,,i])) - e1)/eps
    check[i] <- isTRUE(all.equal(d2, test1$dmat[,,i]))
}
all(check)
    
