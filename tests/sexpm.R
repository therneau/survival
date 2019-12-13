# Test the sexpm functions that are used for fast matrix exponentials
#
library(Matrix)
library(survival)

nfun <- length(sexpm)
times <- c(.1, 1, 5)
test1 <- matrix(FALSE, nfun, length(times), 
                dimnames=list(names(sexpm), paste0("t=", times)))
test2 <- rep(TRUE, nfun)
eps <- 1e-8
for (i in 1:nfun) {
    j <- sexpm[[i]]
    tmat <- matrix(0, j$nstate, j$nstate)
    tmat[j$nonzero] <- runif(length(j$nonzero), .1, 4)
    diag(tmat) <- -rowSums(tmat)
    
    dtemp <- logical(length(j$nonzero))
    for (k in 1:length(times)) {
        m1 <- j$mexp(tmat, times[k])
        m2 <- expm(tmat*times[k])
        test1[i,k] <- isTRUE(all.equal(m1, as.matrix(m2), check.attributes=FALSE))
        # now derivatives
        d1 <- j$deriv(tmat, times[k])
        for (nz in 1:length(j$nonzero)) {
            temp <- tmat
            temp[j$nonzero[nz]] <- temp[j$nonzero[nz]] + eps
            diag(temp) <- diag(temp)- rowSums(temp)
            d2 <- (expm(temp* times[k]) - m2)/eps
browser()
            test2[i] <- test2[i] & isTRUE(all.equal(d1[,,nz], d2, tol=eps/2))
        }   
    }
}   

all(test1)  # should all be TRUE

# Now the derivatives
eps <- 1e-8
test2 <- logical(nfun)
for (i in 1:nfun) {
    j <- sexpm[[i]]
    tmat <- matrix(0, j$nstate, j$nstate)
    tmat[j$nonzero] <- runif(length(j$nonzero), .1, 4)
    diag(tmat) <- -rowSums(tmat)

