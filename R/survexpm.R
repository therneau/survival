# Automatically generated from the noweb directory
survexpmsetup <- function(rmat) {
    # check the validity of the transition matrix, and determine if it
    #  is acyclic, i.e., can be reordered into an upper triangular matrix.
    if (!is.matrix(rmat) || nrow(rmat) != ncol(rmat) || any(diag(rmat) > 0) ||
        any(rmat[row(rmat) != col(rmat)] < 0))
        stop ("input is not a transition matrix")
    if (!is.logical(all.equal(rowSums(rmat), rep(0, ncol(rmat)))))
        stop("input is not a transition matrix")
    nc <- ncol(rmat)
    lower <- row(rmat) > col(rmat)
    if (all(rmat[lower] ==0))  return(0)  # already in order
    
    # score each state by (number of states it follows) - (number it precedes)
    temp <- 1*(rmat >0) # 0/1 matrix
    indx <- order(colSums(temp) - rowSums(temp))
    temp <- rmat[indx, indx]  # try that ordering
    if (all(temp[lower]== 0)) indx  # it worked!
    else -1  # there is a loop in the states
}
survexpm <- function(rmat, time=1.0, setup, eps=1e-6) {
    # rmat is a transition matrix, so the diagonal elements are 0 or negative
    if (length(rmat)==1) exp(rmat[1]*time)  #failsafe -- should never be called
    else {
        nonzero <- (diag(rmat) != 0)
        if (sum(nonzero ==0)) diag(nrow(rmat))  # expm(0 matrix) = identity
        if (sum(nonzero) ==1) {
            j <- which(nonzero)
            emat <- diag(nrow(rmat))
            temp <- exp(rmat[j,j] * time)
            emat[j,j] <- temp
            emat[j, -j] <- (1-temp)* rmat[j, -j]/sum(rmat[j,-j])
            emat
        }
        else if (missing(setup) || setup[1] < 0 ||
                 any(diff(sort(diag(rmat)))< eps)) pade(rmat*time)
        else {
            if (setup[1]==0) .Call(Ccdecomp, rmat, time)$P
            else {
                temp <- rmat
                temp[setup, setup] <- .Call(Ccdecomp, rmat[setup, setup], time)
                temp$P
            }
        }
    }
}
derivative <- function(rmat, time, dR, setup, eps=1e-8) {
    if (missing(setup) || setup[1] <0 || any(diff(sort(diag(rmat)))< eps)) 
        return (pade(rmat*time, dR*time))

    if (setup==0) dlist <- .Call(Ccdecomp, rmat, time)
    else dlist <- .Call(Ccdecomp, rmat[setup, setup], time)
    ncoef <- dim(dR)[3]
    nstate <- nrow(rmat)
    
    dmat <- array(0.0, dim=c(nstate, nstate, ncoef))
    vtemp <- outer(dlist$d, dlist$d,
                   function(a, b) {
                       ifelse(abs(a-b)< eps, time* exp(time* (a+b)/2),
                         (exp(a*time) - exp(b*time))/(a-b))})

    # two transitions can share a coef, but only for the same X variable
    for (i in 1:ncoef) {
        G <- dlist$Ainv %*% dR[,,i] %*% dlist$A
        V <- G*vtemp
        dmat[,,i] <- dlist$A %*% V %*% dlist$Ainv
    }
    dlist$dmat <- dmat
    
    # undo the reordering, if needed
    if (setup[1] >0) {
        indx <- order(setup)
        dlist <- list(P = dlist$P[indx, indx],
                      dmat = apply(dmat,1:2, function(x) x[indx, indx]))
    }
                      
    dlist
}
