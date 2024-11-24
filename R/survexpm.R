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
    # DOI: 10.1007/978-3-662-48971-0_15 show that in general, determining f
    # a matrix can be permuted to upper triangular is hard and give an 
    # exponential time algorithm.  The crude algorithm below can get lucky if
    # the transition matrix is sparse, which many are, but it is best if the
    # user orders states in a way that makes it easy.
    temp <- 1*(rmat >0) # 0/1 matrix
    indx <- order(colSums(temp) - rowSums(temp))
    temp <- rmat[indx, indx]  # try that ordering
    if (all(temp[lower]== 0)) indx  # it worked!
    else -1  # there is a loop in the states
}

survexpm <- function(rmat, time=1.0, setup, eps=1e-6) {
    # rmat is a transition matrix, so the diagonal elements are 0 or negative
    if (length(rmat)==1) exp(rmat*time)  #failsafe -- should never be called
    else {
        nonzero <- (diag(rmat) != 0)
        zcol  <- colSums(rmat >0) >0
        if (sum(nonzero) ==0) diag(nrow(rmat))  # expm(0 matrix) = identity
        if (sum(nonzero) ==1) {   # only one state had departures
            j <- which(nonzero)
            emat <- diag(nrow(rmat))
            temp <- exp(rmat[j,j] * time)
            emat[j,j] <- temp
            emat[j, -j] <- (1-temp)* rmat[j, -j]/sum(rmat[j,-j])
            emat
        }
        else if (sum(zcol)==1) {  # only one state has additions
            k <- which(zcol)
            emat <- diag(exp(-rmat[,k] * time))
            emat[-k,k] <- 1- diag(emat)[-k]
            emat
        }                 
        else if (missing(setup) || setup[1] < 0 ||
                 any(diff(sort(diag(rmat)))< eps)) pade(rmat*time)
        else {
            if (setup[1]==0) .Call(Ccdecomp, rmat, time)$P
            else {
                temp <- .Call(Ccdecomp, rmat[setup, setup], time)
                temp$P[order(setup2), order(setup2)]
            }
        }
    }
}

expmderiv<- function(rmat, time=1.0, dR, setup, eps=1e-8) {
    if (length(rmat)==1) { # should not happen, not a multi-state model
        stop("deriv function called with 1x1 matrix")
    }
    nr <- nrow(rmat)
    nonzero <- (diag(rmat) != 0)
    if (sum(nonzero) ==0) { # expm(0 matrix) = identity
        return(list(P= diag(nr), dmat=array(0, dim=c(nr, nr, length(dR)))))
    }
               
    if (sum(nonzero) ==1) {   # only one state had departures
        stop("not yet filled in")
    }
    zcol <- colSums(rmat>0) >0  # states with at least one inward transition
    if (sum(zcol)==1) {
        stop("not yet filled in")
    }

    if (missing(setup) || setup[1] <0 || any(diff(sort(diag(rmat)))< eps)) 
        return (pade(rmat*time, dR*time))

    # The cholesky case
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
    
    # undo the reordering, if needed
    if (setup[1] >0) {
        indx <- order(setup)
        list(P = dlist$P[indx, indx], dmat = dmat[indx, indx, ])
    } else list(P= dlist$P, dmat=dmat)
}
