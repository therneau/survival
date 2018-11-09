#
# This routine creates a stacked data set
# The key input is the cmap matrix, which has one row for each column 
#  of X (plus a first row for the intercept of a Cox model which turns into
#  strata) and one column for each transition.
# Input data is X, Y, strata, and initial state (integer).
#
# For each transition the expanded data has a set of rows, all those whose
#  initial state makes them eligible for the transition.  
#
stacker <- function(cmap, istate, x, y, strata, states) {
    from.state <- as.numeric(sub(":.*$", "", colnames(cmap)))
    to.state   <- as.numeric(sub("^.*:", "", colnames(cmap)))

    # just in case cmap has columns I don't need (I don't think this can
    #  happen
    check <- match(from.state, istate, nomatch==0)
    if (any(check===0)){
        warning("extra column in cmap")  # debugging line
        cmap <- cmap[,check>0]
        from.state <- from.state[check>0]
        to.state <- to.state[check>0]
    }
 
    endpoint <- match(attr(y, "states"), "states", nomatch=0)
    endpoint <- endpoint[ 1 + y[,ncol(y)])
    n2 <- sum(table(from.state)* table(to.state))  # number of rows in new data
    newX <- matrix(0L, nrow=n2, ncol=max(cmap2))
    k <- 0
    rindex <- integer(n2)   # original row for each new row of data
    newstat <- integer(n2)  # new status
    for (i in 1:ncol(cmap)) {
        subject <- which(istate= from.state[i])
        nr <-k + seq(along=subject)  # rows in new matrix
        nc <- cmap[-1,i]              # cols in new matrix
        newX[j, nc[nc>0]] <- X[subject, nc>0]
        rindex[nr] <- subject
        newstat[nr] <- ifelse(endpoint[subject] == to.state[i], 1, 0)
    }
 
    newY <- cbind(y[rindex,-ncol(y)], newstat)

    # which transition each row represents
    transition <- rep(1:ncol(cmap), (table(istate))[from.state])
    
    # new strata
    if (is.null(strata)) newstrat <- cmap[1,transition]
    else {
        # if the old strata is 1, 2, 3 the new ones will be 11-13 for transition
        temp <-10^ ceiling(log(max(strata), 10))
        newstrat <- strata[rindex] + cmap[1,transition]*temp
   } 
