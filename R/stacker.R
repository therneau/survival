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
stacker <- function(cmap, istate, X, Y, strata, states) {
    from.state <- as.numeric(sub(":.*$", "", colnames(cmap)))
    to.state   <- as.numeric(sub("^.*:", "", colnames(cmap)))

    # just in case cmap has columns I don't need (I don't think this can
    #  happen
    check <- match(from.state, istate, nomatch=0)
    if (any(check==0)){
        # I think that this is impossible
        browser()
        warning("extra column in cmap")  # debugging line
        cmap <- cmap[,check>0]
        from.state <- from.state[check>0]
        to.state <- to.state[check>0]
    }
 
    endpoint <- c(0, match(attr(Y, "states"), states))
    endpoint <- endpoint[ 1 + Y[,ncol(Y)]]
    n2 <- sum(table(from.state)* table(istate))  # number of rows in new data
    newX <- matrix(0L, nrow=n2, ncol=max(cmap))
    k <- 0
    rindex <- integer(n2)   # original row for each new row of data
    newstat <- integer(n2)  # new status
    for (i in 1:ncol(cmap)) {
        subject <- which(istate==from.state[i])
        nr <-k + seq(along=subject)  # rows in new matrix
        nc <- cmap[-1,i]              # cols in new matrix
        newX[nr, nc[nc>0]] <- X[subject, nc>0]
        rindex[nr] <- subject
        newstat[nr] <- ifelse(endpoint[subject] == to.state[i], 1, 0)
        k <- max(nr)
    }

    if (ncol(Y) ==2) newY <- Surv(Y[rindex,1], newstat)
    else newY <- Surv(Y[rindex,1], Y[rindex,2], newstat)

    # which transition each row represents
    transition <- rep(1:ncol(cmap), (table(istate))[from.state])
    
    # new strata
    if (is.null(strata)) newstrat <- cmap[1,transition]
    else {
        # if the old strata is 1, 2, 3 the new ones will be 11-13 for transition
        temp <-10^ ceiling(log(max(strata), 10))
        newstrat <- strata[rindex] + cmap[1,transition]*temp
    } 

    # give variable names to the new data
    vname <- rep("", ncol(X))
    ctemp <- cmap[-1,]
    vname[ctemp[ctemp>0]] <- colnames(X)[row(ctemp)[ctemp>0]]
    colnames(newX) <- vname

    list(X=newX, Y=newY, strata=newstrat, transition=transition, rindex=rindex)
}
