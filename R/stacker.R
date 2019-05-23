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
        warning("extra column in cmap")  # debugging line
        browser()
        cmap <- cmap[,check>0]
        from.state <- from.state[check>0]
        to.state <- to.state[check>0]
    }
 
    endpoint <- c(0, match(attr(Y, "states"), states))
    endpoint <- endpoint[ 1 + Y[,ncol(Y)]]  # endpoint of each row, 0=censor

    # Usually each transition is a separate stratum, but the user can set it
    #  up otherwise.  The first row of cmat gives the strata for each change
    ustrata <- unique(cmap[1,])
    nstrat <- length(ustrata)
    n.perstrat <- integer(nstrat)
    for (i in 1:nstrat) {
        itemp <- unique(from.state[cmap[1,] == ustrata[i]])
        n.perstrat[i] <- sum(istate %in% itemp)
    }
    
    # The constructed X matrix has a block or rows for each ustrata level
    n2 <- sum(n.perstrat)  # number of rows in new data
    newX <- matrix(0L, nrow=n2, ncol=max(cmap[-1,]))
    k <- 0
    rindex <- integer(n2)   # original row for each new row of data
    newstat <- integer(n2)  # new status
    for (i in 1:nstrat) {
        whichcol <- (cmap[1,] == ustrata[i])  # cols of cmap to look at
        subject <- which(istate %in% from.state[whichcol]) #participants in block
        nr <-k + seq(along=subject)  # rows in new matrix
        nc <- cmap[-1,i]             # variables in new matrix
        newX[nr, nc[nc>0]] <- X[subject, nc>0]
        rindex[nr] <- subject
        
        event.that.counts <- (endpoint[subject] %in% to.state[whichcol])
        newstat[nr] <- ifelse(event.that.counts, 1L, 0L)
        k <- max(nr)
    }

    # which transition each row  of newX represents
    transition <- rep(ustrata, n.perstrat)

    # remove any rows where X is missing
    #  these arise when a variable is used only for some transitions
    #  the row of data needs to be tossed for the given ones, but will be
    #    okay for other transitions
    keep <- !apply(is.na(newX), 1, any)
    if (!all(keep)) {
        newX <- newX[keep,, drop=FALSE]
        rindex <- rindex[keep]
        newstat <- newstat[keep]
        transition <- transition[keep]
    }

    if (ncol(Y) ==2) newY <- Surv(Y[rindex,1], newstat)
    else newY <- Surv(Y[rindex,1], Y[rindex,2], newstat)

    # new strata, add on any implicit strata in the data.
    if (is.null(strata)) newstrat <- ustrata[transition]
    else {
        # if the old strata is 1, 2, 3 and states are 1-5, the new ones will 
        #  11-15 for strata 1, 21-25 for strata 2, etc.
        # (this is always called with 1,2,3... for the strata)
        # there are usually <10 transitions, so this makes simple labels
        temp <-10^ ceiling(log(length(ustrata), 10))
        newstrat <- ustrata[transition] + strata[rindex]*temp
    } 

    # give variable names to the new data
    vname <- rep("", ncol(newX))
    ctemp <- cmap[-1,,drop=FALSE]
    vname[ctemp[ctemp>0]] <- colnames(X)[row(ctemp)[ctemp>0]]
    colnames(newX) <- vname

    list(X=newX, Y=newY, strata=as.integer(newstrat), 
         transition= as.integer(transition), rindex=rindex)
}
