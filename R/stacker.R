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
stacker <- function(cmap, smap, istate, X, Y, strata, states) {
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
    #  up otherwise.  The first row of smap gives the strata for each change
    ustrata <- unique(smap[1,])
    nstrat <- length(ustrata)
    n.perstrat <- integer(nstrat)
    for (i in 1:nstrat) {
        itemp <- unique(from.state[smap[1,] == ustrata[i]])
        n.perstrat[i] <- sum(istate %in% itemp)
    }
    
    # The constructed X matrix has a block or rows for each ustrata level
    n2 <- sum(n.perstrat)  # number of rows in new data
    newX <- matrix(0, nrow=n2, ncol=max(cmap))
    k <- 0
    rindex <- integer(n2)   # original row for each new row of data
    newstat <- integer(n2)  # new status
    for (i in 1:nstrat) {
        whichcol <- which(smap[1,] == ustrata[i])  # cols of cmap to look at
        subject <- which(istate %in% from.state[whichcol]) # data rows in strata
        nr <- k + seq(along=subject)  # rows in the newX for this strata
        rindex[nr] <- subject
        # Fill in X one transition at a time
        for (j in whichcol) {
            j1 <- which(istate == from.state[j]) # rows of X in transition
            j2 <- which(istate[subject] == from.state[j]) # new rows
            nc <- cmap[,j]             # variables in this transition
            newX[nr[j2], nc[nc>0]] <- X[j1, nc>0] # rows of cmap = cols of X
        }
        
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

    # newstrat should be an integer vector, which can be used for the interal C calls
    newstrat <- ustrata[transition]
    if (is.matrix(strata)){
        # this is the most complex case.  Some transitions use some columns of
        # istrat, and some use others
        maxstrat <- apply(strata, 2, max)  # max in each colum of strata
        mult <- cumprod(c(1, maxstrat))
        temp <- max(mult) * newstrat
        for (i in 1:ncol(strata)) {
            k <- smap[i+1, transition]
            temp <- temp + ifelse(k ==0, 0L, strata[i, rindex]* temp[i] -1L)
        } 
        newstrat <- match(temp, sort(unique(temp)))
    }
    else if (length(strata) > 0) {
        # strata will be an integer vector, values from 1 to number of strata
        mult <- max(strata)
        temp <- mult * newstrat + ifelse(smap[2,transition]==0, 0L, strata[rindex] -1L)
        newstrat <- match(temp, sort(unique(temp)))
    }       
 
    # give variable names to the new data
    vname <- rep("", ncol(newX))
    vname[cmap[cmap>0]] <- colnames(X)[row(cmap)[cmap>0]]
    colnames(newX) <- vname

    list(X=newX, Y=newY, strata=as.integer(newstrat), 
         transition= as.integer(transition), rindex=rindex)
}
