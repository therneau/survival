#
# This routine creates a stacked data set.
# The key input is the cmap matrix, which has one row for each column 
#  of X and one column per transition.  It may have extra rows if there are
#  proportional baseline hazards
# The first row of smat contains state to strata information.
# Input data is X, Y, strata, and initial state (integer).
#
# For each transition the expanded data has a set of rows, all those whose
#  initial state makes them eligible for the transition.  
# Strata is most often null; it encodes a users strata() addition(s); less
#  often used in multistate.
#
stacker <- function(cmap, smap, istate, X, Y, strata, states) {
    from.state <- as.numeric(sub(":.*$", "", colnames(cmap)))
    to.state   <- as.numeric(sub("^.*:", "", colnames(cmap)))

    # just in case cmap has columns I don't need (I don't think this can
    #  happen
    check <- match(from.state, istate, nomatch=0)
    if (any(check==0)){
        # I think that this is impossible
        warning("extra column in cmap, this is a bug")  # debugging line
        # browser()
        cmap <- cmap[,check>0]
        from.state <- from.state[check>0]
        to.state <- to.state[check>0]
    }

    # Don't create X and Y matrices for transitions with no covariates
    zerocol <- apply(cmap==0, 2, all)
    if (any(zerocol)) {
        cmap <- cmap[,!zerocol, drop=FALSE] 
        smap <- smap[,!zerocol, drop=FALSE]
        smap[,] <- match(smap, sort(unique(c(smap)))) # relabel as 1, 2,...
    }
        
    endpoint <- c(0, match(attr(Y, "states"), states))
    endpoint <- endpoint[ 1 + Y[,ncol(Y)]]  # endpoint of each row, 0=censor

    # Jan 2021: changed from looping once per strata to once per transition.
    #  Essentially, a block of data for each unique column of cmap.  If two
    #  of those columns have the same starting state, it makes me nervous
    #  (statistically), but forge onward and sort the issues out in the
    #  fits.
    # Pass 1 to find the total data set size
    nblock <- ncol(cmap)
    n.perblock <- integer(nblock)
    for (i in 1:nblock) {
        n.perblock[i] <- sum(istate == from.state[i]) # can participate
    }
    
    # The constructed X matrix has a block of rows for each column of cmap
    n2 <- sum(n.perblock)  # number of rows in new data
    newX <- matrix(0, nrow=n2, ncol=max(cmap))
    k <- 0
    rindex <- integer(n2)   # original row for each new row of data
    newstat <- integer(n2)  # new status
    Xcols   <- ncol(X)      # number of columns in X
    for (i in 1:nblock) {
        subject <- which(istate == from.state[i]) # data rows in strata
        nr <- k + seq(along.with =subject)  # rows in the newX for this strata
        rindex[nr] <- subject
        nc <- cmap[,i]  
        if (any(nc > Xcols)) { # constructed PH variables
            newX[nr, nc[nc>Xcols] ] <- 1
            nc <- nc[1:Xcols]
        }
        newX[nr, nc[nc>0]] <- X[subject, which(nc>0)] # row of cmap= col of X
        
        event.that.counts <- (endpoint[subject] == to.state[i])
        newstat[nr] <- ifelse(event.that.counts, 1L, 0L)
        k <- max(nr)
    }

    # which transition each row of newX represents
    transition <- rep(1:nblock, n.perblock)

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

    # newstrat will be an integer vector.
    newstrat <- smap[1, transition]   # start with strata induced by multi-state
    # then add any strata from the users strata() terms
    if (is.matrix(strata)){
        # this is the most complex case. 
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
        # strata will be an integer vector with elements of 1, 2 etc 
        mult <- max(strata)
        temp <- mult * newstrat + ifelse(smap[2,transition]==0, 0L, strata[rindex] -1L)
        newstrat <- match(temp, sort(unique(temp)))
    }       
 
    # give variable names to the new data  (some names get used more than once)
#    vname <- rep("", ncol(newX))
#    vname[cmap[cmap>0]] <- colnames(X)[row(cmap)[cmap>0]]
    first <- match(sort(unique(cmap[cmap>0])), cmap) #first instance of each value
    vname <- rownames(cmap)[row(cmap)[first]]
    colnames(newX) <- vname
    list(X=newX, Y=newY, strata=as.integer(newstrat), 
         transition= as.integer(transition), rindex=rindex)
}
