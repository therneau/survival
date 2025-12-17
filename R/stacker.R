#
# This routine creates a stacked data set.
# The key input is the cmap matrix, which has one row for each column 
#  of X and one column per transition.  It may have extra rows if there are
#  proportional baseline hazards
# The first row of smat contains state to strata information.
# Input data is X, Y, strata, and initial state (integer).
# The model.frame is used only for strata, currently a rare case in multistate
#
# For each transition the expanded data has a set of rows, all those whose
#  initial state makes them eligible for the transition.  
# Strata is most often null; it encodes a users strata() addition(s). Such terms
#  occur less often in multistate models (in my experience so far.)
#
stacker <- function(cmap, smap, istate, X, Y, mf, states, dropzero=TRUE) {
    from.state <- as.numeric(sub(":.*$", "", colnames(cmap)))
    to.state   <- as.numeric(sub("^.*:", "", colnames(cmap)))

    # just in case cmap has columns I don't need (I don't think this can
    #  happen
    check <- match(from.state, istate, nomatch=0)
    if (any(check==0)){
        # I think that this is impossible
        warning("extra column in cmap, this is a bug")  # debugging line
        # browser()
        cmap <- cmap[,check>0, drop=FALSE]
        smap <- smap[,check>0, drop=FALSE]
        from.state <- from.state[check>0]
        to.state <- to.state[check>0]
    }

    # Don't create X and Y matrices for transitions with no covariates, for
    #  coxph calls.  But I need them for survfit.coxph (at present)
    zerocol <- apply(cmap==0, 2, all)
    if (dropzero && any(zerocol)) {
        cmap <- cmap[,!zerocol, drop=FALSE] 
        smap <- smap[,!zerocol, drop=FALSE]
        smap[,] <- match(smap, sort(unique(c(smap)))) # relabel as 1, 2,...
        from.state <- from.state[!zerocol]
        to.state <- to.state[!zerocol]
    }
        
    endpoint <- c(0, match(attr(Y, "states"), states))
    endpoint <- endpoint[ 1 + Y[,ncol(Y)]]  # endpoint of each row, 0=censor

    # Jan 2021: changed from looping once per strata to once per transition.
    #  Essentially, a block of data for each unique column of cmap.  If two
    #  of those columns have the same starting state, it makes me nervous
    #  (statistically), but forge onward and sort the issues out in the
    #  fits.
    # Dec 2024: change to once per stratum.  The issue was multiple
    #  transitions like 1:4 and 1:5 in the same stratum, which led to the
    #  same row of data twice in one stratum; which is not valid.
    # May 2025: further update in conjunction with "Shared coefficients
    #  and shared baselines in multistate models", which forced a re-thinking
    #  of the process.  For the use cases found there, the code below is now
    #  correct.  For others, the data produced may be questionable;
    #  we await an actual use case to work things out.
    # The constructed X matrix will have a block of rows for each column of
    #  cmap such that: the column is not 0, it is not a duplicate of another
    #  cmap column that is in the same stratum.  
    # Within a stratum we retain sort order of the data

    sgrp   <- match(smap[1,], unique(smap[1,]))  # the stratum group
    nblock <- max(sgrp)  # total number of blocks
    # Pass 1 to find the total data set size
    n.perblock <- integer(nblock)
    cgrp <- vector("list", nblock)
    for (i in 1:nblock) {
        j <- which(sgrp == i)  # which cols of cmap & smap apply
        ctemp <- !duplicated(t(cmap[,j, drop=FALSE])) # unique coef sets in sgrp
        cgrp[[i]] <- lapply(j[which(ctemp)],
                       function(i) {
                           temp <- cmap[,j,drop=FALSE] - cmap[,i]
                           j[which(apply(temp, 2, function(x) all(x==0)))]})
        # cgrp[[i]] will be a list, 99% of the time it will be of length 1: 
        #  a vector containing all the columns of cmap that go with this strata
        #  This is the case found in the vignette (with multiple subcases):
        #  everyone at risk appears once in the stratum
        # if cgrp has length >1 the resulting risk sets have dups of some or
        #  all obs, but with different coefficient patterns.
        temp <- 0
        for (k in cgrp[[i]])
            temp <- temp + sum(istate %in% from.state[k])
        n.perblock[i] <- temp
    }
    
    # The constructed X matrix has a block of rows for each stratum
    n2 <- sum(n.perblock)  # number of rows in new data
    newX <- matrix(0, nrow=n2, ncol=max(cmap))
    k <- 0
    rindex <- integer(n2)   # original row for each new row of data
    newstat <- integer(n2)  # new status
    Xcols   <- ncol(X)      # number of columns in X
    for (i in 1:nblock) {
        for (j in cgrp[[i]]) {
            subject <- which(istate %in% from.state[j]) # data rows in strata
            nr <- k + seq(along.with =subject)  # rows in the newX for subblock
            rindex[nr] <- subject
            nc <- cmap[,min(j)]  # all cols j of cmap are the same
            if (any(nc > Xcols)) { # constructed PH variables
                newX[nr, nc[nc>Xcols] ] <- 1
                nc <- nc[1:Xcols]
            }
            newX[nr, nc[nc>0]] <- X[subject, which(nc>0)] #row of cmap= col of X
        
            event.that.counts <- (endpoint[subject] %in% to.state[j])
            newstat[nr] <- ifelse(event.that.counts, 1L, 0L)
            k <- max(nr)
        }
    }

    # which (grouped) transition each row of newX represents
    transition <- rep(1:nblock, n.perblock)
    # create the remaining strata, if needed
    newstrat <- factor(colnames(smap)[transition], colnames(smap))
    if (nrow(smap) >1) {
        tmap <- smap[-1,, drop=FALSE]  #ignore (Baseline) row
        # if all the states have the same strata variables, things are a
        # bit easier
        allsame <- TRUE
        for (i in 1:ncol(tmap)) 
            if (any(tmap[,i] != tmap[,1])) allsame <- FALSE
        if (allsame) temp <- do.call(paste, c(mf[,tmap[,1]], sep='.'))
        else {
            rtemp <- split(rindex, rep(1:nblock, n.perblock)) #rows per trans
            temp <- vector("list", nblock)
            for (i in 1:nblock)
                temp[[i]] <- do.call(paste, 
                                     c(mf[rtemp[[i]], tmap[,i]], sep='.'))
        }
        newstrat <- factor(paste(newstrat, unlist(temp), sep='.'))
    } 

    # remove any rows where X is missing
    #  these arise when a variable is used only for some transitions
    #  the row of data needs to be tossed for the given ones, but will be
    #  okay for other transitions which do not use the offending variable.
    keep <- !apply(is.na(newX), 1, any)
    if (!all(keep)) {
        newX <- newX[keep,, drop=FALSE]
        rindex <- rindex[keep]
        newstat <- newstat[keep]
        transition <- transition[keep]
        newstrat <- newstrat[keep]
    }

    if (ncol(Y) ==2) newY <- Surv(Y[rindex,1], newstat)
    else newY <- Surv(Y[rindex,1], Y[rindex,2], newstat)


    # If there were multiple unique cmap columns within one strata, then an
    #  observation might appear twice, and we might no longer be sorted.
    # (This is also the case where we do not yet know what a 'sensible'
    # result of stacker should be, i.e., very rare.)
    if (any(sapply(cgrp, length)) > 1) {
        indx <- order(newstrat, rindex)
        if (any(diff(indx) !=1)) {
            newX <- newX[indx,, drop=FALSE]
            rindex <- rindex[indx]
            newstat <- newstat[indx]
            transition <- transition[indx]
        }
    }
                      
    #
    # give variable names to the new data  (some names get used more than once)
    #    vname <- rep("", ncol(newX))
    #    vname[cmap[cmap>0]] <- colnames(X)[row(cmap)[cmap>0]]
    first <- match(sort(unique(cmap[cmap>0])), cmap) #first instance of each value
    vname <- rownames(cmap)[row(cmap)[first]]
    colnames(newX) <- vname
    list(X=newX, Y=newY, strata=as.integer(newstrat), 
         transition= as.integer(transition), rindex=rindex)
}
