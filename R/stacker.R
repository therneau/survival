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
    # June 2026: Further thought has shown that for every practical use case
    #  that I know, the resulting shared strata should have exactly 1 copy of
    #  each row that was at risk for one of the shared transitions.  Since 
    #  every row has a single istate, the simple way to handle this is to
    #  select istate in list of "from" states.
    # July 2026: Carrying through a difficult case, carefully, I now think 
    #  that the 2021 argument is the right one. We have come full circle.
    # The methods vignette has more discussion. See the 'shared hazards' secton.

    # Don't create X and Y matrices for transitions with no covariates, for
    #  coxph calls (dropzero TRUE). But I need them for survfit.coxphms.
    if (dropzero) keepcol <- apply(cmap>0, 2, any)
    else keepcol <- rep(TRUE, ncol(cmap))
    baseline <- smap[1,]
    sgrp <- match(baseline[keepcol], baseline)
    ugrp <- unique(sgrp)
    nblock <- length(ugrp)  # total number of blocks

    # Pass 1 to find the total data set size
    n.perblock <- integer(nblock)
    cgrp <- vector("list", nblock)
    icount <- table(istate)
    for (i in 1:nblock) {
        k <- which(baseline == ugrp[i]) 
        cgrp[[i]] <- k
        # cgrp[[i]] will be a vector, 99% of the time there are no shared
        #  baseline hazards and it will be of length 1
        #
        n.perblock[i] <- sum(icount[from.state[k]])
    }
    
    # The constructed X matrix has a block of rows for each stratum
    n2 <- sum(n.perblock)  # number of rows in new data
    newX <- matrix(0, nrow=n2, ncol=max(cmap))
    rindex <- integer(n2)   # original row for each new row of data
    # Often there are no shared baseline hazards (do duplicates in baseline)
    #  and in that case shared and subshare will be identical, and every
    #  element of submap a list of length 1
    shared   <- rep(1:nblock, n.perblock) # the shared baseline number
    subshare <- integer(n2)
    newstat <- integer(n2)  # new status
    Xcols <- 1:ncol(X)
    hasph <- nrow(cmap) > ncol(X)  # there are constructed PH variables
    kk <- 0
    for (i in 1:nblock) {
        j <- cgrp[[i]]
        # do separate push for each unique istate
        #jinit <- unique(from.state[j])
        #for (k in jinit) {
        for (k in j) {
            subject <- which(istate %in% from.state[k]) # data rows in strata
            nr <- kk + seq(along.with =subject)  # rows in the newX for subblock
            rindex[nr] <- subject

            subshare[nr] <- k
            nc <- cmap[Xcols,k] # coefficients for this subshare
            newX[nr, nc] <- X[subject, which(nc>0)]
            if (hasph && any(cmap[-Xcols,k]>0)) {
                # constructed PH variables aren't in X, but act like x=1
                newX[nr, cmap[-Xcols,k]] <- 1
            }
            event.that.counts <- (endpoint[subject] %in% to.state[k])
            newstat[nr] <- ifelse(event.that.counts, 1L, 0L)
        
            kk <- max(nr)
        }
    }

    # newstrat is the strata the the coxph maximizer will use
    newstrat <- shared 
    if (nrow(smap) >1) {
        # there are strata in the call as well, so expand the strata to
        #  "strata from the model".1:2, etc
        sstrat <- function(...) strata(..., shortlabel=TRUE)
        tmap <- smap[-1,, drop=FALSE]  #ignore (Baseline) row
        # rtemp will be list, row numbers for subshare 1, rows for 2,..
        rtemp <- split(rindex, rep(1:nblock, n.perblock)) #rows per trans
        temp <- vector("list", nblock) # one per subshare
        for (i in 1:nblock){
            j <- (tmap[,i]>0) # which stata vars for this subshare?
            if (sum(j) ==1) {
                zz <- mf[[(rownames(tmap))[j]]] # the strata
                temp[[i]] <- zz[rtemp[[i]]]         
            } else if (sum(j)>1) { 
                zz <- mf[rtemp[[i]], (rownames(tmap))[j]] 
                temp[[i]] <- do.call(sstrat, zz)
            } else temp[[i]] <- rep(0, n.perblock[i]) #no prior strata    
        }
        newstrat <- do.call(sstrat, list(newstrat, unlist(temp)))
    } 

    # remove any rows where X or strata is missing
    #  these arise when a variable is used only for some subshares
    #  the row of data needs to be tossed for the given ones, but will be
    #  okay for other subshares which do not use the offending variable.
    keep <- !apply(is.na(newX), 1, any) & !is.na(newstrat)
    if (!all(keep)) {
        newX <- newX[keep,, drop=FALSE]
        rindex <- rindex[keep]
        newstat <- newstat[keep]
        subshare <- subshare[keep]
        shared <- shared[keep]
        newstrat <- newstrat[keep]
    }

    if (ncol(Y) ==2) newY <- Surv(Y[rindex,1], newstat)
    else newY <- Surv(Y[rindex,1], Y[rindex,2], newstat)

    #
    # give variable names to the new data  (some names get used more than once)
    #
    first <- match(sort(unique(cmap[cmap>0])), cmap) #first instance of each 
    vname <- rownames(cmap)[row(cmap)[first]]
    colnames(newX) <- vname
    #  There is a block for each transition (hazard), but if there are
    #  shared hazards they are not necessarily in the order 1, 2, 3, ...
    #  If there are no shared transtions hazard never gets used.
    list(X=newX, Y=newY, strata=as.integer(newstrat), rindex=rindex,
         shared = shared,  hazard=subshare)
}
