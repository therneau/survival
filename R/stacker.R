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
    # June 2026: if 1:3, 2:3 and 1:4 were shared strata, there are 3 blocks
    #  First all those at risk for 1:3: istate==1, then all those at risk
    #  for 2:3: istate==2, then at risk for 1:3: istate==1. If the user included
    #  a +shared arg the resulting model could possibly make sense?
    #  Also return an transition vector containing the transition number of 
    #  each block. This is different from strans = the shared one and strata=
    #  strans + user strata; the likelihood computation needs the last of 
    #  these.
    # Within a given hazard we retain sort order of the data

    # Don't create X and Y matrices for transitions with no covariates, for
    #  coxph calls.  But I need them for survfit.coxphms (at present)
    # June 2026 -- I not longer think so
    zerocol <- apply(cmap==0, 2, all)
    baseline <- smap[1,]
    if (dropzero && any(zerocol)) {
        sgrp <- match(baseline[!zerocol], unique(baseline))
    } else sgrp   <- match(baseline, unique(baseline)) 
#        cmap <- cmap[,!zerocol, drop=FALSE] 
#        smap <- smap[,!zerocol, drop=FALSE]
#        smap[,] <- match(smap, sort(unique(c(smap)))) # relabel as 1, 2,...
#        from.state <- from.state[!zerocol]
#        to.state <- to.state[!zerocol]
#    }
#     sgrp   <- match(smap[1,], unique(smap[1,]))  # the stratum group
    ugrp <- unique(sgrp)
    nblock <- length(ugrp)  # total number of blocks
    # Pass 1 to find the total data set size
    n.perblock <- integer(nblock)
    cgrp <- vector("list", nblock)
    icount <- table(istate)
    for (i in 1:nblock) {
        cgrp[[i]] <- which(baseline == ugrp[i]) 
        # cgrp[[i]] will be a list, 99% of the time there are no shared
        #  baseline hazards and it will be of length 1
        #
        temp <- 0
        for (k in cgrp[[i]])
            temp <- temp + icount[from.state[k]]
        n.perblock[i] <- temp
    }
    
    # The constructed X matrix has a block of rows for each stratum
    n2 <- sum(n.perblock)  # number of rows in new data
    newX <- matrix(0, nrow=n2, ncol=max(cmap))
    k <- 0
    rindex <- integer(n2)   # original row for each new row of data
    newstat <- integer(n2)  # new status
    stran   <- rep(1:nblock, n.perblock) # the shared transition number
    transition <- unlist(lapply(cgrp, function(i) {
        rep(i, icount[from.state[i]])}))
    Xcols <- ncol(X)
    for (i in 1:nblock) {
        for (j in cgrp[[i]]) {
            subject <- which(istate %in% from.state[j]) # data rows in strata
            nr <- k + seq(along.with =subject)  # rows in the newX for subblock
            rindex[nr] <- subject
            nc <- cmap[,j] # coefficients for this transition
            newX[nr, nc] <- X[subject, which(nc>0)]
            if (any(nc > Xcols)) { # constructed PH variables
                newX[nr, nc[nc>Xcols] ] <- 1
                nc <- nc[1:Xcols]
            }
        
            event.that.counts <- (endpoint[subject] %in% to.state[j])
            newstat[nr] <- ifelse(event.that.counts, 1L, 0L)
            k <- max(nr)
        }
    }

    # which (grouped) transition each row of newX represents
    newstrat <- stran
    # newstrat will be 1:2, ... i.e. the from:to pairs
    if (nrow(smap) >1) {
        # there are strata in the call as well, so expand the strata to
        #  "strata from the model".1:2, etc
        sstrat <- function(...) strata(..., shortlabel=TRUE)
        tmap <- smap[-1,, drop=FALSE]  #ignore (Baseline) row
        # rtemp will be list, row numbers for transition 1, rows for 2,..
        rtemp <- split(rindex, rep(1:nblock, n.perblock)) #rows per trans
        temp <- vector("list", nblock) # one per transition
        for (i in 1:nblock){
            j <- (tmap[,i]>0) # which stata vars for this transition?
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
    #  these arise when a variable is used only for some transitions
    #  the row of data needs to be tossed for the given ones, but will be
    #  okay for other transitions which do not use the offending variable.
    keep <- !apply(is.na(newX), 1, any) & !is.na(newstrat)
    if (!all(keep)) {
        newX <- newX[keep,, drop=FALSE]
        rindex <- rindex[keep]
        newstat <- newstat[keep]
        transition <- transition[keep]
        stran <- stran[keep]
        newstrat <- newstrat[keep]
    }

    if (ncol(Y) ==2) newY <- Surv(Y[rindex,1], newstat)
    else newY <- Surv(Y[rindex,1], Y[rindex,2], newstat)

    #
    # give variable names to the new data  (some names get used more than once)
    #    vname <- rep("", ncol(newX))
    #    vname[cmap[cmap>0]] <- colnames(X)[row(cmap)[cmap>0]]
    first <- match(sort(unique(cmap[cmap>0])), cmap) #first instance of each value
    vname <- rownames(cmap)[row(cmap)[first]]
    colnames(newX) <- vname
    list(X=newX, Y=newY, strata=as.integer(newstrat), 
         transition= factor(transition, 1:ncol(smap), colnames(smap)),
         rindex=rindex,
         stran = stran)
}
