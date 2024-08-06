# residuals.survfit for an AJ model.  The file residuals.survfit was broken
#  in two for readability

rsurvpart2 <- function(Y, casewt, istate, times, type, fit, method) {
    ny <- ncol(Y)
    n  <- nrow(Y)
    ntime <- length(times)
    nstate <- length(fit$states)
    
    # ensure that Y, istate, and fit all use the same set of states
    states <- fit$states
    if (!identical(attr(Y, "states"), fit$states)) {
        map <- match(attr(Y, "states"), fit$states)
        Y[,ny] <- c(0, map)[1+ Y[,ny]]    # 0 = censored
        attr(Y, "states") <- fit$states
    }
    if (is.null(istate)) istate <- rep(1L, nrow(Y)) #everyone starts in s0
    else {
        if (is.character(istate)) istate <- factor(istate)
        if (is.factor(istate)) {
            if (!identical(levels(istate), fit$states)) {
                map <- match(levels(istate), fit$states)
                if (any(is.na(map))) stop ("invalid levels in istate")
                istate <- map[istate]
            }       
        } 
    }  # else istate is numeric, we take what we get and hope it is right

    # Compute the initial leverage
    p0 <- fit$p0         
    inf0 <- NULL
    if (is.null(fit$call$p0) && any(p0 < 1)) { 
        #p0 was not supplied by the user, and the intitial states vary
        inf0 <- matrix(0., nrow=nrow(Y), ncol=nstate)
        p0 <- fit$p0
        t0 <- fit$t0
        if (ny==2) at.zero <- 1:n   #everyone at risk
        else  at.zero <- which(Y[,1] < t0 & Y[,2] >= t0)
        for (j in 1:nstate) {
            inf0[at.zero, j] <- (ifelse(istate[at.zero]==states[j], 1, 0) -
                                     p0[j])/sum(casewt[at.zero])
            }
    }

    t0 <- fit$t0
    # get the values of from:to for each transition
    temp <- colnames(fit$cumhaz)
    from <- as.integer(sub(":[0-9]*$", "", temp))
    to   <- as.integer(sub("^[0-9]*:", "", temp))
 
    if (type== "cumhaz") {
        events <- (rowSums(fit$n.event) >0)  # minor speedup, only event times
        hazard <- apply(rbind(0, fit$cumhaz[events,]), 2, diff)
        nrisk  <- fit$n.risk[events,]
        dtime  <- fit$time[events]

        # Create the index of the reporting times into the hazard curve
        #   tindex = largest event time <= reporting time
        #   yindex = largest event time <= per-obs event/censor time
        #   sindex = largest event time <= per-obs entry time
        tindex <- findInterval(times, dtime, left.open=FALSE)
        yindex <- findInterval(Y[, ny-1], dtime, left.open=FALSE)
        if (ny==3) sindex <- findInterval(Y[,1], dtime, left.open=FALSE)
    
        # A common operation is that resid[i,j] is updated using 
        #    xxx[min(yindex[i], tindex[j])] for some vector xxx, or
        #    xxx[min(yindex[i], tindex[j])] - xxx[min(sindex[i], tindex[j])]
        #   for (time1, time2) data.
        # For the dN term the rule is 
        #   if (death and tindex >= yindex) then yindex  else 0
        # i.e. for any row dN applies to all reporting times at or
        #  after the death time for that observation
        status <- Y[,ny]
        ymin <- outer(yindex, tindex, pmin)
        dmin <- outer(ifelse(status ==0, 0L, yindex), tindex,
                      function(a, b) ifelse(a==0 | a>b, 0L, a)) 
        if (ny==3) smin <- outer(sindex, tindex, pmin)

        nhaz <- ncol(fit$cumhaz)  # the number of possible transitions
        D <- array(0, dim=c(nrow(Y), nhaz, ntime))

        # Each column of hsum is cumsum (lambda /nrisk), for that transition
        safe <- ifelse(nrisk==0, 1, nrisk) # avoid 0/1
        hsum <- apply(rbind(0, hazard/safe[,from]), 2, cumsum) 
        for (k in 1:nhaz) { 
            atrisk <- (as.integer(istate)==from[k])
            event <- (atrisk & status==to[k]) # an event of this type
            # The dN part of the residual only apples to 'event' rows,
            #   dmin control which columns it adds into
            # The dLambda part applies to 'atrisk' rows, ymin and smin tell
            #   which portion in which columns
            D[event,k,] <- c(0, 1/safe[,from[k]]) [1+ dmin[event,]]
            term2 <- hsum[1+ ymin[atrisk, ],k]            # the d\Lambda part
            if (ny ==2) D[atrisk,k,] <- D[atrisk,k,] - term2
            else {
                # events happen at the end of an interval, so no dN at the start
                term3 <- hsum[1 + smin[atrisk, ],k]
                D[atrisk,k,] <- D[atrisk,k,] + term3 - term2
            } 
        }
    } else if (TRUE) {
        # Compute the result using the direct method, in C code
        # the routine is called separately for each curve, data in sorted order
        #
        is1 <- as.integer(istate) -1L  # 0 based subscripts for C
        if (is.null(inf0)) inf0 <- matrix(0, nrow=nrow(Y), ncol=nstate)
 
        if (ny==2) asort1 <- 0L else asort1 <- order(Y[,1]) -1L
        asort2 <- order(Y[,ny-1]) -1L
        storage.mode(casewt) <- "double"  # the data set might have had integers
        if (TRUE) {   
            tfit <- .Call(Csurvfitresid, Y, asort1, asort2, is1, 
                          casewt, p0, inf0, times, fit$t0, 
                          type== "auc")
        } else { # more efficient(?) - still in development
            event <- (Y[,3] != 0) # the non-censored rows
            etime <- sort(unique(Y[event,2]))
            hindex <- matrix(0L, nstate, nstate)
            hindex[cbind(from, to)] <- seq_along(from) -1L
            tindex <- cbind(from, to) - 1L
            # the max number at risk for any state, used for the skiplists
            if (is.null(fit$counts)) { # integer data
                maxrisk <- as.integer(apply(fit$n.risk,2, max))
            } else {
                temp <- fit$counts[, grepl("nrisk", colnames(fit$counts))]
                maxrisk <- as.integer(apply(temp, 2, max))
            }
            # CMD check doesn't like the reference to a not-there-yet C
            #tfit <- .Call(Csurvfitresid2, Y, asort1, asort2, etime, is1, 
            #              times, fit$n.transition, fit$n.risk, p0, inf0,
            #              hindex, tindex, maxrisk)
        }
            
        if (ntime==1) {
            if (type=="auc") D <- tfit[[2]] else D <- tfit[[1]]
        }
        else {
            if (type=="auc") D <- array(tfit[[2]], dim=c(nrow(Y), nstate, ntime))
            else         D <- array(tfit[[1]], dim=c(nrow(Y), nstate, ntime))
        }
    }   
    D
}
