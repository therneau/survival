# The routine for Aalen-Johansen estimates
# Major documentation is in the methods vignette
survfitAJ <- function(X, Y, weights, id, cluster, robust, istate,
                      stype=1, ctype=1,
                      se.fit=TRUE,
                      conf.int= .95,
                      conf.type=c('log',  'log-log',  'plain', 'none', 
                                  'logit', "arcsin"),
                      conf.lower=c('usual', 'peto', 'modified'),
                      influence = FALSE, start.time, p0, type, entry=FALSE,
                      time0){

    # This first chunk is for backwards compatability
    if (!missing(type)) {
        if (!missing(ctype) || !missing(stype))
            stop("cannot have both an old-style 'type' argument and the stype/ctype arguments that replaced it")
        if (!is.character(type)) stop("type argument must be character")
        # older style argument is allowed
        temp <- charmatch(type, c("kaplan-meier", "fleming-harrington", "fh2"))
        if (is.na(temp)) stop("invalid value for 'type'")
        type <- c(1,3,4)[temp]
    }
    else {
        if (!(ctype %in% 1:2)) stop("ctype must be 1 or 2")
        if (!(stype %in% 1:2)) stop("stype must be 1 or 2")
        type <- as.integer(2*stype + ctype  -2)
    }

    # The various types are not implemented for AJ, but I want common arguments
    # with survfitKM
    if (type != 1 || ctype !=1)
        warning("only stype=1, ctype=1 implimented for multi-state data")

    conf.type <- match.arg(conf.type)
    conf.lower<- match.arg(conf.lower)
    if (conf.lower != "usual") 
        warning("conf.lower is ignored for multi-state data")
    if (is.logical(conf.int)) {
        # A common error is for users to use "conf.int = FALSE"
        #  it's illegal per documentation, but be kind
        if (!conf.int) conf.type <- "none"
        conf.int <- .95
    }

    # survfitKM allows for 0,1, 2, or 3; we just treat it as yes/no
    if (is.logical(influence)) {
        influence <- ifelse(influence, 1L, 0L)
        # numeric is treated as all or nothing
        if (!influence) influence <- 0L
        else influence <- 3L
    }
    else if (!is.numeric(influence))
        stop("influence argument must be numeric or logical")
    influence <- min(influence, 1L)  # make it 0/1
  
    if (!se.fit) {
        # if the user asked for no standard error, skip any robust computation
        ncluster <- 0L
        influence <- 0L
    }

    type <- attr(Y, "type")
    ny <- ncol(Y)
    # This line should be unreachable, unless they call "surfitAJ" directly
    if (type !='mright' && type!='mcounting')
         stop(paste("multi-state computation doesn't support \"", type,
                          "\" survival data", sep=''))

    # Use survcheck to validate the data   
    if (length(id) ==0) id <- seq(along=X)  # a dummy value for later work
    if (missing(istate) || is.null(istate))
        mcheck <- survcheck2(Y, id)  
    else mcheck <- survcheck2(Y, id, istate)
    if (mcheck$flag["overlap"] > 0)
        stop("a subject has overlapping time intervals")
    if (any(mcheck$flag > 0)) stop("one or more flags are >0 in survcheck")
    # We may eventually remove the line above, then the line below is active
    if (mcheck$flag["gap"] > 0 || mcheck$flag["jump"] > 0)
        warning("subject(s) with time gaps, results are questionable")

    istate <- mcheck$istate  # this copy will have the levels set properly
    states <- mcheck$states
    nstate <- length(states) 
    if (!identical(states, attr(Y, "states"))) { # new states are a superset
        smap <- c(0, match(attr(Y, "states"), states))
        Y[,ncol(Y)] <- smap[Y[,ncol(Y)] +1]  
        attr(Y, 'states') <- states
    }       

    # Now that we know the states, verify that p0 is correct (if present)
    if (!missing(p0) && length(p0) >0) {
        if (length(p0) != nstate) stop("wrong length for p0")
        if (!is.numeric(p0) || !isTRUE(all.equal(sum(p0), 1)))
            stop("p0 must be a numeric vector that adds to 1")
    } else p0 <- NULL

    # If there is a start.time directive, remove any data rows that
    #  are completely prior to that time
    n <- nrow(Y)
    if (!missing(start.time)) {
        if (!is.numeric(start.time) || length(start.time) !=1
            || !is.finite(start.time))
            stop("start.time must be a single numeric value")
        cmax <- tapply(Y[,ny-1], X, max) # max time in each curve
        if (any(cmax <= start.time)) 
            stop("start.time has removed all the observations from at least one curve")
        toss <- which(Y[,ny-1] < start.time) # remove the obs
        if (length(toss) >0 ) {
            # why not just set Y[,3] = censored?  Because a few lines down
            #  survcheck will then be unhappy, and rightfully so.
            Y <- Y[-toss,,drop=FALSE]
            X <- X[-toss]
            weights <- weights[-toss]
            if (length(id) ==n) id <- id[-toss]
            istate <- istate[-toss]
        }
    n <- nrow(Y)
    }

    
    # The user can call with cluster, id, robust or any combination.
    # If only id, treat it as the cluster too
    if (missing(robust) || length(robust)==0) robust <- TRUE
    if (!robust) stop("multi-state survfit supports only a robust variance")

    has.cluster <-  !(missing(cluster) || length(cluster)==0) 
    has.id <-       !(missing(id) || length(id)==0)

    # Residuals will be in the order of integer(id), and we want them to
    #  be in the same order as the data.  So prevent the default sorted levels
    if (has.id) id <- factor(id, unique(id))
    else  {
        if (ncol(Y) ==3) stop("an id statement is required for start,stop data")
        else id <- 1:n
    }
    if (influence && !(has.cluster || has.id)) {
        cluster <- seq(along.with=X)
        has.cluster <- TRUE
    }

    if (has.cluster) {
        cluster <- factor(cluster, unique(cluster)) #see id comment above
        clname <- levels(cluster)
        cluster <- as.integer(cluster)
        ncluster <- length(clname)
        # If there is both id and cluster, check that any given id lies entirely
        # in a single cluster.  The converse is very likely a data error.
        # table(id, cluster) is the quick way, but could generate a huge sparse
        #  matrix that fills memory.  Use a small C routine.
        if (has.id) {
            twoclust <- .Call(Ctwoclust, as.integer(id), as.integer(cluster), 
                              order(id))
            if (twoclust ==1) 
                warning("an id value appears on more than one cluster")
        }
    } else {
        if (has.id) {
            # treat the id as both identifier and clustering
            clname <- levels(id)
            cluster <- as.integer(id)
            ncluster <- length(clname)
        }
        else {
            # assume an id of 1:n  This will often be competing risks
            id <- 1:n
            cluster <- 1:n
            ncluster <- n
            clname <- 1:n
        }
    }   

    status <- Y[,ncol(Y)]

    # Does everyone start in the same state?  Data may not be sorted, of course
    # If the same id is used in multiple curves, survcheck has likely already
    #  complained, but watch out for that as well.
    if (ny ==2) { # ny==2 and samestate== FALSE is very odd, but possible
        samestate <- all(istate == istate[1]) # one row per id
        stemp <- istate[1]  # the common state
    } else {
        indx <- order(X, id, Y[,1])
        first <- !duplicated(data.frame(X[indx], id[indx]))
        itemp <- (istate[indx])[first]  # first state for each id/curve pair
        samestate <- (all(itemp== itemp[1]))
        stemp <- itemp[1]
    }

    if (is.null(p0) & samestate) {
        p0 <- rep(0., nstate)
        p0[stemp] <- 1
    } 
    # I temporarily had a second check, a warning if samestate and
    #  p0 has a 1 but in a different position.  Likely a user error, but
    #  not guarranteed to be so.

    # Decide on t0 = "time 0" = starting time for all the curves
    if (!missing(start.time)) t0 <- start.time
    else if (ny ==2) 
        t0 <- min(c(0, Y[,2])) # negative times are odd but valid
    else if (samestate)  t0 <- min(Y[,1])
    else { 
        idmin <- tapply(Y[,1], id, min)
        if (isTRUE(all.equal(idmin, rep(idmin[1], length(idmin))))) 
            t0 <- idmin[1]  # everyone starts at the same time
        else {
            t0 <- min(Y[Y[,3]!=0, 2])  # first event time
            # This is the hardest case, often curves on an age scale.
            # There is not a compelling default.  You want to delay t0 as long 
            #  as possible so as to start with a reasonable number at risk, 
            #  but not remove any events. p0 will often be different for 
            #  each curve.
            # A user specified (start.time, p0) pair is best; but que sera sera
            cmax <- tapply(Y[,2], X, max) # max end time in each curve
            cmin <- tapply(Y[,1], X, min) # min starting time in each
            msg <- paste("no obs overlap the default start time of", t0,
                         "in at least one curve; specify a start.time")
            if (any(cmax < t0 | cmin >= t0)) stop(msg)
        }
    }

    status <- Y[,ncol(Y)]
    ncurve <- length(levels(X))

    if (ncol(Y) ==3) {
        sort1 <- order(X, Y[,1])
        sort2 <- order(X, Y[,2])
        position <- survflag(Y, id, X) 
    } else {
        # 2 column Y + multilple states = competing risks
        # a less common case, add a dummy 'start' column to Y
        entry <- FALSE  # option makes no sense
        minY <- min(Y[,1])

        if (minY >0) Y <- cbind(0, Y)
        else Y <- cbind(2*minY -1, Y) # intervals have to have width > 0
        sort1 <- order(X, Y[,1])
        sort2 <- order(X, Y[,2])
        position <- rep(3L, nrow(Y)) #every obs is the first and last of subject
    }
  
    # Set up indexing for the C routine
    #  nhaz = number of unique transitions that were observed
    #  hindx[i,j] = the transtion number, 0 if none of that type
    #  trmat = 2 columns containing from:to for the observed transtions
    # remember the C uses 0 based subscripts
    # mcheck returns a transition matrix which has the states as rows, but 
    #  not all states need appear as a column. Use it to create a full matrix
    #  'temp', and from it our indices
    temp <- matrix(0L, nstate, nstate)
    indx1 <- match(rownames(mcheck$transitions), states, nomatch= 0)
    indx2 <- match(colnames(mcheck$transitions), states, nomatch=0)
    temp[indx1,indx2] <- mcheck$transitions[indx1>0 ,indx2>0]
    nhaz <- sum(temp>0)  # all the unique transtions that occur
    hindx <- matrix(0L, nstate, nstate)
    hindx[temp>0] <- 1:nhaz
    trmat <- cbind(from= row(temp)[temp>0], to= col(temp)[temp>0])
    separm <- 2L*influence + 1L*(se.fit)  # for the C routine 1= se, 3=both

    # The C routine is called once per curve (X), with sort1/sort2 pointing to
    #  the relevant subset of points.
    curves <- vector("list", ncurve)
    names(curves) <- levels(X)
    n.per.curve <- as.numeric(table(X))
    n2 <- cumsum(c(0, n.per.curve))
    iX <- as.integer(X)

    for (i in 1:ncurve) {
        indx <- seq(n2[i]+1, n2[i+1])  # the relevant rows of sort1 and sort2
        indx2 <- (iX ==i)              # rows of Y, weight, etc
        
        # utime = set of time points to be reported in the output
        # because we report se(AUC), the utime vector needs to have the
        # starting time
        if (entry) {
            # There is no need to list a time point where nothing happened
            # e.g., a survSplit cutpoint that doesn't appear in the raw data.
            # So ignore rows with neither entry or exit, i.e., position==0
            # which don't end in an event
            ignore <- (position==0 & Y[,3]== 0)
            utime <- unique(sort(Y[indx2 & !ignore, 1:2])) #not sort(unique( !
        } else {
            # count only ending times: position =2 or 3
            ignore <- position <2 & Y[,3]==0
            utime <- unique(sort(Y[indx2 & !ignore, 2]))
        } 

        if (time0) utime <- unique(c(t0, utime[utime > t0]))
        else utime <- utime[utime >= t0]

        # Only compute influence for those who actually appear in this curve.
        # lablel them as 0,1,2, ... in the same order as they appear in the data
        #  (that is what the C routine wants)
        # In c2, we don't care what value is assigned to those not in this curve
        uclust <- unique(cluster[indx2])
        c2 <- match(cluster, uclust, nomatch= 1L)- 1L
        nclust <- length(uclust)
        n.id <- length(unique(id[indx2]))
        
        # If p0 is NULL we need to estimate it per curve
        #  Use the distribution of states at time t0, or if t0 is NULL
        #  the distribution just before the first event.
        #  Also create U0 = influence matrix for p0, which should be weighted
        U0 <- matrix(0, length(indx), nstate)
        sd0 <- NULL
        if (is.null(p0)) {  # a per-curve p0
            if (ny==2) atrisk <- indx2
            else atrisk <- (indx2 & Y[,1] < t0 & Y[,2] >= t0)
            if (!any(atrisk)) 
                stop(paste("no one at risk for one of the curves, at the",
                           "default time 0; specify a start.time"))
            # compute p0 and i0
            wtsum <- sum(weights[atrisk])  # weights at that time
            p00 <- tapply(weights[atrisk], istate[atrisk], sum) / wtsum
            p00 <- ifelse(is.na(p00), 0, p00)  #if a state has no one, tapply =NA
            if (all(p00 <1) & separm>0) {  # compute intitial leverage
                # To deal properly with clustering and case weights, first
                #  get the per obs influence, then weight it, then collapse
                # The C code returns clustered weighted influence
                utemp <- matrix(0, length(indx), nstate)
                for (j in 1:nstate) {
                    temp <- ifelse(istate[indx]==states[j], 1, 0)- p00[j]
                    utemp[,j] <- ifelse(atrisk[indx], weights[indx]*temp,0)/ wtsum
                }   
                U0 <- rowsum(utemp, c2[indx2]) # collapse
                sd0 <- sqrt(colSums(U0^2))
            } 
        } else p00 <- p0

        fit <- .Call(Csurvfitaj, Y, sort1[indx]-1L, sort2[indx]-1L, 
                     utime, as.integer(istate) -1L, weights, c2, nclust, p00,
                     U0, separm, entry, position, hindx-1L, trmat- 1L, t0)
        curves[[i]] <- c(list(n= length(indx), time=utime, uclust= uclust, 
                            n.id = n.id, p0=p00), fit)
        if (!is.null(sd0)) curves[[i]]$sd0 <- sd0

        if (influence >0) {
            # the influence is returned as a matrix (easier for the C code),
            #  with nclust * nstate collapsed as rows, one col for each
            #  unique event time.
            # the output is expected to be an array with id, time, state (that
            #  is how it was set up originally)
            temp <- array(curves[[i]]$influence, 
                          dim=c(nclust, nstate, length(utime)))
            temp <- aperm(temp, c(1,3,2))
            dimnames(temp) <- list(uclust, NULL, state=states)
            curves[[i]]$influence <- temp
            dimnames(curves[[i]]$influence) <- list(uclust, NULL, states)
            if (any(U0 != 0)) curves[[i]]$i0 <- U0
        }
    }

    # Names for the cumulative hazard
    cname <- paste(trmat[,1], trmat[,2], sep='.')

    # Turn the result into a survfit type object
    # The C routine returns weighted and unweighted versions of the counts:
    #  grabit1 returns the first half, grabit2 the second (uweighted), and
    #  grabit retuns all.
    grabit <- function(clist, element) {
        if (length(clist)==1) return(clist[[1]][[element]])
        temp <-(clist[[1]][[element]])
        if (ncurve ==1) return(temp)
        if (is.matrix(temp)) {
            do.call("rbind", lapply(clist, function(x) x[[element]]))
        } else {
            xx <- as.vector(unlist(lapply(clist, function(x) x[element])))
            if (inherits(temp, "table")) 
                    matrix(xx, byrow=T, ncol=length(temp))
            else xx
        }
    }
    grabit1 <- function(clist, element) {
        temp <-(clist[[1]][[element]]) 
        k <- 1:(ncol(temp)/2)
        if (length(clist) == 1) return(temp[,k, drop=FALSE])
        if (is.matrix(temp)) {
            do.call("rbind", lapply(clist, 
                                    function(x) (x[[element]])[,k, drop=FALSE]))
            }
        else {
            xx <- as.vector(unlist(lapply(clist, function(x) x[element])))
            if (inherits(temp, "table")) matrix(xx, byrow=T, ncol=length(temp))
            else xx
        }
    }
    grabit2 <- function(clist, element) {
        temp <-(clist[[1]][[element]]) 
        k <- seq(1L+ ncol(temp)/2, ncol(temp))
        if (length(clist) == 1) return(temp[,k, drop=FALSE])
        if (is.matrix(temp)) {
            do.call("rbind", lapply(clist, 
                                    function(x) (x[[element]])[,k, drop=FALSE]))
        }
        else {
            xx <- as.vector(unlist(lapply(clist, function(x) x[element])))
            if (inherits(temp, "table")) matrix(xx, byrow=T, ncol=length(temp))
            else xx
        }
    }


    kfit <- list(n =  grabit(curves, "n"),
                 time =   grabit(curves, "time"),
                 n.risk=  grabit1(curves, "n.risk"),
                 n.event= grabit(curves, "n.event"), # this has nstate cols
                 n.censor=grabit1(curves, "n.censor"),
                 pstate = grabit(curves, "pstate"),
                 n.transition = grabit1(curves, "n.transition"),
                 n.id   = grabit(curves, "n.id"),
                 cumhaz = grabit(curves, "cumhaz"))
    colnames(kfit$n.risk) <- states
    colnames(kfit$n.event) <- states
    colnames(kfit$n.censor) <- states
    colnames(kfit$pstate) <- states
    colnames(kfit$cumhaz) <- cname
    colnames(kfit$n.transition) <- cname
    
    if (entry) {
        kfit$n.enter <- grabit1(curves, "n.enter")
        colnames(kfit$n.enter) <- states
    }
    if (!isTRUE(all.equal(weights, rep(1.0, length(weights))))) {
        # also return the unweighted counts
        if (entry) {
            counts <- cbind(grabit2(curves, "n.risk"),
                            grabit2(curves, "n.transition"),
                            grabit2(curves, "n.censor"),
                            grabit2(curves, "n.enter"))
            colnames(counts) <- c(paste0("nrisk:", 1:nstate),
                                  paste0("ntrans:", cname),
                                  paste0("ncensor:", 1:nstate),
                                  paste0("nenter:", 1:nstate))
        } else {
            counts <- cbind(grabit2(curves, "n.risk"),
                            grabit2(curves, "n.transition"),
                            grabit2(curves, "n.censor"))
            colnames(counts) <- c(paste0("nrisk:", 1:nstate),
                                  paste0("nevent:", cname),
                                  paste0("ncensor:", 1:nstate))
        }
        kfit$counts <- counts
    }

    # add p0
    temp <- grabit(curves, "p0")
    if (length(curves) >1 ) 
        temp <- matrix(temp, ncol=nstate, byrow=TRUE, 
                       dimnames = list(names(curves), states))
    else names(temp) <- states
    kfit$p0 <- temp

        
    if (length(curves) > 1) {
        kfit$strata <- unlist(lapply(curves, function(x)
                         if (is.null(x$time)) 0L else length(x$time)))
    }
    if (se.fit) {
        kfit$std.err <- grabit(curves, "std.err")
        kfit$std.chaz<- grabit(curves, "std.chaz")
        kfit$std.auc <- grabit(curves, "std.auc")
        kfit$logse <- FALSE
        if (!time0 && !is.null(curves[[1]]$sd0)) 
            kfit$se0 <- grabit(curves, "sd0")
    }

    if (influence) {
        if (ncurve ==1) kfit$influence.pstate <- curves[[1]]$influence
        else kfit$influence.pstate <- lapply(curves, function(x)
             x$influence)
        if (!time0 && !is.null(curves[[i]]$i0)) {
            # this is uncommon
            if (ncurve==1) kfit$i0 <- curves[[1]]$i0
            else kfit$i0 <- lapply(curves, function(x) x$i0)
        }
    }
                              
    if (!missing(start.time)) kfit$start.time <- start.time
    kfit$transitions <- mcheck$transitions

    #       
    # Add in the confidence bands:
    #  
    if (se.fit && conf.type != "none") {
        ci <- survfit_confint(kfit$pstate, kfit$std.err, logse=FALSE, 
                                  conf.type, conf.int)
        kfit <- c(kfit, ci, conf.type=conf.type, conf.int=conf.int)
    }
    #
    kfit$states <- states
    kfit$type   <- type
    kfit$t0 <- t0
    kfit
}
