# Automatically generated from the noweb directory
docurve2 <- function(entry, etime, status, istate, wt, states, id, 
                     se.fit, influence=FALSE, p0) {
    timeset <- sort(unique(etime))
    nstate <- length(states)
    uid <- sort(unique(id))
    index <- match(id, uid)
    # Either/both of id and cstate might be factors.  Data may not be in
    #  order.  Get the initial state for each subject
    temp1 <- order(id, entry)
    temp2 <- match(uid, id[temp1])
    cstate <- (as.numeric(istate)[temp1])[temp2]  # initial state for each

    # The influence matrix can be huge, make sure we have enough memory
    if (influence) {
        needed <- max(nstate * length(uid), 1 + length(timeset))
        if (needed > .Machine$integer.max)
            stop("number of rows for the influence matrix is > the maximum integer")
    }
    storage.mode(wt) <- "double" # just in case someone had integer weights

    # Compute p0 (unless given by the user)
    if (is.null(p0)) {
        if (all(status==0))  t0 <- max(etime)  #failsafe
        else t0 <- min(etime[status!=0])  # first transition event
        at.zero <- (entry < t0 & etime >= t0) 
        wtsum <- sum(wt[at.zero])  # weights for a subject may change
        p0 <- tapply(wt[at.zero], istate[at.zero], sum) / wtsum
        p0 <- ifelse(is.na(p0), 0, p0)  #for a state not in at.zero, tapply =NA
    }
    # initial leverage matrix
    nid <- length(uid)
    i0  <- matrix(0., nid, nstate)
    if (all(p0 <1)) {  #actually have to compute it
        who <- index[at.zero]  # this will have no duplicates
        for (j in 1:nstate) 
            i0[who,j] <- (ifelse(istate[at.zero]==states[j], 1, 0) - p0[j])/wtsum
    }
     
    storage.mode(cstate) <- "integer"
    storage.mode(status) <- "integer"
    # C code has 0 based subscripts
    if (influence) se.fit <- TRUE   # se.fit is free in this case

    fit <- .Call(Csurvfitci, c(entry, etime), 
                 order(entry) - 1L,
                 order(etime) - 1L,
                 length(timeset),
                 status,
                 as.integer(cstate) - 1L,
                 wt,
                 index -1L,
                 p0, i0,
                 as.integer(se.fit) + 2L*as.integer(influence))

    if (se.fit) 
        out <- list(n=length(etime), time= timeset, p0 = p0,
                    sp0= sqrt(colSums(i0^2)),
             pstate = fit$p, std.err=fit$std,
             n.risk = fit$nrisk,
             n.event= fit$nevent,
             n.censor=fit$ncensor,
             cumhaz = fit$cumhaz)
    else out <- list(n=length(etime), time= timeset, p0=p0,
             pstate = fit$p,
             n.risk = fit$nrisk, 
             n.event = fit$nevent, 
             n.censor= fit$ncensor, 
             cumhaz= fit$cumhaz)
    if (influence) {
        temp <-  array(fit$influence, 
                       dim=c(length(uid), nstate, 1+ length(timeset)),
                       dimnames=list(uid, NULL, NULL))
        out$influence.pstate <- aperm(temp, c(1,3,2))
    }
    out
}
survfitCI <- function(X, Y, weights, id, cluster, robust, istate,
                      stype=1, ctype=1,
                      se.fit=TRUE,
                      conf.int= .95,
                      conf.type=c('log',  'log-log',  'plain', 'none', 
                                  'logit', "arcsin"),
                      conf.lower=c('usual', 'peto', 'modified'),
                      influence = FALSE, start.time, p0, type){

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
    if (type != 1) warning("only stype=1, ctype=1 currently implimented for multi-state data")

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

 
    if (is.logical(influence)) {
        # TRUE/FALSE is treated as all or nothing
        if (!influence) influence <- 0L
        else influence <- 3L
    }
    else if (!is.numeric(influence))
        stop("influence argument must be numeric or logical")
    if (!(influence %in% 0:3)) stop("influence argument must be 0, 1, 2, or 3")
    else influence <- as.integer(influence)
 
    if (!se.fit) {
        # if the user asked for no standard error, skip any robust computation
        ncluster <- 0L
        influence <- 0L
    }

    type <- attr(Y, "type")
    # This line should be unreachable, unless they call "surfitCI" directly
    if (type !='mright' && type!='mcounting')
         stop(paste("multi-state computation doesn't support \"", type,
                          "\" survival data", sep=''))
    
    # If there is a start.time directive, start by removing any prior events
    if (!missing(start.time)) {
        if (!is.numeric(start.time) || length(start.time) !=1
            || !is.finite(start.time))
            stop("start.time must be a single numeric value")
        toss <- which(Y[,ncol(Y)-1] <= start.time)
        if (length(toss)) {
            n <- nrow(Y)
            if (length(toss)==n) stop("start.time has removed all observations")
            Y <- Y[-toss,,drop=FALSE]
            X <- X[-toss]
            weights <- weights[-toss]
            if (length(id) ==n) id <- id[-toss]
            if (!missing(istate) && length(istate)==n) istate <- istate[-toss]
            }
    }
    n <- nrow(Y)
    status <- Y[,ncol(Y)]
    ncurve <- length(levels(X))
    
    # The user can call with cluster, id, robust or any combination.
    # If only id, treat it as the cluster too
    if (missing(robust) || length(robust)==0) robust <- TRUE
    if (!robust) stop("multi-state survfit supports only a robust variance")

    has.cluster <-  !(missing(cluster) || length(cluster)==0) 
    has.id <-       !(missing(id) || length(id)==0)
    if (has.id) id <- as.factor(id)
    else  {
        if (ncol(Y) ==3) stop("an id statement is required for start,stop data")
        id <- 1:n  # older default, which could lead to invalid curves
    }
    if (influence && !(has.cluster || has.id)) {
        cluster <- seq(along=X)
        has.cluster <- TRUE
    }

    if (has.cluster) {
        if (is.factor(cluster)) {
            clname <- levels(cluster)
            cluster <- as.integer(cluster)
        } else {
            clname  <- sort(unique(cluster))
            cluster <- match(cluster, clname)
        }
        ncluster <- length(clname)
    } else {
        if (has.id) {
            # treat the id as both identifier and clustering
            clname <- levels(id)
            cluster <- as.integer(id)
            ncluster <- length(clname)
        }
        else {
            ncluster <- 0  # has neither
            clname <- NULL
        }
    }

    if (missing(istate) || is.null(istate))
        mcheck <- survcheck2(Y, id)  
    else mcheck <- survcheck2(Y, id, istate)
    if (any(mcheck$flag > 0)) stop("one or more flags are >0 in survcheck")
    states <- mcheck$states
    istate <- mcheck$istate
    nstate <- length(states) 
    smap <- c(0, match(attr(Y, "states"), states))
    Y[,ncol(Y)] <- smap[Y[,ncol(Y)] +1]      # new states may be a superset
    status <- Y[,ncol(Y)]

    if (mcheck$flag["overlap"] > 0)
        stop("a subject has overlapping time intervals")
#    if (mcheck$flag["gap"] > 0 || mcheck$flag["jump"] > 0)
#        warning("subject(s) with time gaps, results may be questionable")

    # The states of the status variable are the first columns in the output
    #  any extra initial states are later in the list. 
    # Now that we know the names, verify that p0 is correct (if present)
    if (!missing(p0) && !is.null(p0)) {
        if (length(p0) != nstate) stop("wrong length for p0")
        if (!is.numeric(p0) || abs(1-sum(p0)) > sqrt(.Machine$double.eps))
            stop("p0 must be a numeric vector that adds to 1")
    } else p0 <- NULL
    curves <- vector("list", ncurve)
    names(curves) <- levels(X)
                            
    if (ncol(Y)==2) {  # 1 transition per subject
        # dummy entry time that is < any event time
        t0 <- min(0, Y[,1])
        entry <- rep(t0-1, nrow(Y))
        for (i in levels(X)) {
            indx <- which(X==i)
            curves[[i]] <- docurve2(entry[indx], Y[indx,1], status[indx], 
                                    istate[indx], weights[indx], 
                                    states, 
                                    id[indx], se.fit, influence, p0)
         }
    }
    else {
        # extra censors
        indx <- order(id, Y[,2])   # in stop order
        extra <- (survflag(Y[indx,], id[indx]) ==0 & (Y[indx,3] ==0))
        # If a subject had obs of (a, b)(b,c)(c,d), and c was a censoring
        #  time, that is an "extra" censoring/entry at c that we don't want
        #  to count.  Deal with it by changing that subject
        #  to (a,b)(b,d).  Won't change S(t), only the n.censored/n.enter count.
        if (any(extra)) {
            e2 <- indx[extra]
            Y <- cbind(Y[-(1+e2),1], Y[-e2,-1])
            status <- status[-e2]
            X <- X[-e2]
            id <- id[-e2]
            istate <- istate[-e2]
            weights <- weights[-e2]
            indx <- order(id, Y[,2])
        }
        # Now to work
        for (i in levels(X)) {
            indx <- which(X==i)
        #    temp <- docurve1(Y[indx,1], Y[indx,2], status[indx], 
        #                          istate[indx], weights[indx], states, id[indx])
            curves[[i]] <- docurve2(Y[indx,1], Y[indx,2], status[indx], 
                                    istate[indx],
                                    weights[indx], states, id[indx], se.fit, 
                                    influence, p0)
        }
    }

    # Turn the result into a survfit type object
    grabit <- function(clist, element) {
        temp <-(clist[[1]][[element]]) 
        if (is.matrix(temp)) {
            do.call("rbind", lapply(clist, function(x) x[[element]]))
            }
        else {
            xx <- as.vector(unlist(lapply(clist, function(x) x[element])))
            if (inherits(temp, "table")) matrix(xx, byrow=T, ncol=length(temp))
            else xx
        }
    }

    # we want to rearrange the cumulative hazard to be in time order
    #   with one column for each observed transtion.  
    nstate <- length(states)
    temp <- matrix(0, nstate, nstate)
    indx1 <- match(rownames(mcheck$transitions), states)
    indx2 <- match(colnames(mcheck$transitions), states, nomatch=0) #ignore censor
    temp[indx1, indx2[indx2>0]] <- mcheck$transitions[,indx2>0]
    ckeep <- which(temp>0)
    names(ckeep) <- outer(1:nstate, 1:nstate, paste, sep='.')[ckeep]
    #browser()

    if (length(curves) ==1) {
        keep <- c("n", "time", "n.risk", "n.event", "n.censor", "pstate",
                  "p0", "cumhaz", "influence.pstate")
        if (se.fit) keep <- c(keep, "std.err", "sp0")
        kfit <- (curves[[1]])[match(keep, names(curves[[1]]), nomatch=0)]
        names(kfit$p0) <- states
        if (se.fit) kfit$logse <- FALSE
        kfit$cumhaz <- t(kfit$cumhaz[ckeep,,drop=FALSE])
        colnames(kfit$cumhaz) <- names(ckeep)
    }
    else {    
        kfit <- list(n =      as.vector(table(X)),  #give it labels
                     time =   grabit(curves, "time"),
                     n.risk=  grabit(curves, "n.risk"),
                     n.event= grabit(curves, "n.event"),
                     n.censor=grabit(curves, "n.censor"),
                     pstate = grabit(curves, "pstate"),
                     p0     = grabit(curves, "p0"),
                     strata= unlist(lapply(curves, function(x)
                         if (is.null(x$time)) 0L else length(x$time))))
        kfit$p0 <- matrix(kfit$p0, ncol=nstate, byrow=TRUE,
                          dimnames=list(names(curves), states))
        if (se.fit) {
            kfit$std.err <- grabit(curves, "std.err")
            kfit$sp0<- matrix(grabit(curves, "sp0"),
                              ncol=nstate, byrow=TRUE)
            kfit$logse <- FALSE
        }

        # rearrange the cumulative hazard to be in time order, with columns
        #  for each transition
        kfit$cumhaz <- do.call(rbind, lapply(curves, function(x)
            t(x$cumhaz[ckeep,,drop=FALSE])))
        colnames(kfit$cumhaz) <- names(ckeep)
     
        if (influence) kfit$influence.pstate <- 
            lapply(curves, function(x) x$influence.pstate)
    }                         

    if (!missing(start.time)) kfit$start.time <- start.time
    kfit$transitions <- mcheck$transitions

    #       
    # Last bit: add in the confidence bands:
    #  
    if (se.fit) {
        kfit$conf.int <- conf.int
        kfit$conf.type <- conf.type
        if (conf.type != "none") {
            ci <- survfit_confint(kfit$pstate, kfit$std.err, logse=FALSE, 
                                  conf.type, conf.int)
            kfit <- c(kfit, ci, conf.type=conf.type, conf.int=conf.int)
        }
    }
    kfit$states <- states
    kfit$type   <- attr(Y, "type")
    kfit
}
