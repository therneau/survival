#
# The multi-state routines need to check their input data
#  y = 3 column survival object
#  id = subject identifier
#  istate = starting state for each row, this can be missing.
# 
multicheck <- function(y, id, istate=NULL, dummy="( )") {
    n <- length(id)
    # this next line is a debug for me, since multicheck is not visible
    #  to users
    if (!is.Surv(y) || is.null(attr(y, "states")) ||
        any(y[,ncol(y)] > length(attr(y, "states"))))
        stop("multicheck called with an invalid y argument")
    to.names <- c("(censored)", attr(y, "states"))
 
    if (length(istate)==0) {
        inull<- TRUE
        istate <- factor(rep(dummy, n))
    }
    else {
        if (length(istate) !=n) stop ("wrong length for istate")
        istate <- as.factor(istate)
        inull <- FALSE
    }

    flags <- c(overlap=0, gaps=0, inconsistent=0, jumps=0)
    # Check 0: if y has only 2 colums there isn't much to do
    if (ncol(y)==2) {
        flags["overlap"] <- sum(duplicated(id))
        overlap <- list(row=duplicated(id), id=unique(id[duplicated(id)]))
        transitions <- table(istate, y[,2])
        dimnames(transitions)[[2]] <- to.names
        rval <- list(istate=istate, transitions=transitions,
                    states=c(levels(istate), to.names[-1]), flags=flags)
        if (length(overlap)) rval$overlap <- overlap
        return(rval)
    }
    
    # first check: no one is two places at once
    #  sort by stop time within subject, start time as the tie breaker
    index <- order(id, y[,2], y[,1])  # by start time, within stop time
    indx1 <- index[-1]
    indx2 <- index[-n]
    gap <- sign(y[indx1,1] - y[indx2,2])
    oldid <- duplicated(id[index])[-1]
    
    # If no one has more than one obs our work is done: overlap and gaps are
    #  impossible
    to.names <- c("", attr(y, "states"))
    if (all(!oldid)) {
        transitions <- table(istate, y[,2])
        dimnames(transitions)[[2]] <- to.names
        return(list(istate=istate, transitions=transitions,
                    states=c(levels(istate), to.names[-1]), flags=flags))
    }

    # gap = 0   (0, 10], (10, 15]
    # gap = 1   (0, 10], (12, 15]  # a hole in the time
    # gap = -1  (0, 10], (9, 15]   # overlapping times
    overlap <- which(gap <0) & oldid
    if (length(overlap)) {
        flags["overlap"] <- length(overlap)
        overlap <- list(id= (id[indx2])[overlap],
                        row = indx2[overlap])
    } 

        
    # check 2: if istate is present, someone can only "jump" states if they
    #  have a gap.  If status is '2' at time 10 (a change to state 2), then
    #  for a next obs that starts at 10 it must have an istate of 2.  
    # If an obs has status=0 that means "nothing happened" and we
    #  retain the prior state.  The C routine returns a "carry forward" state
    #  vector, which starts out as istate for each new subject (or 0 if the
    #  istate vector is NULL) and marches forward.
    id.int <- match(id, unique(id))
    if (inull) { # construct an istate
        gaps <- which(oldid & gap==1)
        if (length(gaps)) {
            # data has gaps but no istate  
            error["gaps"] <- length(gaps) # data has gaps but no istate
            gaps <- list(id= (id[indx2])[gaps],
                        row = indx2[gaps])
        }
        pstate <- .Call(Cmulticheck, y, id.int, integer(n), index-1L) 
        istate <- factor(pstate, 0:max(pstate), c("", to.names[-1])) 
        states <- c("", to.names[-1])
        jumps <- NULL
    }
    else {
        pstate <- .Call(Cmulticheck, y, id.int, as.integer(istate), index-1L)
        # the result is a mixed factor, give it proper levels
        itemp <- length(levels(istate))
        ptemp <- ifelse(c(FALSE, oldid), pstate, pstate + itemp)
        states <- unique(c(levels(istate), to.names[-1]))
        pstate <- factor(pstate, seq_len(itemp + length(to.names) -1), states)
                                         
        bad <- ((pstate != istate)[-1] & gap==0 & oldid)
        if (any(bad)) {
            temp <- which(bad)
            flags["inconsistent"] <- length(temp)
            inconsistent <- list(id=  (id[index])[temp], row= index[temp])
        }
        # jumps are unobserved transitions.  A subject changed to state 1 at 
        #  time 10, next obs is at time 20 in state 2 for instance.  
        jumps <- which((pstate != istate)[-1] & gap==1 & oldid)
        if (length(jumps)) {
            flags["jumps"] <- length(jumps)
            jumps <- list(id = (id[index])[jumps], row= index[jumps])
        }       
    }   
                           
    # create the table of transitions
    trtable <- table(istate, y[,3], useNA= "ifany")
    dimnames(trtable) <- list(from = levels(istate), to= to.names)   
    
    rval <- list(istate=istate, transitions =trtable, states=states, 
                 flags =flags)
    if (length(overlap)) rval$overlap <- overlap
    if (lengh(gaps))     rval$gaps <- gaps
    if (any(bad))     rval$inconsistent <- inconsistent
    if (length(jumps)) rval$jumps <- table(pstate[jumps], istate[jumps])
    rval
}
    
