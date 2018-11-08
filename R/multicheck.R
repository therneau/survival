#
# The multi-state routines need to check their input data
#  y = 3 column survival object
#  id = subject identifier
#  istate = starting state for each row, this can be missing.
# 

multicheck <- function(y, id, istate=NULL, nerror=6) {
    n <- length(id)
    # Check 0: if y has only 2 colums there isn't much to do
    if (ncol(y)==2) {
        if (any(duplicated(id))) 
            stop("multiple observations per subject requires (start, end) times")
        istate <- rep(0L, n)
        return(list(istate=istate, transitions=table(istate, y[,2])))
    }

    # first check: no one is two places at once
    #  sort by stop time within subject, start time as the tie breaker
    index <- order(id, y[,2], y[,1])  # by start time, within stop time
    indx1 <- index[-1]
    indx2 <- index[-n]
    gap <- sign(y[indx1,1] - y[indx2,2])
    oldid <- duplicated(id[index])[-1]

    # gap = 0   (0, 10], (10, 15]
    # gap = 1   (0, 10], (12, 15]  # a hole in the time
    # gap = -1  (0, 10], (9, 15]   # overlapping times
    overlap <- (gap <0) & oldid
    if (any(overlap) & nerror>0) {
        temp <- which(overlap) 
        if (length(temp) > nerror) temp <- temp[1:nerror]
        stop("observation(s) that overlap in time, id=", (id[indx2])[temp],
             " rows = ", indx2[temp])
    }
    # If no one has more than one obs our work is done
    if (all(!oldid)) {
        istate <- rep(0L, n)
        return(list(istate=istate, transitions = table(istate, y[,3])))
    }

    # check 2: if istate is present, someone can only "jump" states if they
    #  have a gap.  If status is '2' at time 10 (a change to state 2), then
    #  for a next obs that starts at 10 it must have an istate of 2.  
    # If an obs has status=0 that means "nothing happened" and we
    #  retain the prior state.  The C routine returns a "carry forward" state
    #  vector, which starts out as istate for each new subject (or 0 if the
    #  istate vector is NULL) and marches forward.
    id.int <- match(id, unique(id))
    if (length(istate) ==0) { # construct an istate
        if (any(oldid & gap==1)) 
            warning("data has gaps and istate was not specified")
        pstate <- .Call(Cmulticheck, y, id.int, integer(n), index-1L) 
        istate <- pstate  # this gets returned
        jumps <- NULL
    }
    else {
        if (length(istate) != n) stop("wrong length for istate")
        pstate <- .Call(Cmulticheck, y, id.int, istate, index-1L) 
        bad <- ((pstate != istate)[-1] & gap==0 & oldid)
        if (any(bad) & nerror>0) {
            temp <- which(bad)
            if (length(temp) > nerror) temp <- temp[1:nerror]
            stop("observations with inconsistent istate ", index[temp],
                 " id= ", (id[index])[temp])
        }
        # jumps are unobserved transitions.  A subject changed to state 1 at 
        #  time 10, next obs is at time 20 in state 2 for instance.  
        jumps <- which((pstate != istate)[-1] & gap==1 & oldid)
    }   
                           
    # create the table of transitions
    keep <- y[,3]>0
    trtable <- table(istate[keep], y[keep,3], useNA= "ifany")
    
    rval <- list(istate=istate, transitions =trtable)
    if (length(jumps)) rval$jumps <- table(pstate[jumps], istate[jumps])
    if (any(overlap)) rval$overlap <- index[which(overlap)]
    if (any(bad))     rval$badstate <- index[which(bad)]
    rval
}
    
