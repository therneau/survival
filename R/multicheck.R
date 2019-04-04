multicheck <- function(formula, data, id, istate, ...) {
    Call <- match.call()
    indx <- match(c("formula", "data", "id", "istate"),
                  names(Call), nomatch=0) 
    if (indx[1] ==0) stop("A formula argument is required")
    tform <- Call[c(1,indx)]  # only keep the arguments we wanted
    tform[[1L]] <- quote(stats::model.frame)  # change the function called

    mf <- eval(tform, parent.frame())
    
    Y <- model.response(mf)
    if (!inherits(Y, "Surv")) stop("response must be a survival object")
    type <- attr(Y, "type")
    if (!(type %in% c("mright", "mcounting")))
        stop("response must be a multi-state survival object")
    n <- nrow(Y)
    
    id <- model.extract(mf, "id")
    if (is.null(id)) stop("an id argument is required")
    else if (length(id) !=n) stop("wrong length for id")
    
    istate <- model.extract(mf, "istate")
    if (!is.null(istate) && length(istate) !=n) stop("wrong length for istate")

    fit <- multicheck2(Y, id, istate)
    fit$istate <- NULL   # used by coxph, but not part of user printout
    na.action <- attr(mf, "na.action")
    if (!is.null(na.action)) {
        fit$na.action <- na.action
        # make any numbering match the input data, not the retained data
        lost <- unclass(na.action)  # the observations that were removed
        dummy <- seq.int(1, n+ length(lost))[-lost]
        for (i in c("overlap", "gap", "teleport", "jump")){
            if (!is.null(fit[[i]]$row)) {
                temp <- fit[[i]]
                temp$row <- dummy[temp$row]
                fit[[i]] <- temp
            }
        }
    }
    fit$call <- Call
    class(fit) <- "multicheck"
    fit
}

#
# The multi-state routines need to check their input data
#  y = survival object
#  id = subject identifier
#  istate = starting state for each row, this can be missing.
# The routine creates a proper current-state vector accounting for censoring
#  (which should be called cstate for 'current state', but istate is retained)
#  If a subject started in state A for instance, and had observations
#  of (0, 10, B), (10, 15, censor), (15, 20, C) the current state is
#  (A, B, B).  It checks that against the input 'istate' if it is present
#  to generate checks.
# Multiple other checks as well
#  
multicheck2 <- function(y, id, istate=NULL, dummy="()") {
    n <- length(id)
    ny <- ncol(y)
    # this next line is a debug for my code, since multicheck2 is not visible
    #  to users
    if (!is.Surv(y) || is.null(attr(y, "states")) ||
        any(y[,ncol(y)] > length(attr(y, "states"))))
        stop("multicheck2 called with an invalid y argument")
    to.names <- c(attr(y, "states"), "(censored)")
 
    if (length(istate)==0) {
        inull<- TRUE
        cstate <- factor(rep(dummy, n))
    }
    else {
        if (length(istate) !=n) stop ("wrong length for istate")
        cstate <- as.factor(istate)
        inull <- FALSE
    }

    # The vector of all state names is put in a nice printing order:
    #   initial states that are not destination states, then
    #   the destination states.  This keeps destinations in the order the
    #   user chose, while still putting initial states first.
    index <- match(levels(cstate), attr(y, "states"), nomatch=0)
    states <- c(levels(cstate)[index==0], attr(y, "states"))
    cstate2 <- factor(cstate, states)
    # we keep a form with all the levels for returning to the parent (cstate2)
    #  one without this (cstate) to make a smaller transitions table

    # initialize counts
    flag <- c(overlap=0, gap=0, teleport=0, jump=0)

    # first check: no one is two places at once
    #  sort by stop time within subject, start time as the tie breaker
    index <- order(id, y[,ny-1], y[,1])  # by start time, within stop time
    indx1 <- index[-1]
    indx2 <- index[-n]
    oldid <- duplicated(id[index])

    # if someone has an intermediate visit, i.e., (0,10, 0)(10,20,1), don't
    #  report the false 'censoring' in the transitions table
    yfac <- factor(y[,ny], c(seq(along=attr(y, "states")), 0), to.names)
    keep <- (y[,ny]!=0 | !duplicated(id[index], fromLast=TRUE))
    transitions <- table(from=cstate[keep], to= yfac[keep], useNA="ifany")

    # If no one has more than one obs our work is done: overlap and gap are
    #  impossible
    if (all(!oldid)) {
        return(list(istate=cstate2, transitions=transitions,
                    states=states, flag=flag))
    }

    # Check 0: if y has only 2 colums there isn't much to do
    if (ncol(y)==2) {
        flag["overlap"] <- sum(duplicated(id))
        overlap <- list(row= which(duplicated(id)), 
                        id=unique(id[duplicated(id)]))
        rval <- list(istate=cstate2, transitions=transitions,
                    states=states, flag=flag, overlap= overlap)
        return(rval)
    }
        
    # gap = 0   (0, 10], (10, 15]
    # gap = 1   (0, 10], (12, 15]  # a hole in the time
    # gap = -1  (0, 10], (9, 15]   # overlapping times
    gap <- sign(y[indx1,1] - y[indx2,2])
    overlap <- which(c(FALSE, gap <0) & oldid)
    if (length(overlap)) {
        flag["overlap"] <- length(overlap)
        overlap <- list(id= unique((id[indx2])[overlap]),
                        row = indx2[overlap])
    } 

    # check 2: if istate is present, someone can only "jump" states if they
    #  have a gap.  If status is '2' at time 10 (a change to state 2), then
    #  for a next obs that starts at 10 it must have an istate of 2.  
    # If an obs has status=0 that means "nothing happened" and we
    #  retain the prior state.  The C routine returns a "carry forward" state
    #  vector, which starts out as istate for each new subject (or 0 if the
    #  istate vector is NULL) and marches forward.
    # Jumps and gaps are mutually exclusive.  gap = a hole without an istate
    #  jump = a gap where they come back in another state
    id.int <- match(id, unique(id))
    jump <- tgap <- bad <- NULL
    if (inull) { # construct an istate
        gap <- which(oldid & c(FALSE, gap==1))
        if (length(gap)) {
            # data has gap but no istate  
            flag["gap"] <- length(gap) 
            tgap <- list(id= unique((id[indx2])[gap]),
                         row = indx2[gap])
        }
        pstate <- .Call(Cmulticheck, y, id.int, integer(n), index-1L) 
        cstate2 <- factor(pstate, seq(along=states)-1, states) 
        transitions <- table(from=cstate2[keep,drop=TRUE], to= yfac[keep], 
                             useNA="ifany")
   }
    else {
        ptemp <- .Call(Cmulticheck, y, id.int, -as.integer(cstate), index-1L)
        # the result is a mixed bag, using 1,2,.., length(y$state) for the
        #  responses and -1, -2, etc for starting states.
        #  subsequent ones.  Normalize it.
        temp1 <- match(to.names[pmax(1, ptemp)], states)
        temp2 <- match(levels(cstate)[pmax(1, -ptemp)], states)
        pstate<- factor(ifelse(ptemp>0 , temp1, temp2), 
                        seq(along=states), states)
                                           
        bad <- ((pstate != cstate2) & oldid & c(FALSE, gap==0))
        if (any(bad)) {
            temp <- which(bad)
            flag["teleport"] <- length(temp)
            teleport <- list(id= unique((id[index])[temp]), row= index[temp])
        }
        # jump are unobserved transitions.  A subject changed to state 1 at 
        #  time 10, next obs is at time 20 in state 2 for instance.  
        jump <- which((pstate != cstate2) & oldid & c(FALSE, gap==1))
        if (length(jump)) {
            flag["jump"] <- length(jump)
            jump <- list(id = unique((id[index])[jump]), row= index[jump])
        }       
    }   

    # we return the current state that was handed to us, and let the routine
    #  complain
    rval <- list(istate=cstate2, transitions =transitions, states=states, 
                 flag =flag)
    if (length(overlap)) rval$overlap <- overlap
    if (length(tgap))     rval$gap <- tgap
    if (any(bad))     rval$teleport <- teleport
    if (length(jump)) rval$jump <- jump

    rval
}
    
