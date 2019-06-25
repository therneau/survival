survcheck <- function(formula, data, id, istate, istate0="(s0)", ...) {
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
    if (type=="right") Y <- Surv(Y[,1], factor(Y[,2]))  # pretend its multi
    else if (type=="counting")  Y <- Surv(Y[,1], Y[,2], factor(Y[,3]))
    else if (!(type %in% c("mright", "mcounting")))
        stop("response must be right censored")
    n <- nrow(Y)
    
    id <- model.extract(mf, "id")
    if (is.null(id)) stop("an id argument is required")
    else if (length(id) !=n) stop("wrong length for id")
     
    istate <- model.extract(mf, "istate")
    if (!is.null(istate) && length(istate) !=n) stop("wrong length for istate")

    fit <- survcheck2(Y, id, istate, istate0)
#    fit$istate <- NULL   # used by coxph, but not part of user printout
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
    fit$Y <- Y      # used by the summary function
    fit$id <- unname(id)    #   ""  ""
    fit$call <- Call
    class(fit) <- "survcheck"
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
survcheck2 <- function(y, id, istate=NULL, istate0="(s0)") {
    n <- length(id)
    ny <- ncol(y)
    # this next line is a debug for my code, since survcheck2 is not visible
    #  to users
    if (!is.Surv(y) || is.null(attr(y, "states")) ||
        any(y[,ncol(y)] > length(attr(y, "states"))))
        stop("survcheck2 called with an invalid y argument")
    to.names <- c(attr(y, "states"), "(censored)")
 
    if (length(istate)==0) {
        inull<- TRUE
        cstate <- factor(rep(istate0, n))
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

    # Calculate the counts per id for each state, e.g., 10 subjects had
    #  3 visits to state 2, etc.  
    # Don't count "censored" as an endpoint, nor any missings.  But we can't
    #  just omit censored rows, or those with 0 transitions don't get counted!
    #  Instead remove the 'censored' column after making tab1
    tab1 <- table(id, factor(y[,ncol(y)], 0:length(attr(y, 'states'))))[,-1, drop=FALSE]
    tab1 <- cbind(tab1, rowSums(tab1))
    tab1.levels <- sort(unique(c(tab1)))  #unique counts
    events <- apply(tab1, 2, function(x) table(factor(x, tab1.levels)))
    dimnames(events) = list("count"= tab1.levels,
                                "state"= c(attr(y, "states"), "(any)"))
    
    # first check: no one is two places at once
    #  sort by stop time within subject, start time as the tie breaker
    index <- order(id, y[,ny-1], y[,1])  # by start time, within stop time
    indx1 <- index[-1]
    indx2 <- index[-n]
    oldid <- duplicated(id[index])

    # if someone has an intermediate visit, i.e., (0,10, 0)(10,20,1), don't
    #  report the false 'censoring' in the transitions table
    yfac <- factor(y[,ny], c(seq(along=attr(y, "states")), 0), to.names)
    keep <- (y[index,ny]!=0 | !duplicated(id[index], fromLast=TRUE))
    transitions <- table(from=cstate2[keep,drop=TRUE], 
                         to= yfac[keep, drop=TRUE], 
                         useNA="ifany")

    # If no one has more than one obs our work is done: overlap and gap are
    #  impossible
    if (all(!oldid)) {
        return(list(states=states, transitions=transitions,
                    events = t(events), flag=flag, istate=cstate2))
    }

    # Check 0: if y has only 2 colums there isn't much to do
    if (ncol(y)==2) {
        flag["overlap"] <- sum(duplicated(id))
        overlap <- list(row= which(duplicated(id)), 
                        id=unique(id[duplicated(id)]))
        rval <- list(states=states, transitions=transitions,
                     events= t(events), flag=flag, istate=cstate2,
                     overlap= overlap)
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
    id.int <- match(id, unique(id))  # for the C routine
    y3 <- factor(yfac, states)
    jump <- tgap <- bad <- NULL
    if (inull) { # construct an istate
        gap <- which(oldid & c(FALSE, gap==1))
        if (length(gap)) {
            # data has gap but no istate  
            flag["gap"] <- length(gap) 
            tgap <- list(id= unique((id[indx2])[gap]),
                         row = indx2[gap])
        }
        pstate <- .Call(Cmulticheck, as.integer(y[,ny]), id.int, 
                        integer(n), index-1L) 
        pstate <- factor(pstate, seq(along=states)-1, states) 
        # line above: pstate will be 0 for an entry, non-zero for carried
        #  forward state, and states will be (s0) followed by attr(y, states)
        # since no istate was passed in, make transitions using our
        #  constructed state
        transitions <- table(from=pstate[keep,drop=TRUE], to= yfac[keep], 
                             useNA="ifany")
   }
    else {
        # for this call we want to map y[,ny] to the states, but 0 for censored
        y3 <- ifelse(y[,ny]==0, 0L, match(to.names[pmax(1, y[,ny])], states))
        pstate <- .Call(Cmulticheck, y3, id.int, as.integer(cstate2), index-1L)
        pstate <- factor(pstate, seq(along=states), states)
                                           
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
    rval <- list(states=states, transitions =transitions, events = t(events),
                 flag = flag, istate=pstate)
    if (length(overlap)) rval$overlap <- overlap
    if (length(tgap))     rval$gap <- tgap
    if (any(bad))     rval$teleport <- teleport
    if (length(jump)) rval$jump <- jump

    rval
}
    
