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
    fit$istate <- NULL
    na.action <- attr(mf, "na.action")
    if (!is.null(na.action)) {
        # make any numbering match the input data, not the retained data
        lost <- unclass(na.action)  # the observations that were removed
        dummy <- seq.int(1, n+ length(lost))[-lost]
        for (i in c("overlap", "gap", "teleport", "jump")){
            if (!is.null(fit[[i]])$row) {
                temp <- fit[[i]]
                temp$row <- dummy[temp$row]
                fit[[i]] <- temp
            }
        }
    }
    fit
}

#
# The multi-state routines need to check their input data
#  y = survival object
#  id = subject identifier
#  istate = starting state for each row, this can be missing.
# 
multicheck2 <- function(y, id, istate=NULL, dummy="( )") {
    n <- length(id)
    # this next line is a debug for my code, since multicheck2 is not visible
    #  to users
    if (!is.Surv(y) || is.null(attr(y, "states")) ||
        any(y[,ncol(y)] > length(attr(y, "states"))))
        stop("multicheck2 called with an invalid y argument")
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

    # normalize the state names
    states <- unique(c(levels(istate), attr(y, "states")))

    # initialize counts
    flag <- c(overlap=0, gap=0, teleport=0, jump=0)

    # Check 0: if y has only 2 colums there isn't much to do
    if (ncol(y)==2) {
        flag["overlap"] <- sum(duplicated(id))
        overlap <- list(row=duplicated(id), id=unique(id[duplicated(id)]))
        transitions <- table(istate, y[,2])
        dimnames(transitions)[[2]] <- to.names
        rval <- list(istate=istate, transitions=transitions,
                    states=states, flag=flag)
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
    
    # If no one has more than one obs our work is done: overlap and gap are
    #  impossible
    to.names <- c("", attr(y, "states"))
    if (all(!oldid)) {
        transitions <- table(istate, y[,2])
        dimnames(transitions)[[2]] <- to.names
        return(list(istate=istate, transitions=transitions,
                    states=states, flag=flag))
    }

    # gap = 0   (0, 10], (10, 15]
    # gap = 1   (0, 10], (12, 15]  # a hole in the time
    # gap = -1  (0, 10], (9, 15]   # overlapping times
    overlap <- which((gap <0) & oldid)
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
    id.int <- match(id, unique(id))
    if (inull) { # construct an istate
        gap <- which(oldid & gap==1)
        if (length(gap)) {
            # data has gap but no istate  
            error["gap"] <- length(gap) 
            gap <- list(id= unique((id[indx2])[gap]),
                        row = indx2[gap])
        }
        pstate <- .Call(Cmulticheck, y, id.int, integer(n), index-1L) 
        istate <- factor(pstate, 0:max(pstate), c("", to.names[-1])) 
        states <- c("", to.names[-1])
        jump <- NULL
    }
    else {
        pstate <- .Call(Cmulticheck, y, id.int, as.integer(istate), index-1L)
        # the result is a mixed factor, give it proper levels
        itemp <- length(levels(istate))
        ptemp <- ifelse(c(FALSE, oldid), pstate, pstate + itemp)
        pstate <- factor(pstate, seq_len(itemp + length(to.names) -1), states)
                                         
        bad <- ((pstate != istate)[-1] & gap==0 & oldid)
        if (any(bad)) {
            temp <- which(bad)
            flag["teleport"] <- length(temp)
            teleport <- list(id= unique((id[index])[temp]), row= index[temp])
        }
        # jump are unobserved transitions.  A subject changed to state 1 at 
        #  time 10, next obs is at time 20 in state 2 for instance.  
        jump <- which((pstate != istate)[-1] & gap==1 & oldid)
        if (length(jump)) {
            flag["jump"] <- length(jump)
            jump <- list(id = unique((id[index])[jump]), row= index[jump])
        }       
    }   
                           
    # create the table of transitions
    trtable <- table(istate, y[,3], useNA= "ifany")
    dimnames(trtable) <- list(from = levels(istate), to= to.names)   
    
    rval <- list(istate=istate, transitions =trtable, states=states, 
                 flag =flag)
    if (length(overlap)) rval$overlap <- overlap
    if (lengh(gap))     rval$gap <- gap
    if (any(bad))     rval$teleport <- teleport
    if (length(jump)) rval$jump <- table(pstate[jump], istate[jump])
    rval
}
    
