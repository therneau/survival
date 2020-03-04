survcheck <- function(formula, data, id, istate, istate0="(s0)", 
                      timefix=TRUE, ...) {
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
    if (!is.logical(timefix) || length(timefix) > 1)
        stop("invalid value for timefix option")
    if (timefix) Y <- aeqSurv(Y)
    
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
    # the next few line are a debug for my code, since survcheck2 is not visible
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
        if (is.factor(istate)) cstate <- istate[, drop=TRUE] #drop unused levels
        else cstate <- as.factor(istate)
        inull <- FALSE
    }

    ystate <- attr(y, "states")
    # The vector of all state names is put in a nice printing order:
    #   initial states that are not destination states, then
    #   the destination states.  This keeps destinations in the order the
    #   user chose, while still putting initial states first.
    index <- match(levels(cstate), ystate, nomatch=0)
    states <- c(levels(cstate)[index==0], ystate)
    cstate2 <- factor(cstate, states)
    # we keep a form with all the levels for returning to the parent (cstate2)
    #  one without this (cstate) to make a smaller transitions table

    # Calculate the counts per id for each state, e.g., 10 subjects had
    #  3 visits to state 2, etc.  
    # Don't count "censored" as an endpoint, nor any missings.  But we can't
    #  just omit censored rows, or those with 0 transitions don't get counted!
    #  Instead remove the 'censored' column after making tab1
    tab1 <- table(id, factor(y[,ncol(y)], 0:length(ystate)))[,-1, drop=FALSE]
    tab1 <- cbind(tab1, rowSums(tab1))
    tab1.levels <- sort(unique(c(tab1)))  #unique counts
    events <- apply(tab1, 2, function(x) table(factor(x, tab1.levels)))
    dimnames(events) = list("count"= tab1.levels,
                                "state"= c(ystate, "(any)"))
    
    # check for errors
    sindx <- match(ystate, states)
    if (ncol(y)==2) y <- cbind(0,y) # make it 3 cols for the C routine
    stat2 <- ifelse(y[,3]==0, 0L, 
                    sindx[pmax(1, y[,3])])  # map the status
    index <- order(id, y[,ny-1], y[,1])  # by start time, within stop time
    id2 <- match(id, id)  # we need unique integers, but don't need 1, 2,...
    check <- .Call(Cmulticheck, y, stat2, id2, as.integer(cstate2), index-1L)
    if (inull) {
        # if there was no istate entered in, use the constructed one from
        # the check routine
        cstate2 <-factor(check$cstate, seq(along=states), states)
    }       

    # create the transtions table
    # if someone has an intermediate visit, i.e., (0,10, 0)(10,20,1), don't
    #  report the false 'censoring' in the transitions table
    yfac <- factor(y[,3], c(seq(along=ystate), 0), to.names)
    keep <- (y[index,3]!=0 | !duplicated(id[index], fromLast=TRUE))
    keep[index] <- keep  # put it back into data order
    transitions <- table(from=cstate2[keep], 
                         to= yfac[keep], 
                         useNA="ifany")

    # now continue with error checks
    # A censoring hole in the middle, such as happens with survSplit,
    #  is correctly ignored by the C routine. 
    #
    mismatch <- (as.numeric(cstate2) != check$cstate)
                 

    # gap = 0   (0, 10], (10, 15]
    # gap = 1   (0, 10], (12, 15]  # a hole in the time
    # gap = -1  (0, 10], (9, 15]   # overlapping times
    flag <- c(overlap= sum(check$gap < 0), 
              gap =    sum(check$gap > 0 & !mismatch),
              jump =   sum(check$gap > 0 & mismatch),
              teleport = sum(check$gap==0 & mismatch & check$dupid==1))

    rval <- list(states=states, transitions=transitions,
                 events= t(events), flag=flag, 
                 istate= factor(check$cstate, seq(along=states), states))
 
    # add error details, if necessary
    if (flag["overlap"] > 0) {
        j <- which(check$gap < 0)
        rval$overlap <- list(row=j, id= unique(id[j]))
    }       
    if (flag["gap"] > 0) {
        j <- which(check$gap > 0 & !mismatch)
        rval$gap <- list(row=j, id= unique(id[j]))
    }
    if (flag["jump"] > 0) {
        j <- which(check$gap > 0 & mismatch)
        rval$jump <- list(row=j, id= unique(id[j]))
    }
    if (flag["teleport"] > 0) {
        j <- (check$gap==0 & mismatch)
        rval$teleport <- list(row=j, id= unique(id[j]))
    }

    rval
}
    
