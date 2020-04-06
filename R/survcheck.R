survcheck <- function(formula, data, subset, na.action,  id, istate, 
                      istate0="(s0)", timefix=TRUE, ...) {
    Call <- match.call()
    indx <- match(c("formula", "data", "id", "istate", "subset", "na.action"),
                  names(Call), nomatch=0) 
    if (indx[1] ==0) stop("A formula argument is required")
    tform <- Call[c(1,indx)]  # only keep the arguments we wanted
    tform[[1L]] <- quote(stats::model.frame)  # change the function called

    mf <- eval(tform, parent.frame())
    Terms <- terms(mf)

    Y <- model.response(mf)
    isSurv2 <- inherits(Y, "Surv2")
    if (isSurv2) {
        # this is Surv2 style data
        # if there were any obs removed due to missing, remake the model frame
        if (length(attr(mf, "na.action"))) {
            tform$na.action <- na.pass
            mf <- eval.parent(tform)
        }
        if (!is.null(attr(Terms, "specials")$cluster))
            stop("cluster() cannot appear in the model statement")
        new2 <- surv2data(mf, check=TRUE)
        mf <- new2$mf
        istate <- new2$istate
        id <- new2$id
        Y <- new2$y
        Ydup <- which(Y[,1] == Y[,2])
        iddup <- id[Ydup]
        if (anyNA(mf[-1]) || length(Ydup)) { 
            if (missing(na.action)) temp <- get(getOption("na.action"))(mf[-1])
            else temp <- na.action(mf[-1])
            omit <- c(Ydup, unclass(attr(temp, "na.action")))
            mf <- mf[-omit,]
            Y <- Y[-omit]
            id <- id[-omit]
            istate <- istate[-omit]
        }  else omit <- Ydup                   
        n <- nrow(mf)
    }       
    else {
        if (!is.Surv(Y)) stop("Response must be a survival object")
        id <- model.extract(mf, "id")
        istate <- model.extract(mf, "istate")
        omit <- attr(mf, "na.action")
        n <- nrow(Y)
    }
    if (n==0) stop("No (non-missing) observations")

    type <- attr(Y, "type")
    if (type=="right") Y <- Surv(Y[,1], factor(Y[,2]))  # pretend its multi
    else if (type=="counting")  Y <- Surv(Y[,1], Y[,2], factor(Y[,3]))
    else if (!(type %in% c("mright", "mcounting")))
        stop("response must be right censored")
    n <- nrow(Y)
    if (!is.logical(timefix) || length(timefix) > 1)
        stop("invalid value for timefix option")
    if (timefix) Y <- aeqSurv(Y)
    
    if (is.null(id)) stop("an id argument is required")
    else if (length(id) !=n) stop("wrong length for id")
     
    if (!is.null(istate) && length(istate) !=n) stop("wrong length for istate")

    fit <- survcheck2(Y, id, istate, istate0)
    fit$n <- c(id = length(unique(id)), observations =length(id), 
               transitions = sum(fit$transitions))

    fit$flag <- c(fit$flag, "duplicate"=0)
    if (isSurv2) {
        # make any numbering match the input data, not the retained data
        toss1 <- new2$isort[new2$last]  # original obs numbers that were the
                                       # last for a subject
        dummy <- seq(along=new2$isort)[-toss1]  # rows we kept
        if (length(Ydup)) {
            fit$flag["duplicate"] <- length(Ydup)
            # if rows i and i+1 are duplicate times, we see it as i, the
            #  R duplicated function as i+1.  Mimic that rule.
            fit$duplicate <- list(id= unname(iddup), row=dummy[Ydup]+1)
        } 

        if (length(omit)>0) dummy <- dummy[-omit]
        for (i in c("overlap", "gap", "teleport", "jump")){
            if (!is.null(fit[[i]]$row)) {
                temp <- fit[[i]]
                temp$row <- dummy[temp$row]
                fit[[i]] <- temp
            }
        }
    }
    else if (length(omit)) {
        dummy <- seq.int(1, n+ length(omit))[-omit]
        for (i in c("overlap", "gap", "teleport", "jump")){
            if (!is.null(fit[[i]]$row)) {
                temp <- fit[[i]]
                temp$row <- dummy[temp$row]
                fit[[i]] <- temp
            }
        }
    }

    fit$Y <- Y      # used by the summary function for details
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
    # the next few line are a debug for my code; survcheck2 is not visible
    #  to users so only survival can call it directly
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

    # Calculate the counts per id for each state, e.g., 10 subjects had
    #  3 visits to state 2, etc.  
    # Count the censors, so that each subject gets a row in the table,
    #  but then toss that column
    tab1 <- table(id, factor(y[,ncol(y)], 0:length(ystate)))[,-1, drop=FALSE]
    tab1 <- cbind(tab1, rowSums(tab1))
    tab1.levels <- sort(unique(c(tab1)))  #unique counts
    events <- apply(tab1, 2, function(x) table(factor(x, tab1.levels)))
    dimnames(events) = list("count"= tab1.levels,
                                "state"= c(ystate, "(any)"))

    # Use a C routine to create 3 variables: a: an index of whether this is
    #   the first (1) or last(2) observation for a subject, 3=both, 0=neither,
    #  b. current state, and 
    #  c. sign of (start of this interval - end of prior one)
    # start by making stat2 = status re-indexed to the full set of states
    ny <- ncol(y)
    sindx <- match(ystate, states)
    stat2 <- ifelse(y[,ny]==0, 0L, sindx[pmax(1L, y[,ny])])
    id2 <- match(id, unique(id))  # we need unique integers
    if (ncol(y)==2) {
        index <- order(id, y[,1])
        check <- .Call(Cmulticheck, rep(0., n), y[,1], stat2, id2,
                       as.integer(cstate2), index- 1L)
    } else {
        index <- order(id, y[,2], y[,1])
        check <- .Call(Cmulticheck, y[,1], y[,2], stat2, id2,
                       as.integer(cstate2), index- 1L)
    }       

    if (inull && ny> 2) {
        # if there was no istate entered in, use the constructed one from
        # the check routine
        # if ny=2 then every row starts at time 0
        cstate2 <-factor(check$cstate, seq(along=states), states)
    }       

    # create the transtions table
    # if someone has an intermediate visit, i.e., (0,10, 0)(10,20,1), don't
    #  report the false 'censoring' in the transitions table
    # make it compact by removing any cols that are all 0, and rows of
    #  states that never occur (sometimes the starting state is a factor
    #  with unused levels)
    keep <- (stat2 !=0 | check$dupid > 1)  # not censored or last obs of this id
    transitions <- table(from=cstate2[keep], 
                         to= factor(stat2[keep], c(seq(along=states), 0),
                                    c(states, "(censored)")),
                         useNA="ifany")
    nr <- nrow(transitions)
    never <- (rowSums(transitions) + colSums(transitions[,1:nr]))==0
    transitions <- transitions[!never, colSums(transitions)>0, drop = FALSE]

    # now continue with error checks
    # A censoring hole in the middle, such as happens with survSplit,
    #  uses "last state carried forward" in Cmultistate, which also
    #  sets the "gap" to 0 for the first obs of a subject
    mismatch <- (as.numeric(cstate2) != check$cstate)
                 
    # gap = 0   (0, 10], (10, 15]
    # gap = 1   (0, 10], (12, 15]  # a hole in the time
    # gap = -1  (0, 10], (9, 15]   # overlapping times
    flag <- c(overlap= sum(check$gap < 0), 
              gap =    sum(check$gap > 0 & !mismatch),
              jump =   sum(check$gap > 0 & mismatch),
              teleport = sum(check$gap==0 & mismatch & check$dupid%%2 ==0))

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
        j <- which(check$gap==0 & mismatch)
        rval$teleport <- list(row=j, id= unique(id[j]))
    }

    rval
}
    
