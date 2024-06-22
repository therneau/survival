#
# Functions to go back and forth from timeline form
#
totimeline <- function(formula, data, id, istate) {
    if (missing(formula) || missing(data) || missing(id))
        stop("formula, data, and id arguments are required")
    Call <- match.call()
    tcall <- Call
    tcall[[1L]] <- quote(stats::model.frame)

    mf <- eval(tcall, parent.frame())

    Y <- model.response(mf)
    id <- model.extract(mf, "id")

    if (!inherits(Y, "Surv")) stop("response must be a Surv object")
    type <- attr(Y, "type")
    if (type== "left" || type=="interval")
        stop("not valid for interval censored or left censored data")
    if (ncol(Y) != 3) stop("initial data is not in (time1, time2) form")
    
    # get a name for the resulting time and state variables, by parsing
    #  the formula.  Don't try too hard, though.
    tname <- "(time)"
    sname <- "(state)"  # backup defaults
    lhs <- formula[[2]]
    if ((is.name(lhs[[1]]) && lhs[[1]]== as.name("Surv")) ||
        deparse(lhs[[1]]) == "survival::Surv")  {
        # the lhs if of length 4
        if (is.name(lhs[[3]])) tname <- as.character(lhs[[3]])
        if (is.name(lhs[[4]]))  sname <- as.character(lhs[[4]])
        else if (is.call(lhs[[4]])){
            temp <- lhs[[4]]
            if (deparse(temp[[1]])== "factor" && is.name(temp[[2]]))
                sname <- as.character(temp[[2]])
        }
    }   
    
    if (is.name(Call$id)) {
        idname <- as.character(Call$id)
        i <- match(idname, names(mf))
        if (is.na(i)) names(mf)[match("(id)", names(mf))] <- idname
    } else stop("id must be a simple variable name")

    # Get the list of states and istate
    if (missing(istate)) check <- survcheck2(Y, id)
    else {
        istate <- model.extract(mf, "istate")
        check <- survcheck2(Y, id, istate)
        }
    if (any(check$states == "censor")) states <- c("(censor)", check$states)
    else states <- c("censor", check$states)
    nstate <- length(check$states)
    # In the new data, there is 1 more row per subject.
    # newtime for subject i = c(first time1, time2) 
    # newstate is c(initial, Y[i,3])
    # for covariates, the last row of each subject is a repeat
    first <- !duplicated(id)
    last  <- !duplicated(id, fromLast=TRUE)
    n <- nrow(mf)

    indx1 <- rep(1:n, ifelse(first, 2, 1))
    indx2 <- rep(1:n, ifelse(last,  2, 1))
    newtime <- Y[indx1,2]
    newstat <- c(0L, match(attr(Y, "states"), check$states))[1L+ Y[indx1,3]]

    row1 <- duplicated(indx1, fromLast=TRUE) # first row of each subject
    newtime[row1] <- Y[first, 1]
    newstat[row1] <- as.numeric(check$istate[first])
    newdata <- cbind(data.frame("(time)"= newtime, 
                                "(state)"= factor(newstat, 0:nstate, states)),
                     mf[indx2,-1])
    row.names(newdata) <- NULL  # they are annoying and useless
    names(newdata)[1:2] <- c(tname, sname)
    indx <- match(c("(id)", "(istate)"), names(newdata), nomatch=0)
    newdata[, indx==0]  # remove the redundant (id) and (istate) columns
}
   
# The function from timeline to counting process is likely more useful
# A basic assumption is that covariate values are filled in when they are
#  observed, but may be NA when not.  We should skip NA values, just as tmerge
#  does.  The response is different.
fromtimeline <- function(formula, data, id, istate="istate") {
    Call <- match.call()
    keep <- match(c("formula", "data", "id"), names(Call), nomatch=0) 
    tcall <- Call[c(1, keep)]
    tcall[[1]] <- quote(stats::model.frame)
    mf <- eval(tcall, parent.frame())
    id <- model.extract(mf, "id")
    Y <- model.response(mf)

    if (!inherits(Y, "Surv")) stop("response must be a Surv object")
    type <- attr(Y, "type")
    if (!(type== "right" || type=="mright"))
        stop("only valid for a right censored response")

    # get a name for the resulting time and state variables, by parsing
    #  the formula.  Don't try too hard, though.
    tname <- c("tstart", "tstop")
    sname <- "state"  
    lhs <- formula[[2]]
    if ((is.name(lhs[[1]]) && lhs[[1]]== as.name("Surv")) ||
        deparse(lhs[[1]]) == "survival::Surv")  {
        # the lhs if of length 3
        if (is.name(lhs[[2]])) {
            temp <- as.character(lhs[[2]])
            tname <- paste0(temp, 1:2)
        }
        if (is.name(lhs[[3]]))  sname <- as.character(lhs[[3]])
        else if (is.call(lhs[[3]])){
            temp <- lhs[[3]]
            if (deparse(temp[[1]])== "factor" && is.name(temp[[2]]))
                sname <- as.character(temp[[2]])
        }
    } else if (is.name(lhs[[1]])) tname <- paste0(lhs[[1]], 1:2)
    
    # The id variable should be simple
    if (is.name(Call$id)) {
        idname <- as.character(Call$id)
        i <- match(idname, names(mf))
        if (is.na(i)) names(mf)[match("(id)", names(mf))] <- idname
    } else stop("id must be a simple variable name")
    
    # Use tmerge to do the work, since it understands missings
    tmin <- tapply(Y[,1], id, min)
    tmax <- tapply(Y[,1], id, max)
    if (any(tmin == tmax)) {
        warning("identifiers with only 1 row were removed")
        keep <- (tmin < tmax)
        mf <- mf[keep,]
        tmin <- tmin[keep]
        tmax <- tmax[keep]
        Y <- Y[,keep]
        id <- id[keep]
    }
    
    # mark variables like sex, treament, etc that are static
    static <- sapply(mf[,-1], 
                     function(x)(all(!is.na(x)) &&
                                 all(tapply(x, id, function(z) all(z==z[1])))))
    temp <- static
    temp["(id)" == names(temp)] <- FALSE  # don't include in new df
    d1 <- mf[!duplicated(id), c(FALSE, temp)] 

    # eval(parse(text)) is usually a sign of bad programming.  In this case
    #  I need to plug in the actual variable name of the id variable, however.
    #  The other choice of the non-syntactic name .. id = `(id)` doesn't
    #  seem to work.
    tcall <- c("tmerge(d1, cbind(d1, tmin=tmin, tmax=tmax), id=", idname,
               ",tstart=tmin, tstop=tmax,", 
               "options=list(tstartname=tname[1], tstopname= tname[2]))")
    new <- eval(parse(text=tcall))
    
    # Make sure the event name is not in the data set, as some other variable
    i <- 1
    while(!is.na(match(sname, names(mf)[-1]))){
        sname <- paste0(sname, i)
        i <- i+1
    }
    # add in the event, we assume no td variable has name of .y1. or .y2.
    temp <- !static
    temp[names(temp)== idname] <- TRUE  # force this variable into d2
    d2 <- data.frame(mf[, c(FALSE, temp), drop=FALSE],  
                     .y1. =Y[,1], .y2. = Y[,2])
    tcall <- paste("tmerge(new, d2, id=", idname, ", ", sname, 
                    "=event(.y1., .y2.)")
    
    if (any(!static)) { # there are tdc variables
        dname <- names(mf)[c(FALSE, !static)]
        temp <-  paste(paste(dname, "=tdc(.y1.," , dname,")"), collapse=', ')
        tcall <- paste(tcall, temp, sep=', ')
        }
    tcall <- paste(tcall, ")")

    new <- eval(parse(text=tcall))
    
    # add the istate variable, survcheck knows what to do
    if (any(istate == names(new))) 
        stop("istate option duplicates an existing variable")
    n <- nrow(new)
    itemp <- as.integer(new[[sname]])[c(1,seq(1, n-1))]
    itemp[!duplicated(new[[idname]])] <- Y[!duplicated(id),2]
    if (any(Y[!duplicated(id),2]==0)) 
        stop("no observation should start in a censored state")
    states <- attr(Y, "states")
    if (!is.null(states)) {
        nstate <- length(states)
        itemp <- factor(itemp, 1:nstate, states)
        new[[sname]] <- factor(new[[sname]], 0:nstate, c("censor", states))
    }
    new[[istate]] <- itemp
    new
}   
    
