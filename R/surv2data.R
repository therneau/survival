#
# Routine to turn a "timeline" model frame into a "counting process" model
#  frame. 
# The mf arg is almost certainly the result of a call to coxph or survfit (or a
#  survfit.coxphms or residuals.coxphms call that had to recreate the data), and
#  that fcn is using this function to transform to CP style. As such, there will 
#  normally be a column named "(id)", and perhaps "(istate)" and/or "(cluster)".
# The task is fairly simple:
#   1. An id with 10 rows will have 9 rows in the new data set, 2-10 
#     contribute the (time2, outcome) part of the survival response, and
#     original rows 1-9 everything else.
#   2. Missing covariates are filled in using last-value-carried forward.
#   3. If there is no '(istate)', and normally there won't be because a
#     coxph call doesn't need that argument for Surv2 data, then create it.
#
# The repeat attribute of the Surv2 object determines whether events can
#   "stutter", i.e., two of the same type that are adjacent, or only have
#    censored times between them, count as a second event.
#
surv2counting <- function(mf, repeated=FALSE, lvcf=TRUE) {
    Terms <- terms(mf)
    y <- model.response(mf)
    n <- nrow(y)
    states <- attr(y, "states")
    # the next line supports Surv2 objects
    if (missing(repeated) && !is.null(attr(y, "repeated"))) 
        repeated <- attr(y, "repeated")

    id <- model.extract(mf, "id")
    if (length(id) != n) stop("id statement is required")

    # relax this some later day (or not?)
    if (any(is.na(id)) || any(is.na(y[,1])))
        stop("id and time cannot be missing")
    
    # The data isn't necessarity sorted, don't reorder it
    isort <- order(id, y[,1])
    id2 <- id[isort]
    y2  <- y[isort,]
    first <- !duplicated(id2)
    last  <- !duplicated(id2, fromLast=TRUE)
    # check for duplcate times
    temp <- data.frame(id2, y2[,1]) 
    if (any(duplicated(temp)))
        stop("duplicated time for an id")
        
    mf2 <- mf[isort,][!last,]  # starting point for the output

    # The LVCF operation replaces any missing values with the most recent non
    #  missing.  This is easiest to do with the tmerge3 C routine, which expects
    #  the data to be in time within id order
    # Use LVCF on all colums except: doesn't have any missings (not needed),
    #  column 1= response, those with () names, normally (id), (cluster)
    id3 <- id2[!last]  # I need idi again later
    if (is.factor(id3)) idi <- as.integer(id3)
    else    idi <- match(id3, unique(id3))  # tmerge3 wants an integer id

    if (lvcf) {
        skip <- c(1, which(substring(names(mf), 1, 1) == "("))
        for (i in (1:ncol(mf2))[-skip]) {
            miss <- is.na(mf2[[i]])
            if (any(miss)) {
                k <- .Call(Ctmerge3, idi, miss)
                # k is the row number of the value to be carried forward
                update <- (miss & k>0)  # some values will be updated
                if (any(update))  { # some successful replacements
                    if (is.matrix(mf2[[i]])) 
                        mf2[[i]][update,] <-mf2[[i]][k[update],]
                    else mf2[[i]][update] <- (mf2[[i]])[k[update]]
                }
            }
        }
    }
    
    # If there are missing status values, assume that NA is actually
    #  the code for censored.
    if (any(is.na(y2[,2]))) {
        if (is.null(states)) {
            # Really? The user had NA, 1, 2 and didn't use a factor?
            # I'm going to guess not
        } else { # the user input was a factor
            newstate <- attr(y, "inputAttributes")$event$levels
            if (is.null(newstate)) {cat("surv2data bug "); browser()}
            y2[,2] <- ifelse(is.na(y2[,2]), 0, 1+ y2[,2])
            attr(y2, "states") <- newstate
            states <- newstate
        }
    }
    # Create a current state vector, which is just LVCF using the state.
    # NA or censored  treated as "not present" for the tmerge3 call.  
    # Y[,2]=0 is censored for either single state or multistate timeline
    y3 <- y2[!last]
    censored <- ifelse(is.na(y2[,2]), TRUE, y2[,2] ==0) # sorted version
    if (all(censored[first])) {
        # special case -- no initial state for anyone, don't create istate
        # this occurs for ordinary 0/1 status, or perhaps competing risks
        istate2 <- NULL
    } 
    else if (any(censored[first]))
        stop("everyone or no one should have an initial state")
    else {
        istate2 <- y2[!last,2]
        censored <- censored[!last]
        if (any(censored)) { # 0 censored rows is possible, but rare
            k <- .Call(Ctmerge3, idi, censored) # replace censored rows
            istate2[censored] <- istate2[k[censored]]
        }
        mf2[, "(istate)"] <- factor(istate2, 1:length(states), states)

        istate <- model.extract(mf, "istate")
        if (!is.null(istate)) {
            istate <- model.extract
            # The user had an istate= argument in their multistate (coxph) or 
            #  Aalen-Johansen (survfit) call.  It should agree, perfectly.  
            # Well, they could legally put the states in a different order
            xx <-"istate argument does not agree with initial Surv2 values"
            if (is.numeric(istate) && any(istate2 != istate)) stop(xx)
            else {
                istate <- as.factor(istate)  # in case is is type character
                itemp <- match(levels(istate), states, nomatch=0)
                if (any(itemp==0)) stop(xx)
                if (any(istate2[itemp] != as.numeric(istate))) stop(xx)
            }
        }            
    }
    # Now replace the Surv2 response with one of the standard forms
    #    "right"   :  (time, status) 
    #    "counting":  (time1, time2, status)
    #    "mright"  :  (time, state)
    #    "mcounting": (time1, time2, state)
    status <- y2[!first, 2] # integer version
    if (!repeated  && !is.null(istate2)) 
        status[status== istate2] <- 0 # no repeat of current state
    if (!is.null(states)) 
        status <- factor(status, 0:length(states), c("censor", states))

    tstart <- y2[!last, 1]
    tstop  <- y2[!first,1]
    # right or mright
    if (!any(duplicated(id3)))  # one obs per id
        mf2[[1]] <- Surv(tstop, status)
    else
        mf2[[1]] <- Surv(tstart, tstop, status)

    #put the data back into the original order
    jj <- order(isort[!last])  # this line is not obvious, but it works!
    mf2[jj,]
}
   
# User callable version
Surv2data <- function(formula, data, subset, id){
    .Deprecated("timeline2counting", old="Surv2data")
    Call <- match.call()
    indx <- match(c("formula", "data", "weights", "subset", "na.action",
                    "cluster", "id"),
                  names(Call), nomatch=0)
    if (indx[1] ==0) stop("A formula argument is required")
    tform <- Call[c(1, indx)]  # only keep arguments we wanted
    tform$na.action <- stats::na.pass
    tform[[1L]] <- quote(stats::model.frame)
    mf <- eval(tform, parent.frame())
    if (!inherits(model.response(mf), "Surv2"))
        stop("the response must be a Surv2 object")
    surv2counting(mf)
}

timeline2counting <- function(formula, data, subset, id, repeated= FALSE,
                              lvcf = TRUE, yname=c("tstart", "tstop", "status")){
    # Get the formula for the response
    Call <- match.call()
    indx <- match(c("formula", "data", "subset", "id"),
                  names(Call), nomatch=0)
    if (indx[1] ==0) stop("a formula argument is required")
    if (indx[2] ==0) stop("the data argument is required")
    if (indx[4] ==0) stop("the id argument is required")
    tform <- Call[c(1, indx)]  # only keep arguments we wanted
    tform$na.action <- stats::na.pass
    tform[[1L]] <- quote(stats::model.frame)
    mf <- eval(tform, parent.frame())
    
    Y <- model.response(mf)
    id <- model.extract(mf, "id")
    if (!(is.Surv(Y) || inherits(Y, "Surv2")))
        stop("response must be a survival object")
    if (is.Surv(Y) && attr(Y, "type") %in% c("counting", "mcounting"))
        stop("response cannot be of counting process type")
    if (is.null(id) || !any(duplicated(id)))
        stop("data does not appear to be timeline data")
    new <- surv2counting(mf, repeated=repeated, lvcf=lvcf)

    # Remove (id) and rename (istate) to istate
    ii <- match("(istate)", names(new)) # added by surv2counting
    if (!is.na(ii)) { # not added in some cases
        if (!any(names(mf)== "istate")) names(new)[ii] <- "istate"
    }
    new["(id)"] <- NULL

    #
    # Give new names to response, so that the data set is more useable
    #
    yy <- new[[1]] # the response is always variable 1
    states <- attr(yy, "states")
    ny <- ncol(yy)  # ny=2 for competing risks, for instance
    if (is.null(states)) status <- yy[,ny]     
    else status <- factor(yy[,ny], 0:length(states), c("censor", states))
    if (ncol(yy) ==3) 
        tdata <- data.frame(tstart= yy[,1], tstop = yy[,2], status)
    else tdata <- data.frame(tstart=yy[,1], status)

    if (!missing(yname)) {
        if (any(!is.na(match(yname, names(new)))))
            stop("element of yname conflicts with an existing name in the data")
        if (ny==2) {
            if (!(length(yname) %in% 2:3)) stop("wrong length for yname")
            names(tdata) <- c(yname[1], yname[length(yname)])
        } else {
            if (length(yname) != 3) stop("wrong length for yname") 
            names(tdata) <- yname
        }
    } else {
        # Use a best guess, if possible
        resp <- formula[[2]] # "Surv(age, dead)", say
        if (is.name(resp[[2]])) {
            if (ny==2) yname <- c(as.character(resp[2]), yname[3])
            else yname[1:2] <- paste0(as.character(resp[2]), 1:2)
        }
        if (is.name(resp[[3]])) yname[length(yname)] <- as.character(resp[[3]])
        
        conflict <- match(yname, names(new), nomatch=0)
        if (any(conflict >0))
            yname[conflict>0] <- paste0("_", yname[conflict>0], "_")
        names(tdata) <- yname
    }
    cbind(new[,-1,drop=FALSE], tdata)
}



