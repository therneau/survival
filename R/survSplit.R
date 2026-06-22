survSplit <- function(formula, data, subset, na.action=na.pass, id,
                              cut, zero=0, episode, start="tstart",
                              end="tstop", event="event", added, timefix=TRUE) {
    Call <- match.call()
    if (missing(formula) || is.data.frame(formula)) {
        # an old style call
        # match arguments and build a formula
        warning("using survSplit without a formula is no longer supported")
        #stop("using survSplit without a formula is no longer supported")
        if (missing(data)) {
            if (!missing(formula)) {
                names(Call)[[2]] <- "data"
                # The line above is sneaky: it makes model.frame() work later
                data <- formula
            }
            else  stop("a data frame is required")
        }
        if (missing(end) || missing(event))
            stop("either a formula or the end and event arguments are required")

        if (!(is.character(event) && length(event) ==1 &&
              event %in% names(data)))
            stop("'event' must be a variable name in the data set")

        if (!(is.character(end) && length(end) ==1 &&
              end %in% names(data)))
            stop("'end' must be a variable name in the data set")
    
        if (!(is.character(start) && length(start)==1))
            stop("'start' must be a variable name")
        if (start %in% names(data)) temp <- paste(start, end, event, sep=',')
        else temp <- paste(end, event, sep=',')
        
        formula <- as.formula(paste("Surv(", temp, ")~ ."))
    }
    else if (missing(formula) || !inherits(formula, "formula")) 
        stop("either a formula or the end and event arguments are required")
 
    # create a call to model.frame() that contains the formula (required)
    #  and any other of the relevant optional arguments
    #  then evaluate it in the proper frame
    # one catch is that for Surv2 we need the id as an argument in the
    #  same way as id in coxph or survfit, but for Surv have historically
    #  allowed a character name, used to add an id for Surv(time, status) data. 
    # Checking for a character form of id is subtle, it works for id="charlie"
    #  but not for id=zed where zed is a length 1 character variable.  We have
    #  removed any mention of this from the help file, but
    if (!missing(id) && is.character(Call$id)){
        charid <- TRUE
        indx <- match(c("formula", "data", "subset", "na.action"), 
                      names(Call),nomatch=0)
    } else {
        charid <- FALSE
        indx <- match(c("formula", "data", "subset", "na.action", "id"), 
                      names(Call), nomatch=0)
    }

    if (indx[1]==0) temp$formula <- formula
    if (indx[2]==0) stop("a data argument is required")
    temp <- Call[c(1L,indx)]  # only keep the arguments we wanted
    if (missing(na.action)) temp$na.action <- quote(stats::na.pass)
    temp[[1L]] <- quote(stats::model.frame)  # change the function called
    mf <- eval.parent(temp)      

    Y <- model.response(mf)
    if (is.Surv(Y)) {
        if (!(attr(Y, "type") %in% c("right", "mright","counting","mcounting")))
            stop("not valid for ", attr(Y, "type"), " censored survival data")
    } else if (!inherits(Y, "Surv2"))
        stop ("the model must have a Surv or Surv2 object as the response")
    states <- attr(Y, "states")
    nY <- ncol(Y)
    ymiss <- is.na(Y)  # these pass through unchanged

    if (!is.numeric(cut) || any(!is.finite(cut)))
        stop("cut must be a vector of finite numbers")
    cut <- unique(sort(cut))
    ntimes <- length(cut)
    n <- nrow(data)

    # Deal with the near-ties problem
    if (!is.logical(timefix) || length(timefix) > 1)
        stop("invalid value for timefix option")
    if (timefix) {
        # to make this work, we have to include the cutpoints in the aeqSurv
        #  process
        if (inherits(Y, "Surv2")) {
            ty <- rbind(cbind(cut,0), Y)
            class(ty) <- "Surv2"  # not quite true, but okay for aeqSurv
            ty <- aeqSurv(ty)
            cut <- ty[1:ntimes,1]
            Y[,1] <- ty[-(1:ntimes),1]
        } else if (ncol(Y)==2) {
            ty <- rbind(cbind(cut,0), as.matrix(Y))
            class(ty) <- "Surv" # fool aeqSurv
            ty <- aeqSurv(ty)
            cut <- ty[1:ntimes, 1]
            Y[,1] <- ty[-(1:ntimes), 1]
        } else {
            ty <- rbind(cbind(min(Y[!ymiss,1]), cut, 0), as.matrix(Y))
            class(ty) <- "Surv" # fool aeqSurv
            ty <- aeqSurv(ty)
            cut <- ty[1:ntimes,2]
            Y[,1:2] <- ty[-(1:ntimes), 1:2]
        }
    }

    # If the RHS was ~. and no rows were removed by an na.action,
    #  then modify data, not the model frame. It makes "id"
    #  pass through properly, and all variables end stay in the same position
    rightdot <- (is.name(formula[[3]]) && formula[[3]]== as.name('.') &&
                 nrow(mf) == nrow(data))
    # If charid and simple survival and the id variable name is not already in
    # the data, then add it. For anything other than simple survival we don't
    # have the information to add one.
    if (charid && is.na(match(id, names(data))) && is.Surv(Y) &&
        attr(Y, "type") == "right") {
        data[[id]] <- 1:nrow(data)
        mf[[id]] <- 1:nrow(data)
    }

    # Now for the actual work
    if (inherits(Y, "Surv2")) {
        id <- model.extract(mf, "id")
        if (is.null(id)) stop("an id statement is required")
        idord <- order(id, Y[,1])
        # Create a fake 'stop' for each Y[,1] element, which is the next
        #  time for this id, or a dup of Y[j,1] if j is the last for the id
        storage.mode(Y) <- "double"
        yfake <- Y[,1]
        notlast <- duplicated(id[idord], fromLast=TRUE)
        yfake[idord[notlast]] <- yfake[c(idord[-1], 0)[notlast]]
        # We make it look like (time1, time2) data for the .C call,
        # last interval has length 0
        index <- .Call(Csurvsplit, Y[,1], yfake, as.double(cut))
 #       browser()
        status <- Y[index$row,2]
        # The C routine marks the parent of each insertion as censored, which
        #  is correct for (time1, time2) data, but for Surv2 the insertion is
        #  censored row
        status[1L + which(index$censor)] <- 0  
    }

    else {
        if (nY ==2) {
            if (any(Y[!ymiss,1] <= zero))
                stop("'zero' parameter must be less than any observed times")
            Y <- cbind(zero, Y)
        }
        temp <- (Y[!ymiss,1] >= Y[!ymiss,2])
        if (any(temp)) stop("start time must be < stop time")
        
        storage.mode(Y) <- "double"
        index <- .Call(Csurvsplit, Y[,1], Y[,2], as.double(cut))
        
        status <- Y[index$row,3]
        status[index$censor] <- 0
    }

    if (rightdot) newdata <- data[index$row,, drop=FALSE]
    else {
        newdata <- mf[index$row, -1, drop=FALSE]
        attr(newdata, "terms") <- NULL  # look less like a model frame
    }

    # the factor line needs 0:length(states) in case there is a state
    #  with no representatives, then status will be missing that category
    if (!is.null(states))  
        status <- factor(status, 0:length(states),labels=c("censor", states))

    # The remaining work is for sensible names for time/status in the output
    # 1. Easy case: the user has a Surv or Surv2 with sensible variable names
    if (inherits(formula[[2]], "call") && (formula[[2]][[1]]== as.name("Surv")
        || formula[[2]][[1]]== as.name("Surv2"))) {
        # it was a call, figure out the names
        # The user might have used something like Surv(status=abc, time=fred),
        #  so use match.call to find "abc" and "fred".  But give up if there
        #  is anything complex. 
        if (inherits(Y, "Surv2")) {
            temp <- match.call(Surv2, formula[[2]])
            if (!is.null(temp[["time"]]) && is.name(temp[["time"]]))
                start <- as.character(temp$time)
            if (!is.null(temp[["event"]]) && is.name(temp[["event"]]))
                event <- as.character(temp$event)
            newdata[[start]] <- index$start
            newdata[[event]] <- status
        } else { 
            temp <- match.call(Surv, formula[[2]])
            if (nY==2) { 
                if (missing(end) && !is.null(temp[["time"]]) 
                    && is.name(temp[["time"]]))
                    end <- as.character(temp$time) # time might match 'time2'
                if (missing(event) && !is.null(temp$time2) &&
                    is.name(temp$time2)) event <- as.character(temp$time2)
                if (missing(event) && !is.null(temp$event) && 
                    is.name(temp$event)) event <- as.character(temp$event)
            }
            else {
                if (missing(end) && !is.null(temp[["time"]]) 
                    && is.name(temp["time"]))
                    start <- as.character(temp[["time"]])
                if (missing(end) && !is.null(temp$time2) && is.name(temp$time2)) 
                    end <- as.character(temp$time2)
                if (missing(event) && !is.null(temp$event) && 
                    is.name(temp$event)) event <- as.character(temp$event)
                if (missing(start) && !is.null(temp$time) && is.name(temp$time))
                    start <- as.character(temp$time)
            }       
            newdata[[start]] <- index$start
            newdata[[end]]   <- index$end
            newdata[[event]] <- status
        }
    } else {
        # They have a Surv or Surv2 object
        if (!inherits(formula[[2]], "name"))
            stop("trouble making new names, left hand side not recognized")
        temp <- as.character(formula[[2]])
        if (is.Surv(Y))
            newdata[[temp]] <- Surv(index$start, index$end, status)
        else newdata[[temp]]<- Surv2(Ynew[,1], status)
    }

    if (!missing(episode)) {
        if (!is.character(episode)) stop("episode must be a character string")
        newdata[[make.names(episode)]] <- index$interval +1
    }
    if (!missing(added)) {
        if (!is.character(added)) stop("added must be a character string")
        newdata[[make.names(added)]] <- index$censor
    }
        
    row.names(newdata) <- NULL    # erase R's manufactured row names
    newdata
}
