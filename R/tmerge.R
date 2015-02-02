# Automatically generated from all.nw using noweb
tmerge <- function(data1, data2, id, ..., tstart, tstop, options) {
    Call <- match.call()
    # The function wants to recognize special keywords in the
    #  arguments, so define a set of functions which will be used to
   #  mark objects
    new <- new.env(parent=parent.frame())
    assign("tdc", function(...) {x <- list(...); class(x) <- "tdc"; x},
           envir=new)
    assign("cumtdc", function(...) {x <- list(...); class(x) <-"cumtdc"; x},
           envir=new)
    assign("event", function(...) {x <- list(...); class(x) <-"event"; x},
           envir=new)
    assign("cumevent", function(...) {x <- list(...); class(x) <-"cumevent"; x},
           envir=new)

    if (missing(data1)) stop("the data1 argument is required")
    if (missing(id)) stop("the id argument is required")
    
    defopt <- function(tname, defer =0, ...) {
        temp <- c(id="id", tstart="tstart", tstop="tstop")
        if (missing(tname)) tname <- temp
        else {
            if (is.list(tname)) tname <- unlist(tname) # a common mistake, unlist it
            if (!is.character(tname)) 
                stop("tname option must contain variable names")
            indx <- match(names(tname), temp)
            if (any(is.na(indx))) stop("tname option contains an unrecogized name")
            temp[names(tname)] <- tname 
            tname <- temp
            }
        if (length(defer) !=1 || !is.numeric(defer) || defer <0)
            stop("defer option must be a non-negative number")
        list(tname=tname, defer=defer)
    }
    if (!missing(options)) {
        if (!is.list(options)) stop("options must be a list")
        topt <- do.call(defopt, options)
    }
    else topt <- defopt()

    # id, tstart, tstop are found in data2, if it is present
    if (!missing(data2)) {
        id <- eval(Call[["id"]], data2)
        if (!missing(tstart)) tstart <- eval(Call[["tstart"]], data2)
        if (!missing(tstop))  tstop <-  eval(Call[["tstop"]],  data2)
    }

    if (!missing(tstart) && length(tstart) != length(id))
        stop("tstart and id must be the same length")
    if (!missing(tstop) && length(tstop) != length(id))
        stop("tstop and id must be the same length")
    # grab the... arguments
    notdot <- c("data1", "data2", "id", "tstart", "tstop", "tname")
    dotarg <- Call[is.na(match(names(Call), notdot))]
    dotarg[[1]] <- as.name("list")  # The as-yet dotarg arguments
    if (missing(data2)) args <- eval(dotarg, envir=new)
    else  args <- eval(dotarg, data2, enclos=new)
        
    argclass <- sapply(args, function(x) (class(x))[1])
    argname <- names(args)
    if (any(argname== "")) stop("all additional argments must have a name")
           
    check <- match(argclass, c("tdc", "cumtdc", "event", "cumevent"))
    if (any(is.na(check)))
        stop(paste("argument(s)", argname[is.na(check)], 
                       "not a recognized type"))

    # The tcount matrix is useful for debugging
    if (length(argname)) {
        tcount <- matrix(0L, length(argname), 8)
        dimnames(tcount) <- list(argname, c("early","late", "gap", "within", 
                                      "tied edge", "front edge", "back edge",
                                      "tied times"))
        }
    tevent <- attr(data, "tevent") # event type variables
    newdata <- data1  # make a copy
    n <- length(id)   # number of things to add
    if (!missing(tstop) || !inherits(data1, "tmerge")){ 
        # This is a first call
        tname <- c(id="id", tstart="tstart", tstop="tstop") 
        if (n != nrow(data1)) 
            stop("expect nrow(data1)=length(id), on a first call")
        newdata$id <- id
        if (!missing(tstop)) newdata$tstop <- tstop
        else stop("the tstop argument is required on the first call")

        if (!missing(tstart))  newdata$tstart <- tstart
        else newdata$tstart <- rep(0, nrow(newdata))
        if (any(tstart >= tstop)) 
                stop("tstop must be > tstart for all observations")

        # If there are duplicated ids, then we need to ensure that each subject
        #  is a sequential set of observations, sorted by time
        # It tstart was not supplied, we need to correct our "0" created above
        if (any(duplicated(id))) {
            indx <- order(newdata$id, newdata$tstop)
            newdata <- newdata[indx,]
            id <- id[indx]
            if (length(tstart) ==0) 
                newdata$tstart <- c(0, ifelse(duplicated(id), newdata$tstop[-n],0))
            else {
                ok <- (duplicated(id) | newdata$tstart[-1] >= newdata$tstop[-n])
                if (any(!ok)) stop("overlapping time intervals for a subject")
            }
        }
    }
    else { #not a first call
        tname <- attr(data1, "tname")
        if (length(tname) != 3 || any(names(tname) != c("id", "tstart", "tstop")))
            stop("corrupted tname attribute in data1")
        indx <- match(tname, names(data1))
        if (any(is.na(indx))) 
            stop("data1 does not match its own tname attribute")
        if (any(is.na(match(id, data1[[tname[1]]]))))
            stop("id values found in data2 which are not in data1")
    }
    for (i in 1:length(args)) {
        baseid <- newdata[[tname["id"]]]
        dstart <- newdata[[tname["tstart"]]]
        dstop  <- newdata[[tname["tstop"]]]
        
        # if an event time is missing then skip that obs
        etime <- args[[i]][[1]]
        keep <- !is.na(etime)
        etime <- etime[keep]
        id <- id[keep]
        if (length(etime) != n) 
            stop("argument", argname[i], "is not the same length as id")
        
        # For an event or cumevent, one of the later steps becomes much
        #  easier if we sort the new data by id and time
        indx <- order(id, etime)
        id <- id[indx]
        etime <- etime[indx]
        if (length(args[[i]]) > 1)
            yinc <- (args[[i]][2])[indx]
            
        # indx1 points to the closest start time in the baseline data (data1)
        #  that is <= etime.  indx2 to the closest end time that is >=etime.
        # If etime falls into a (tstart, tstop) interval, indx1 and indx2
        #   will match
        # If the "defer" argument is set and this event is of type tdc, then
        #   any event times are artificially moved left by "defer" amount wrt
        #   doing the indx2 match.  This will cause an insertion that is too close
        #   to an event to be labeled as itype=3 (or itype=2 if this was the last
        #   interval for the subject) and so map later in time.
        deftime <- rep(0., n)
        if (dopt$defer >0 && length(event) && 
            argtype[i] %in% c("tdc", "cumtdc")) {
            for (k in event) {
                temp <- edata[[event[k]]]
                if (is.logical(temp)) deftime[temp] <- dopt$defer
                else deftime[temp!=0] <- dopt$defer
            }
        }
        indx1 <- neardate(id, baseid, etime, dstart, best="prior")
        indx2 <- neardate(id, baseid, etime, dstop+deftime, best="after")

        # The event times fall into one of 5 categories
        #   1. Before the first interval
        #   2. After the last interval
        #   3. Outside any interval but with time span, i.e, it falls into
        #       a gap in follow-up
        #   4. Strictly inside an interval (does't touch either end)
        #   5. Inside an interval, but touching.
        itype <- ifelse(is.na(indx1), 1,
                        ifelse(is.na(indx2), 2, 
                               ifelse(indx2 > indx1, 3,
                                      ifelse(etime== dstart[indx1] | 
                                             etime== dstop[indx2], 5, 4))))

        # Subdivide the events that touch on a boundary
        #  1: intervals of (a,b] (b,d], new count at b  "tied edge"
        #  2: intervals of (a,b] (c,d] with c>b, new count at c, "front edge"
        #  3: intervals of (a,b] (c,d] with c>b, new count at b, "back edge"
        #
        subtype <- ifelse(itype!=5, 0, 
                          ifelse(indx1 == indx2+1, 1,
                                 ifelse(etime==dstart[indx1], 2, 3)))
        tcount[i,1:7] <- table(factor(itype+subtype, levels=c(1:4, 6:8)))

        # count ties.  id and etime are not necessarily sorted
        tcount[i,8] <- sum(sapply(etime, id, function(x) sum(duplicated(x))))
        # Look to see if this term has one or two arguments.  If one arg
        #  then the increment is 1, else it is the second arg.  The myfun()
        #  function will later compute totals by unique subject/time pair
        #
        if (length(args[[i]]) >1) {
            if (length(args[[i]]) > 2) 
                stop("too many variables in an", argclass[i], "call")
            if (diff(sapply(args[[i]], length)) !=0)
                stop("different lengths in an", argclass[i], "call")
            if (!is.numeric(yinc) && argclass[i] != "event") 
                stop("non numeric increment in an", argclass[i], "call")
            myfun <- function(x, grp) {
                temp <- tapply(yinc[grp], x[grp], sum)
                ifelse(is.na(temp), 0, temp)
            }
        }
        else {
            myfun <- function(x, grp) table(x[grp])
            yinc <- rep(1L, length(etime))     # each counts as 1
        }
        indx4 <- which(itype==4)
        n4 <- length(indx4)
        icount <- tapply(etime[indx4], indx1[indx4], function(x) sort(unique(x)))
        n.add <- sapply(icount, length)  #number of rows to add

        # expand the data 
        irep <- rep.int(1L, nrow(newdata))
        irep[as.integer(names(icount))] <- 1+ n.add # number of rows in new data
        newdata <- newdata[irep,]  #expand it out
        dstart <- newdata[[tname[2]]]
        dstop <-  newdata[[tname[3]]]

        #fix up times
        temp <- vector("list", n4)
        iend <- (cumsum(irep))[irep >1]  #end row of each duplication set
        for (i in 1:n4) temp[[i]] <-  -(seq(n.add[i] -1, 0)) + iend[i]
        newrows <- unlist(temp)
        dstart[newrows] <- dstop[newrows-1] <- unlist(icount)
        newdata[[tname["tstart"]]] <- dstart
        newdata[[tname["tstop"]]]  <- dstop
        for (i in events) newdata[newrows-1, event[i]] <- 0
        if (dopt$defer > 0) {
            deftime <- deftime[irep]
            deftime[newrows] <- 0
        }

        # refresh indices
        indx1 <- neardate(id, baseid, etime, dstart, best="prior")
        indx2 <- neardate(id, baseid, etime, dstop+deftime , best="after")
        subtype[itype==4] <- 1  #all the "insides" are now on a tied edge
        itype[itype==4]   <- 5  
        n <- nrow(newdata)
        if (argclass[i] %in% c("cumtdc", "cumevent")) 
            yinc <- unlist(tapply(yinc, id, cumsum))

        newvar <- rep(0, n) 
        if (argclass[i] %in% c("event", "cumevent")) {
            newvar <- rep(0, n) 
            keep <- (subtype==1 | subtype==3) # all other events are thrown away
            newvar[indx1[keep]] <- yinc[keep]
        }
        else {
            keep <- itype != 2  # changes after the last interval are ignored
            indx <- ifelse(subtype==1, indx1, 
                           ifelse(subtype==2, indx2, indx2 +1))

            # id can be any data type; feed integers to the C routine
            newvar <- .Call("tmerge", match(baseid, baseid), dstart, dstop, 
                            match(id, baseid)[keep], etime[keep], 
                            yinc[keep], indx[keep])
        }
        newdata[[argname[i]]] <- newvar
        if (argclass %in% c("event", "cumevent")) tevent <- c(tevent, argnames[[i]])
    }
    attr(newdata, "tmerge") <- tmerge
    attr(newdata, "tcount") <- tcount
    if (length(tevent)) attr(newdata, "tevent") <- tevent
    row.names(newdata) <- NULL
    newdata
}
