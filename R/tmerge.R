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
    if (!inherits(data1, "data.frame")) stop("data1 must be a data frame")
    
    tmerge.control <- function(id="id", tstart="tstart", tstop="tstop",  defer =0) {
        if (length(defer) !=1 || !is.numeric(defer) || defer <0)
            stop("defer option must be a non-negative number")
        if (!is.character(id)) stop("id option must be a character string")
        if (!is.character(tstart)) stop("tstart option must be a character string") 
         if (!is.character(tstop)) stop("tstop option must be a character string") 
       list(id=id, tstart=tstart, tstop=tstop, defer=defer)
    }

    tname <- attr(data1, "tname")
    if (!is.null(tname) && any(is.null(match(unlist(tname), names(data1)))))
        stop("data1 does not match its own tname attribute")
    if (!missing(options)) {
        if (!is.list(options)) stop("options must be a list")
        if (!is.null(tname)) {
            # Changing a name partway through a set of calls?
            if (any(!is.na(match(names(options), names(tname)))))
                stop("cannot change names in mid-stream")
            topt <- do.call(tmerge.control, c(tname, options))
            }
        else topt <- do.call(tmerge.control, options)
    }
    else if (length(tname)) topt <- do.call(tmerge.control, tname)
    else  topt <- tmerge.control()

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
    notdot <- c("data1", "data2", "id", "tstart", "tstop", "topt")
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
    tcount <- matrix(0L, length(argname), 8)
    dimnames(tcount) <- list(argname, c("early","late", "gap", "within", 
                                      "boundary", "leading", "trailing",
                                      "tied times"))
    tevent <- attr(data1, "tevent") # event type variables
    newdata <- data1  # make a copy
    if (!missing(tstop) || length(tname)==0) { 
        # This is a first call
        indx <- match(c(topt$id, topt$tstart, topt$tstop), names(data1), nomatch=0)
        if (all(indx[1:2] >0) && missing(tstop)) {
            # case 1 above, just some data checks
            if (indx[3]==0) newdata[[topt$tstart]] <- 0
            if (!is.numeric(newdata[[topt$tstart]]) || 
                !is.numeric(newdata[[topt$tstop]]))
                stop("start and end variables must be numeric")
            if (any(newdata[[topt$tstart]]  >= newdata[[topt$tstop]])) 
                stop("stop time must be > start time for all observations")

            # If there are duplicated ids, then we need to ensure that each subject
            #  is a sequential set of observations, sorted by time (just sort it)
            # If tstart was not supplied, we need to correct our "0" created above
            baseid <- data1[[topt$id]]
            if (any(duplicated(baseid))) {
                indx <- order(newdata$id, newdata$tstop)
                if (any(indx != seq(along=index))) {
                    # sort the data
                    newdata <- newdata[indx,]
                    baseid <- baseid[indx]
                }
                newid <- !duplicated(baseid)   #is this row a new id value?
                n <- nrow(newdata)
                if (length(tstart) ==0) 
                    newdata$tstart <- ifelse(newid, 0, c(0, newdata$tstop[-n]))
                else {
                    ok <- (newid[-1] | newdata$tstart[-1] >= newdata$tstop[-n])
                    if (any(!ok)) stop("overlapping time intervals for a subject")
                }
            }
        } else {
            if (missing(tstop)) {
                # case 3 above
                if (length(argclass)==0 || argclass[1] != "event")
                    stop("neither a tstop argument nor an initial event argument was found")
                tstop <- args[[1]][[1]]
                }
            # case 2 and case 3
            if (any(is.na(tstop))) 
                stop("missing time value, when that variable defines the span")
            if (missing(tstart)) tstart <- rep(0, length(id))
            if (any(tstart >= tstop)) 
                stop("stop time must be > start time for all observations")

            if (indx[1] >0) { # the id variable is in data1
                baseid <- data1[[topt$id]]
                if (any(duplicated(baseid))) 
                    stop("duplicate identifiers in data1")
                indx2 <- match(id, baseid)
                if (any(is.na(indx2)))
                    stop("'id' has values not in data1")
                }
            else {
                if (nrow(data1) != nrow(data2))
                    stop("nrow(data1) != nrow(data2) and data1 is missing the id")
                indx2 <- seq.int(along=id)
                newdata[topt$id] <- id
                }
            newdata[indx2, topt$tstart] <- tstart
            newdata[indx2, topt$tstop]  <- tstop
        }
    }
    else { #not a first call
        if (any(is.na(match(id, data1[[topt$id]]))))
            stop("id values found in data2 which are not in data1")
    }
    saveid <- id
    for (ii in seq(along=args)) {
        argi <- args[[ii]]
        baseid <- newdata[[topt$id]]
        dstart <- newdata[[topt$tstart]]
        dstop  <- newdata[[topt$tstop]]
        
        # if an event time is missing then skip that obs
        etime <- argi[[1]]
        keep <- !is.na(etime)
        etime <- etime[keep]
        id <- saveid[keep]
        if (length(etime) != length(id)) 
            stop("argument", argname[ii], "is not the same length as id")
        
        # For an event or cumevent, one of the later steps becomes much
        #  easier if we sort the new data by id and time
        indx <- order(id, etime)
        id <- id[indx]
        etime <- etime[indx]
        if (length(argi) > 1)
            yinc <- (argi[[2]])[indx]
            
        # indx1 points to the closest start time in the baseline data (data1)
        #  that is <= etime.  indx2 to the closest end time that is >=etime.
        # If etime falls into a (tstart, tstop) interval, indx1 and indx2
        #   will match
        # If the "defer" argument is set and this event is of type tdc, then
        #   any event times are artificially moved left by "defer" amount wrt
        #   doing the indx2 match.  This will cause an insertion that is too close
        #   to an event to be labeled as itype=3 (or itype=2 if this was the last
        #   interval for the subject) and so map later in time.
        defer <- rep(0., nrow(newdata))
        if (topt$defer >0 && length(event) && 
            argtype[ii] %in% c("tdc", "cumtdc")) {
            for (ename in tevent) {
                temp <- newdata[[ename]]
                if (is.logical(temp)) defer[temp] <- topt$defer
                else defer[temp!=0] <- topt$defer
            }
        }
        indx1 <- neardate(id, baseid, etime, dstart, best="prior")
        indx2 <- neardate(id, baseid, etime, dstop+defer, best="after")

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
        tcount[ii,1:7] <- table(factor(itype+subtype, levels=c(1:4, 6:8)))

        # count ties.  id and etime are not necessarily sorted
        tcount[ii,8] <- sum(tapply(etime, id, function(x) sum(duplicated(x))))
        # Look to see if this term has one or two arguments.  If one arg
        #  then the increment is 1, else it is the second arg.  The myfun()
        #  function will later compute totals by unique subject/time pair
        #
        if (length(argi) >1) {
            if (length(argi) > 2) 
                stop("too many variables in an", argclass[ii], "call")
            if (diff(sapply(argi, length)) !=0)
                stop("different lengths in an", argclass[ii], "call")
            if (!is.numeric(yinc) && argclass[ii] != "event") 
                stop("non numeric increment in an", argclass[ii], "call")
            myfun <- function(x, grp) {
                temp <- tapply(yinc[grp], x[grp], sum)
                ifelse(is.na(temp), 0, temp)
            }
        }
        else {
            myfun <- function(x, grp) table(x[grp])
            yinc <- rep(1.0, length(etime))     # each counts as 1
        }
        indx4 <- which(itype==4)
        n4 <- length(indx4)
        if (n4 > 0) {
            icount <- tapply(etime[indx4], indx1[indx4], function(x) sort(unique(x)))
            n.add <- sapply(icount, length)  #number of rows to add 
            
            # expand the data 
            irep <- rep.int(1L, nrow(newdata))
            erow <- unique(indx1[indx4])  # which rows of the data to expand
            irep[erow] <- 1+ n.add # number of rows in new data
            jrep <- rep(1:nrow(newdata), irep)  #stutter the duplicated rows
            newdata <- newdata[jrep,]  #expand it out
            dstart  <- dstart[jrep]
            dstop   <- dstop[jrep]

            #fix up times
            nfix <- length(erow)
            temp <- vector("list", nfix)
            iend <- (cumsum(irep))[irep >1]  #end row of each duplication set
            for (j in 1:nfix) temp[[j]] <-  -(seq(n.add[j] -1, 0)) + iend[j]
            newrows <- unlist(temp)
            dstart[newrows] <- dstop[newrows-1] <- unlist(icount)
            newdata[[topt$tstart]] <- dstart
            newdata[[topt$tstop]]  <- dstop
            for (ename in tevent) newdata[newrows-1, ename] <- 0
            if (topt$defer > 0) {
                defer <- defer[jrep]
                defer[newrows] <- 0
            } else defer <- rep(0, nrow(newdata))

            # refresh indices
            baseid <- newdata[[topt$id]]
            indx1 <- neardate(id, baseid, etime, dstart, best="prior")
            indx2 <- neardate(id, baseid, etime, dstop+ defer , best="after")
            subtype[itype==4] <- 1  #all the "insides" are now on a tied edge
            itype[itype==4]   <- 5  
        }
        # add it in
        if (argclass[ii] %in% c("cumtdc", "cumevent")) 
            yinc <- unlist(tapply(yinc, id, cumsum))

        newvar <- newdata[[argname[ii]]]  #does the variable exist? 
        if (argclass[ii] %in% c("event", "cumevent")) {
            if (is.null(newvar)) newvar <- rep(0, nrow(newdata)) 
            keep <- (subtype==1 | subtype==3) # all other events are thrown away
            newvar[indx2[keep]] <- yinc[keep]
            tevent <- unique(c(tevent, argname[ii]))
        }
        else {
            keep <- itype != 2  # changes after the last interval are ignored
            indx <- ifelse(subtype==1, indx1, 
                           ifelse(subtype==3, indx2+1, indx2))

            if (is.null(newvar)) {
                if (length(argi)==1) newvar <- rep(0.0, nrow(newdata))
                else newvar <- rep(NA_real_, nrow(newdata))
            }
            # id can be any data type; feed integers to the C routine
            storage.mode(dstop) <- storage.mode(newvar) <- "double"
            storage.mode(etime) <- storage.mode(yinc) <- "double"
            newvar <- .Call("tmerge", match(baseid, baseid), dstart, dstop, newvar, 
                            match(id, baseid)[keep], etime[keep], 
                            yinc[keep], indx[keep])
        }

        newdata[[argname[ii]]] <- newvar
    }
    attr(newdata, "tname") <- topt[c("id", "tstart", "tstop")]
    attr(newdata, "tcount") <- tcount
    if (length(tevent)) attr(newdata, "tevent") <- tevent
    row.names(newdata) <- NULL
    class(newdata) <- c("data.frame")
    newdata
}
