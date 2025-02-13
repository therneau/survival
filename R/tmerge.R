#
# See the section on tmerge in the methods document for a discusion of
#  how this function is designed.  There are subtleties
#
tmerge <- function(data1, data2, id, ..., tstart, tstop, options) {
    Call <- match.call()
    # The function wants to recognize special keywords in the
    #  arguments, so define a set of functions which will be used to
    #  mark objects
    new <- new.env(parent=parent.frame())
    assign("tdc", function(time, value=NULL, init=NULL, na.rm=NULL) {
        x <- list(time=time, value=value, default= init, na.rm=na.rm); 
        class(x) <- "tdc"; x},
           envir=new)
    # To allow na.rm in cumtdc I would need to change my algorithm, which
    # currently relies on cumsum, plus I can't think of a use case for it.
    assign("cumtdc", function(time, value=NULL, init=NULL) {
        x <- list(time=time, value=value, default= init); 
        class(x) <-"cumtdc"; x},
           envir=new)
    assign("event", function(time, value=NULL, censor=NULL, rep = TRUE) {
        x <- list(time=time, value=value, censor=censor, rep = rep); 
        class(x) <-"event"; x},
           envir=new)
    assign("cumevent", function(time, value=NULL, censor=NULL) {
        x <- list(time=time, value=value, censor=censor); 
        class(x) <-"cumevent"; x},
           envir=new)

    if (missing(data1) || missing(data2) || missing(id)) 
        stop("the data1, data2, and id arguments are required")
    if (!inherits(data1, "data.frame")) stop("data1 must be a data frame")
    tmerge.control <- function(idname="id", tstartname="tstart", tstopname="tstop",
                               delay =0, na.rm=TRUE, tdcstart=NA_real_, ...) {
        extras <- list(...)
        if (length(extras) > 0) 
            stop("unrecognized option(s):", paste(names(extras), collapse=', '))
        if (length(idname) != 1 || make.names(idname) != idname)
            stop("idname option must be a valid variable name")
        if (!is.null(tstartname) && 
            (length(tstartname) !=1 || make.names(tstartname) != tstartname))
            stop("tstart option must be NULL or a valid variable name")
        if (length(tstopname) != 1 || make.names(tstopname) != tstopname)
            stop("tstop option must be a valid variable name") 
        if (length(delay) !=1 || !is.numeric(delay) || delay < 0)
            stop("delay option must be a number >= 0")
        if (length(na.rm) !=1 || ! is.logical(na.rm))
            stop("na.rm option must be TRUE or FALSE")
        if (length(tdcstart) !=1) stop("tdcstart must be a single value")
        list(idname=idname, tstartname=tstartname, tstopname=tstopname, 
             delay=delay, na.rm=na.rm, tdcstart=tdcstart)
    }

    if (!inherits(data1, "tmerge") && !is.null(attr(data1, "tname"))) {
        # old style object that someone saved!
        tm.retain <- list(tname = attr(data1, "tname"),
                          tevent= list(name=attr(data1, "tevent"),
                                       censor= attr(data1, "tcensor")),
                          tdcvar = attr(data1, "tdcvar"),
                          n = nrow(data1))
        attr(data1, "tname") <- attr(data1, "tevent") <- NULL
        attr(data1, "tcensor") <- attr(data1, "tdcvar") <- NULL
        attr(data1, "tm.retain") <- tm.retain
        class(data1) <- c("tmerge", class(data1))
    }
                          
    if (inherits(data1, "tmerge")) {
        tm.retain <- attr(data1, "tm.retain")
        firstcall <- FALSE
        # check out whether the object looks legit:
        #  has someone tinkered with it?  This won't catch everything
        tname <- tm.retain$tname
        tevent <- tm.retain$tevent
        tdcvar <- tm.retain$tdcvar
        if (nrow(data1) != tm.retain$n)
            stop("tmerge object has been modified, size")
        if (any(is.null(match(unlist(tname), names(data1)))) ||
            any(is.null(match(tm.retain$tcdname, names(data1)))) ||
            any(is.null(match(tevent$name, names(data1)))))
            stop("tmerge object has been modified, missing variables")
        for (i in seq(along=tevent$name)) {
            ename <- tevent$name[i]
            if (is.numeric(data1[[ename]])) {
                if (!is.numeric(tevent$censor[[i]]))
                    stop("event variable ", ename, 
                         " no longer matches it's original class")
            }
            else if (is.character(data1[[ename]])) {
                if (!is.character(tevent$censor[[i]]))
                    stop("event variable ", ename, 
                         " no longer matches it's original class")
            }
            else if (is.logical(data1[[ename]])) {
                if (!is.logical(tevent$censor[[i]]))
                    stop("event variable ", ename,
                         " no longer matches it's original class")
            }
            else if (is.factor(data1[[ename]])) {
                if (levels(data1[[ename]])[1] != tevent$censor[[i]])
                    stop("event variable ", ename,
                         " has a new first level")
            }
            else stop("event variable ", ename, " is of an invalid class")
        }
    } else {
        firstcall <- TRUE
        tname <- tevent <- tdcvar <- NULL
        if (is.name(Call[["id"]])) {
            idx <- as.character(Call[["id"]])
            if (missing(options)) options <-list(idname= idx)
            else if (is.null(options$idname)) options$idname <- idx
        }
    }

    if (!missing(options)) {
        if (!is.list(options)) stop("options must be a list")
        if (!is.null(tname)) {
            # If an option name matches one already in tname, don't confuse
            #  the tmerge.control routine with duplicate arguments
            temp <- match(names(options), names(tname), nomatch=0)
            topt <- do.call(tmerge.control, c(options, tname[temp==0]))
            if (any(temp >0)) {
                # A variable name is changing midstream, update the
                # variable names in data1
                varname <- tname[c("idname", "tstartname", "tstopname")]
                temp2 <- match(varname, names(data1))
                names(data1)[temp2] <- varname
            }
        }
        else topt <- do.call(tmerge.control, options)
    }
    else if (length(tname)) topt <- do.call(tmerge.control, tname)
    else  topt <- tmerge.control()

    # id, tstart, tstop are found in data2
    if (missing(id)) stop("the id argument is required")
    if (missing(data1) || missing(data2))
        stop("two data sets are required")
    id <- eval(Call[["id"]], data2, enclos=emptyenv()) #don't find it elsewhere
    if (is.null(id)) stop("id variable not found in data2")
    if (any(is.na(id))) stop("id variable cannot have missing values")

    if (firstcall) {
        if (!missing(tstop)) {
             tstop <-  eval(Call[["tstop"]],  data2)
             if (length(tstop) != length(id))
                 stop("tstop and id must be the same length")
             # The neardate routine will check for legal tstop data type
          }
        if (!missing(tstart)) {
            tstart <- eval(Call[["tstart"]], data2)
            if (length(tstart)==1) tstart <- rep(tstart, length(id))
            if  (length(tstart) != length(id))        
                stop("tstart and id must be the same length")
            if (any(tstart >= tstop))
                stop("tstart must be < tstop")
             }
    }
    else {
        if (!missing(tstart) || !missing(tstop))
            stop("tstart and tstop arguments only apply to the first call")
    }
    # grab the... arguments
    notdot <- c("data1", "data2", "id", "tstart", "tstop", "options")
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
    tcount <- matrix(0L, length(argname), 9)
    dimnames(tcount) <- list(argname, c("early","late", "gap", "within", 
                                        "boundary", "leading", "trailing",
                                        "tied", "missid"))
    tcens <- tevent$censor
    tevent <- tevent$name
    if (is.null(tcens)) tcens <- vector('list', 0)
    newdata <- data1 #make a copy
    if (firstcall) {
        # We don't look for topt$id.  What if the user had id=clinic, but their
        #  starting data set also had a variable named "id".  We want clinic for
        #  this first call.
        idname <- Call[["id"]]
        if (!is.name(idname)) 
            stop("on the first call 'id' must be a single variable name")
     
        # The line below finds tstop and tstart variables in data1
        indx <- match(c(topt$idname, topt$tstartname, topt$tstopname), names(data1), 
                      nomatch=0)
        if (any(indx[1:2]>0) && FALSE) {  # warning currently turned off. Be chatty?
            overwrite <- c(topt$tstartname, topt$tstopname)[indx[2:3]]
            warning("overwriting data1 variables", paste(overwrite, collapse=' '))
            }
        
        temp <- as.character(idname)
        if (!is.na(match(temp, names(data1)))) {
                data1[[topt$idname]] <- data1[[temp]]
                baseid <- data1[[temp]]
                }
        else stop("id variable not found in data1")

        if (any(duplicated(baseid))) 
            stop("for the first call (that establishes the time range) data1 must have no duplicate identifiers")

        if (missing(tstop)) {
            if (length(argclass)==0 || argclass[1] != "event")
                stop("neither a tstop argument nor an initial event argument was found")
            # this is case 2 -- the first time value for each obs sets the range
            last <- !duplicated(id)
            indx2 <- match(unique(id[last]), baseid)
            if (any(is.na(indx2)))
                stop("setting the range, and data2 has id values not in data1")
            if (any(is.na(match(baseid, id))))
                stop("setting the range, and data1 has id values not in data2")
            newdata <- data1[indx2,]
            tstop <- (args[[1]]$time)[last]
        }
        else {
            if (length(baseid)== length(id) && all(baseid == id)) newdata <- data1
            else {  # Note: 'id' is the idlist for data 2
                indx2 <- match(id, baseid)
                if (any(is.na(indx2)))
                    stop("setting the range, and data2 has id values not in data1")
                if (any(is.na(match(baseid, id))))
                    stop("setting the range, and data1 has id values not in data2")
                newdata <- data1[indx2,]
            }
        }
          
        if (any(is.na(tstop))) 
            stop("missing time value, when that variable defines the span")
        if (missing(tstart)) {
            indx <- which(tstop <=0)
            if (length(indx) >0) stop("found an ending time of ", tstop[indx[1]],
                                      ", the default starting time of 0 is invalid")
            tstart <- rep(0, length(tstop))
        }
        if (any(tstart >= tstop)) 
            stop("tstart must be < tstop")
        newdata[[topt$tstartname]] <- tstart
        newdata[[topt$tstopname]] <- tstop
        n <- nrow(newdata)
        if (any(duplicated(id))) {
            # sort by time within id
            indx1 <- match(id, unique(id))
            newdata <- newdata[order(indx1, tstop),]
         }
        temp <- newdata[[topt$idname]]
        if (any(tstart >= tstop)) stop("tstart must be < tstop")
        if (any(newdata$tstop[-n] > newdata$tstart[-1] &
                temp[-n] == temp[-1]))
            stop("first call has created overlapping or duplicated time intervals")
        idmiss <- 0  # the tcount table should have a zero
    }
    else { #not a first call
        idmatch <- match(id, data1[[topt$idname]], nomatch=0)
        if (any(idmatch==0)) idmiss <- sum(idmatch==0)
        else idmiss <- 0
    }
    saveid <- id  # each ii might use a different subset of id, taken from this
    for (ii in seq(along.with=args)) {
        argi <- args[[ii]]
        baseid <- newdata[[topt$idname]]
        dstart <- newdata[[topt$tstartname]]
        dstop  <- newdata[[topt$tstopname]]

        # Toss rows that won't have an impact, before updating the count matrix
        #  - argi$time is missing
        #  - id does not match anyone in data1
        #  - missing yinc value, unless a tdc with na.rm=FALSE
        #  - events or cumevents that are censored
        #  - repeated events, if rep=FALSE
        #  This can make the resulting data set smaller. 
        # 
        etime <- argi$time
        if (idmiss ==0) keep <- rep(TRUE, length(etime))
        else keep <- (idmatch > 0)
        if (length(etime) != length(saveid))
            stop("argument ", argname[ii], " is not the same length as id")
        else keep <- keep & !is.na(etime)

        if (!is.null(argi$value)) {
           if (length(argi$value) != length(saveid))
                stop("argument ", argname[ii], " is not the same length as id")
           if (argclass[ii] != "tdc" ||
               (is.null(argi$na.rm) && topt$na.rm) || argi$na.rm)
               keep <- keep & !is.na(argi$value)
           
           if (!all(keep)) argi$value <- argi$value[keep]
        }
        if (!all(keep)) {
            etime <- etime[keep]
            argi$time <- argi$time[keep]
        }
        id <- saveid[keep]  #defer the event checks until after the sort
        
        # Later steps become easier if we sort the new data by time within
        # id.  Why not have done this once and for all in an earlier step?
        # A user could have different time variables in different tdc or
        #  event calls.  (Not that I expect this, but anything that can be
        #  done eventually will be done.)  The etime part of the sort would
        #  then change from one ii value to the nex.
        # 
        #  The match() is critical when baseid is not in sorted order.  
        #  We made a choice early on to preserve the order of the data found
        #  in data1 of the first call.  id is the id value in data2, baseid that
        #  for data1.  
        indx <- order(match(id, baseid), argi$time)
        id <- id[indx]
        etime <- argi$time[indx]
        if (!is.null(argi$value)) 
            yinc <- argi$value[indx]
        else yinc <- NULL
 
        # For an event or cumevent type we need to know the censor code, in
        #  order to know which rows to keep. At the same time, do transformations
        #  a. if the variable already exists, complain
        #  b. the output of cumevent() is always an integer, the input can be
        #    logical or numeric (preprocessed with floor), or a 2-level factor
        #    (converted to 0-1)
        #  c. an event() will transformed to numeric if it is logical or numeric
        #   with only 2 values, otherwise it becomes a factor. In all cases,
        #   the internal value during processing will be an integer with 0=
        #   censor. 
        #  d. for factors, the optional censor argument can set the code for
        #    censoring
        if (argclass[ii] %in% c("event", "cumevent")) {
            if (is.null(yinc)) yinc <- rep(1L, length(id))
            cval <- NULL
            if (is.logical(yinc)) {
                yinc <- as.integer(yinc)  # becomes 0/1
                if (!is.null(argi$censor) && argi$censor) yinc <- 1L- yinc
                cval <- 0
            }
            if (argclass[ii] == "cumevent" && !is.numeric(yinc))
                stop("argument for cumevent must be numeric or logical")
            if (is.numeric(yinc)) {
                if (any(yinc != floor(yinc)))
                    stop("the argument for an event must be logical, integer or a factor")
                if (any(yinc <0)) 
                    stop("numeric event values must be non-negative")
                if (!is.null(argi$censor)) cval <- argi$censor 
                else cval <- 0
                if (cval !=0 || length(unique(yinc)) >2) {
                    # turn it into a factor
                    if (cval==0 && is.null(argi$censor)) cval <- "censor"
                    ylev <- unique(c(cval, sort(unique(yinc))))
                    yinc <- factor(yinc, ylev)
                }
            }
            if (is.character(yinc)) yinc <- factor(yinc)
            
            if (is.factor(yinc)) {
                if (is.null(cval)) {
                    if (!is.null(argi$censor)) cval <- argi$censor
                    else cval <- "censor"
                }
                newlev <- unique(c(cval, levels(yinc))) # save for later
                    yinc <- as.integer(factor(yinc, newlev)) -1L
                }
            # Done with the recode, yinc will now be an integer with 0
            #  for the censoring code.
            keep <- which(yinc !=0)  # all non-censored rows
            if (!argi$rep) { 
                # remove 'stutterings' from the  data, i.e. a repeated
                #  non-censored event type, as well
                keep.n <- length(keep)
                stutter <- (id[keep[-1]]== id[keep[-keep.n]]) &
                    (yinc[keep[-1]] == yinc[keep[-keep.n]])
                
                keep <- keep[c(TRUE, !stutter)]
            } 
            # remove censored rows from the additions (why add them?)
            #  along with any stutters
            id <- id[keep]
            etime <- etime[keep]
            yinc <- yinc[keep]
        } #done with pre-processing an event

        if (argclass[ii] %in% c("tdc", "cumtdc") && !is.null(yinc) &&
            ((is.null(argi$na.rm) && topt$na.rm) || argi$na.rm)) {
            # Also remove "stuttered" event() or tdc() rows, unless na.rm=FALSE
            #  Removing repeats in the presence of NA is too much work, and
            # this action only saves space. (Truth be told, the savings may 
            # often be trivial with multiple tcd vars; at any given time
            # *something* will change, and thus a new row.
            ny <- length(id)
            stutter <- (id[-1] == id[-ny]) & (yinc[-1] == yinc[-ny])
            if (any(stutter)) {
                keep <- c(TRUE, !stutter)
                id <- id[keep]
                etime <- etime[keep]
                yinc <- yinc[keep]
            }
        }

        # indx1 points to the closest start time in the baseline data (data1)
        #  that is <= etime.  indx2 to the closest end time that is >=etime.
        # If etime falls into a (tstart, tstop) interval, indx1 and indx2
        #   will match
        # If the "delay" argument is set and this event is of type tdc, then
        #   move any etime that is after the entry time for a subject.
        if (topt$delay >0 && argclass[ii] %in% c("tdc", "cumtdc")) {
            mintime <- tapply(dstart, baseid, min)
            index <- match(id, names(mintime))
            etime <- ifelse(etime <= mintime[index], etime, etime+ topt$delay)
        }
        indx1 <- neardate(id, baseid, etime, dstart, best="prior")
        indx2 <- neardate(id, baseid, etime, dstop, best="after")

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

        # count ties.  id is in the users data order, etime sorted within id
        tcount[ii,8] <- sum(tapply(etime, id, function(x) sum(duplicated(x))))
        tcount[ii,9] <- idmiss
        indx4 <- which(itype==4)
        n4 <- length(indx4)
        if (n4 > 0) {
            # we need to eliminate duplicate times within the same id, but
            #  do so without changing the class of etime: it might
            #  be a Date, an integer, a double, ... 
            # Using unique on a data.frame does the trick
            icount <- data.frame(irow= indx1[indx4], etime=etime[indx4])
            icount <- unique(icount)   
            # the icount data frame will be sorted by second column within first
            #  so rle is faster than table
            n.add <- rle(icount$irow)$length # number of rows to add for each id

            # expand the data 
            irep <- rep.int(1L, nrow(newdata))
            erow <- unique(indx1[indx4])   # which rows in newdata to be expanded
            irep[erow] <- 1+ n.add # number of rows in new data
            jrep <- rep(1:nrow(newdata), irep)  #stutter the duplicated rows
            newdata <- newdata[jrep,]  #expand it out
            dstart <- dstart[jrep]
            dstop <-  dstop[jrep]

            #fix up times
            nfix <- length(erow)
            temp <- vector("list", nfix)
            iend <- (cumsum(irep))[irep >1]  #end row of each duplication set
            for (j in 1:nfix) temp[[j]] <-  -(seq(n.add[j] -1, 0)) + iend[j]
            newrows <- unlist(temp)
            
            # this should not be necessary
            #    if (inherits(dstart, "Date"))
            #        icount$etime <- as.Date(icount$etime, origin= "1970-01-01")
            dstart[newrows] <- dstop[newrows-1] <- icount$etime
            newdata[[topt$tstartname]] <- dstart
            newdata[[topt$tstopname]]  <- dstop
            for (ename in tevent) 
                newdata[newrows-1, ename] <- tcens[[ename]][[1]]

            # refresh indices
            baseid <- newdata[[topt$idname]]
            indx1 <- neardate(id, baseid, etime, dstart, best="prior")
            indx2 <- neardate(id, baseid, etime, dstop, best="after")
            subtype[itype==4] <- 1  #all the "insides" are now on a tied edge
            itype[itype==4]   <- 5  
        }

        newvar <- newdata[[argname[ii]]]  # prior value (sequential tmerge calls)
        # add a tdc variable
        if (argclass[ii] %in% c("tdc", "cumtdc")){
            if (argname[[ii]] %in% tevent)
                stop("attempt to turn event variable", argname[[ii]], "into a tdc")
            if (!(argname[[ii]] %in% tdcvar)){
                tdcvar <- c(tdcvar, argname[[ii]])
                if (!is.null(newvar) && argclass[ii] == "tdc") {
                    warning(paste0("replacement of variable '", argname[ii], "'"))
                    newvar <- NULL
                }
            }
        }
        if (argclass[ii] == "tdc") {
            default <- argi$default   # default value
            if (is.null(default)) default <- topt$tdcstart
            else if (length(default) !=1)
                stop("default tdc value must be of length 1")

            # id can be any data type; feed integers to the C routine
            storage.mode(dstart) <- storage.mode(etime) <- "double"  #if time is integer
            uid <- unique(baseid)
            index <- .Call(Ctmerge2, match(baseid, uid), dstart, 
                                       match(id, uid),  etime)

            if (is.null(newvar)) {  # create new variable 
                if (is.null(yinc)) newvar <- ifelse(index==0, 0L, 1L) #add a 0/1 variable
                else {
                    newvar <- yinc[pmax(1L, index)]
                    if (any(index==0)) {
                        if (is.na(default)) is.na(newvar) <- (index==0L)                
                        else {
                            if (is.numeric(newvar)) newvar[index==0L] <- as.numeric(default)
                            else {
                                if (is.factor(newvar)) {
                                    # special case: if default isn't in the set of levels,
                                    #   add it to the levels
                                    if (is.na(match(default, levels(newvar))))
                                        levels(newvar) <- c(levels(newvar), default)
                                }
                                newvar[index== 0L] <- default
                            }
                        }
                    }
                }
            } else if (is.null(yinc)) {
                # Existing variable, no yinc, so update will be 0 or 1
                # Okay only if the current variable is 0/1
                if (all(is.na(newvar) | newvar==0L | newvar==1L))
                    newvar[index!=0L] <- 1L
                else stop("tdc update does not match prior variable type: ", argname[ii])
            } else {
                # attempt to update an existing tdc with new values from a new variable
                # We know how to handle a few special cases, but the basic strategy
                #  is "don't try to be clever".  Remember that class() can return
                #  a vector of length > 1
                if (inherits(newvar, "factor") && (!inherits(yinc, "factor") ||
                                            !identical(levels(newvar), levels(yinc))))
                    stop("tdc update does not match prior factor: ", argname[ii])
                clnew <- class(yinc)
                clold <- class(newvar)
                if (identical(clnew, clold) || (length(clnew)==1 && length(clold)==1  &&
                                         class(newvar) %in% c("integer", "numeric") && 
                                         class(yinc)   %in% c("integer", "numeric")))
                    newvar[index != 0L] <- yinc[index]
                else stop("tdc update does not match prior variable type: ", argname[ii])
            }
            tdcvar <- unique(c(tdcvar, argname[[ii]]))
        }

        # add events
        if (argclass[ii] %in% c("cumtdc", "cumevent")) {
            if (is.null(yinc)) yinc <- rep(1L, length(id))
        }   
        if (argclass[ii] == "cumevent"){
            ykeep <- (yinc !=0)  # ignore the addition of a censoring event
            yinc <- unlist(tapply(yinc, match(id, baseid), cumsum))
        }

        if (argclass[ii] %in% c("event", "cumevent")) {
            if (is.null(yinc)) yinc <- rep(1L, length(id))
            if (!is.null(newvar)) {
                if (!argname[ii] %in% tevent) {
                    #warning(paste0("non-event variable '", argname[ii], "' replaced by an event variable"))
                    newvar <- NULL
                }
                if 
            } else newvar <- rep(0L, nrow(newdata))
         
            if (is.null(yinc)) yinc <- rep(1L, length(id))
            keep <- (subtype==1 | subtype==3) # all other events are thrown away
            if (argclass[ii] == "cumevent") keep <- (keep & ykeep)
            newvar[indx2[keep]] <- yinc[keep]
            # add this into our list of 'this is an event type variable'
            if (!(argname[ii] %in% tevent)) {
                tevent <- c(tevent, argname[[ii]])
                if (is.factor(yinc)) tcens <- c(tcens, list(levels(yinc)))
                else if (is.logical(yinc)) tcens <- c(tcens, list(FALSE))
                else if (is.integer(yinc))   tcens <- c(tcens, list(0L))
                else tcens <- c(tcens, list(0))
                names(tcens) <- tevent
            }
        }
        else if (argclass[ii] == "cumtdc") {  # process a cumtdc variable
            # I don't have a good way to catch the reverse of this user error
            if (argname[[ii]] %in% tevent)
                stop("attempt to turn event variable", argname[[ii]], "into a cumtdc")

            keep <- itype != 2  # changes after the last interval are ignored
            indx <- ifelse(subtype==1, indx1, 
                           ifelse(subtype==3, indx2+1L, indx2))
            
            # we want to pass the right kind of NA to the C code
            default <- argi$default
            if (is.null(default)) default <- as.numeric(topt$tdcstart)
            else {
                if (length(default) != 1) stop("tdc initial value must be of length 1")
                if (!is.numeric(default)) stop("cumtdc initial value must be numeric")
            }       
            if (is.null(newvar)) {  # not overwriting a prior value
                if (is.null(argi$value)) newvar <- rep(0.0, nrow(newdata))
                else newvar <- rep(default, nrow(newdata))
            }
            
            # the increment must be numeric
            if (!is.numeric(newvar)) 
                stop("data and starting value do not agree on data type")
            # id can be any data type; feed integers to the C routine
            storage.mode(yinc) <- storage.mode(dstart) <- "double"
            storage.mode(newvar) <- storage.mode(etime) <- "double"
            newvar <- .Call(Ctmerge, match(baseid, baseid), dstart, newvar, 
                            match(id, baseid)[keep], etime[keep], 
                            yinc[keep], indx[keep])
        }  

        newdata[[argname[ii]]] <- newvar
    }
    tm.retain <- list(tname = topt[c("idname", "tstartname", "tstopname")],
                      n= nrow(newdata))
    if (length(tevent)) {
        tm.retain$tevent <- list(name = tevent, censor=tcens)
        for (i in 1:length(tevent)) {
            if (length(tcens[[i]]) > 1){ # this was a factor
                temp <- newdata[[tevent[i]]]
                temp <- factor(temp, seq(along.with(tcens[[i]])) -1L, tcens[[i]])
                newdata[[tevent[i]]] <- temp
            }
        }
    }
        
    if (length(tdcvar)>0) tm.retain$tdcvar <- tdcvar
    attr(newdata, "tm.retain") <- tm.retain
    attr(newdata, "tcount") <- rbind(attr(data1, "tcount"), tcount)
    attr(newdata, "call") <- Call

    row.names(newdata) <- NULL  #These are a mess; kill them off.
    # Not that it works: R just assigns new row names.
    class(newdata) <- c("tmerge", "data.frame")
    newdata
}
summary.tmerge <- function(object, ...) {
    if (!is.null(cl <- attr(object, "call"))) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
    }

    print(attr(object, "tcount"))
}

# This could be smarter: if you only drop variables that are not known 
# to tmerge then it would be okay.  But I currently like the "touch it
#  and it dies" philosophy
"[.tmerge" <- function(x, ..., drop=TRUE){
    class(x) <- "data.frame"
    attr(x, "tm.retain") <- NULL
    attr(x, "tcount") <- NULL
    attr(x, "call") <- NULL
    NextMethod(x)
    }
