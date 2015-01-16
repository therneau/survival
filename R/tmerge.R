tmerge <- function(data, id, ..., tname) {
    Call <- match.call()
    # The function wants to recognize special keywords in the
    #  arguments, so define a set of functions which will be used to
    #  mark objects
    new <- new.env(parent=parent.frame())
    assign("count", function(...) {x <- list(...); class(x) <- "count"; x},
           envir=new)
    assign("cumcount", function(...) {x <- list(...); class(x) <-"cumcount"; x},
           envir=new)
    assign("event", function(...) {x <- list(...); class(x) <-"event"; x},
           envir=new)
    assign("cumevent", function(...) {x <- list(...); class(x) <-"cumevent"; x},
           envir=new)
 
    if (missing(data)) stop("a data argument is required")
    dname <- names(data)

    if (is.null(attr(data, "tname"))) {
        # The identites of the key variables are not yet known
        if (missing(tname))
            stop("the key variables must be identified")
        if (is.null(names(tname)) || !is.character(tname))
            stop("tname must be a named list of character strings")
        indx <- match(names(tname),  c("id", "start", "end", "status"))
        if (any(is.na(indx)))
            stop("unrecognized element in 'tname' argument")
        if (all(indx!=1)) stop("an id variable must be identified")
        if (all(indx!=3)) stop("a time variable for the intervals must be identified")
        dindx <- match(tname, dname)
        if (any(is.na(indx)))
            stop(paste("variable", tname[is.na(indx)], "not found in data"))

        if (all(indx!=2)) {
            tname <- c(tname, start="tstart")
            data <- cbind(data, tstart =0)
        }

        if (all (indx !=4)) # no "status" variable
            tname <- tname[c("id", "start", "end")]
        else tname <- tname[c("id", "start", "end", "status")]

        dindx <- match(tname, names(data))
        if (is.na(dindx[2])) {
            # No start variable, initialize it to zeros
            data[[tname[2]]] <- rep(0., nrow(data))
        }
        dindx <- match(tname, names(data))
        if (any(is.na(dindx))) {
            stop(names(tname[is.na(dindx)]), "variable not found")
        }

        if (any(data[[tname[2]]] >= data[[tname[3]]]))
            stop("start time >= end time for at least one observation")
    }
    else tname <- attr(data, 'tname') # data is a prior tmerge object

    # The ... arguments, but don't evaluate them yet
    unused <- Call[is.na(match(names(Call), c("data", "id", "tname")))]
    if (length(unused) ==1) {
        # An initial call, usually, with nothing to add.  Rather than indent
        #  the entire remainder of the code put a return here.
        attr(data, "tname") <- tname
        attr(data, "tcount") <- NULL  #remove tcount if it exists
        return(data)
    }

    # Now for the actual work, adding a variable into the mix
    # Each of the newvars should be a time variable which fits into the
    # time scale of the starter data set
    unused[[1]] <- as.name("list")  # The as-yet unused arguments
    args <- eval(unused, envir=new)
    argclass <- sapply(args, function(x) (class(x))[1])
    argname <- names(args)
    if (any(argname== "")) stop("all additional argments must have a name")
       
    check <- match(argclass, c("count", "cumcount", "event", "cumevent"))
    if (any(is.na(check)))
        stop(paste("argument(s)", argname[is.na(check)], 
                   "not a recognized type"))

    dname <- match(tname, names(data))
    names(dname) <- names(tname)
    if (any(is.na(dname))) 
        stop("data set does not match its own tname attribute")
                   
    indx <- match(id, data[[dname["id"]]])
    if (any(is.na(indx))) stop("new data has subjects not in the base data set")


    # The tcount matrix is useful for debugging
    tcount <- matrix(0L, length(argname), 7)
    dimnames(tcount) <- list(argname, c("early","late", "gap", "within", 
                                        "tied edge", "front edge", "back edge"))

    newdata <- data
    row.names(newdata) <- NULL

    for (i in 1:length(args)) {
        baseid <- newdata[[dname["id"]]]
        dstart <- newdata[[dname["start"]]]
        dstop  <- newdata[[dname["end"]]]

        # if an event time is missing then skip that obs
        etime <- args[[i]][[1]]
        keep <- !is.na(etime)
        etime <- etime[keep]
        id <- id[keep]

        indx1 <- neardate(id, baseid, etime, dstart, best="prior")
        indx2 <- neardate(id, baseid, etime, dstart, best="after")
browser()
        # The event times fall into one of 5 categories
        #   1. Before the first interval
        #   2. After the last interval
        #   3. Outside any interval but with time span, i.e, it falls into
        #       a gap in follow-up
        #   4. Strictly inside an interval (don't touch either end)
        #   5. Inside an interval, but touching.
        itype <- ifelse(is.na(indx1), 1,
                        ifelse(is.na(indx2), 2, 
                               ifelse(indx2 > indx1, 3,
                                      ifelse(etime== dstart[indx1] | 
                                             etime== dstop[indx2], 5, 4))))

        # Subdivide the events that touch on a boundary
        #   Common: e.g. the subject has time intervals of
        #      (a,b] and (b,c] with a new count at b.
        #  Start: an interval (a,b], new count at a, subject not at risk at a-0
        #  End: similar to start
        #  
        subtype <- ifelse(itype!=5, 0, 
                          ifelse(indx1 == indx2+1, 1,
                                 ifelse(etime==dstart[indx1], 2, 3)))
        tcount[i,] <- table(factor(itype+subtype, levels=c(1:4, 6:8)))

        if (length(args[[i]]) >1) {
            if (length(args[[i]]) > 2) 
                stop("too many variables in an", argclass[i], "call")
            if (diff(sapply(args[[i]], length)) !=0)
                stop("different lengths in an", argclass[i], "call")
            yinc <- args[[i]][[2]]
            if (!is.numeric(istep)) 
                stop("non numeric increment in an", argclass[i], "call")
            mfun <- function(x, grp) {
                temp <- tapply(yinc[grp], x[grp], sum)
                ifelse(is.na(temp), 0, temp)
            }
        }
        else {
            mfun <- function(x, grp) table(x[grp])
            yinc <- rep(1L, length(etime))     # each counts as 1
        }
   
        # Now fold it in
        increment <- rep(0, nrow(newdata))
        eflag <- (argclass[i] %in% c("event", "cumevent")) # 'event' type
        if (eflag) {
            if (any(subtype==1)) { #subtype 1 events to the earlier interval
                count1 <- myfun(indx2, subtype==1)
                itemp <- as.numeric(names(count1))
                increment[itemp] <- increment[itemp] + c(count1)
            }
            if (any(subtype==3)) { # go to the matched interval
                count3 <- myfun(indx2, subtype==3)
                itemp <- as.numeric(names(count3))
                increment[itemp] <- increment[itemp] + c(count3)
            }
            #subtype 2 and type 3 events are ignored
        }
        else {
            if (any(itype==1)) {
                count <- myfun(id, itype==1)  #there might be multiples
                itemp <- match(names(count), baseid)
                increment[itemp] <- increment[itemp] + c(count)
            }

            if (any(subtype==1)) { #subtype 1 events to the later interval
                count1 <- myfun(indx1, subtype==1)
                itemp <- as.numeric(names(count1))
                increment[itemp] <- increment[itemp] + c(count1)
            }
  
            if (any(subtype==2)) { # subtype 2 go to the matched interval
                count2 <- myfun(indx1, subtype==2)
                itemp <- as.numeric(names(count2))
                increment[itemp] <- increment[itemp] + c(count2)
            }
            # subtype 3 are ignored (as are type 3 events)
        } 


        # Type 4 forces us to split rows of the data set
        #  A single subject might have muliple jumps on a single day,
        #   or a single day have a jump for multiple subjects, so we need
        #   to count up by day and subject
        #  The "increment" variable has all the correct values at this
        #   point, we just have to assign all parts to the righ row
        if (any(itype==4)) {
            indx4 <- which(itype==4)
            n4 <- length(indx4)
            # The simple approach at this point starts with a table
            #  of subject by time, showing who gets split at what times.
            #  But such a table can be sparse and HUGE.
            # Instead do a sequential approach.  We know the data is
            #  sorted by time within each subject.  The set of pairs
            #  cbind(id, etime)[itype==4,] tells me the id/etime pairs to
            #  expand.  Create vector indices that point to those rows.
            # "first" will be true for first of each unique (id, etime) set
            #  etime[first] is the set of new times, 
            #  indx1[first] = the obs number where it inserts
            #  each of these represents a new value to insert
            firstid <- c(TRUE, diff(indx1[indx4]) !=0)
            first <- (firstid | c(TRUE, diff(etime[indx4])!=0))
            irows<- (indx1[indx4])[first] #which rows in base to expand
            newrows <- c(table(irows))  #single interval may have >1 insert
            irows <- unique(irows)
            # Next line a bit tricky.  I want to know the total count
            #  for each unique id/etime pair. If someone had a prior interval
            #  of (3,10] say, with 2 new obs on day 5 and 1 on day 7, then
            #  newcount will have two values, one for day 5 and one for 7
            #  cumsum(yinc[indx4]) is a running sum of all I want to count.  
            #  The right "diff" gives the correct totals.
            newcount <- diff((c(0, cumsum(yinc[indx4])))[c(first, TRUE)]) 
            etemp <- rep(1L, nrow(newdata))
            etemp[irows] <- 1+ newrows
            newindx <- rep(1:nrow(newdata), etemp)  
            newdata <- newdata[newindx,]
            increment <- increment[newindx]  #expand increment

            #
            # Now fix up the new data set
            #  For each subject the new set of start times = c(old,newtimes)
            #                                stop = c(newtimes, old)
            #                                status=c(rep(0,newtimes), old)
            #                                increment=c(old, count)
            #  rindx is the index of rows for each changed interval
            rindx <- which(diff(newindx)==0)  #the added rows
            newtimes <- (etime[indx4])[first]
            newdata[rindx,   dname["end"]] <- newtimes
            newdata[rindx+1, dname["start"]] <- newtimes
            if (!is.na(dname["status"]))
                newdata[rindx,   dname["status"]] <- 0
            if (eflag) increment[rindx] <- newcount
            else       increment[rindx+1] <- newcount
        }

        # Now update the count variable
        if (argclass[i] %in% c("cumcount", "cumevent")) {
            # Cumulative within person
            temp <- cumsum(increment)
            idname <- dname["id"]
            indx <- match(newdata[[idname]], newdata[[idname]])
            newdata[[argname[i]]] <- temp + increment[indx] - temp[indx]
        }
        else newdata[[argname[i]]] <- increment
    }
    
    attr(newdata, "tcount") <- tcount
    row.names(newdata) <- NULL
    newdata
}

                         
    
