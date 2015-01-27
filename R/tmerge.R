tmerge <- function(data1, data2, id, ..., tstart) {
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
    if (missing(data1)) stop("the data1 argument is required")
    if (missing(id)) stop("the id argument is required")
                          
    newdata <- data1
    row.names(newdata) <- NULL

    # The id variable is found in data2, if present
    if (!missing(data2)) {
        id <- eval(Call[["id"]], data2)
        if (!missing(tstart)) tstart <- eval(Call[["tstart"]], data2)
    }

    temp <- attr(data, "tname")
    if (is.null(temp)) { 
        #first call to the function.  Add variables to the base data set
        tname <- c(id="id", tstart="tstart", tstop="tstop") 
        if (length(id) != nrow(newdata)) stop("wrong length for id")
        newdata$id <- id
        if (missing(tstart)) {
            newdata$tstart <- newdata$tstop <- rep(0, nrow(newdata))
        }
        else {
            if (length(tstart) != nrow(newdata))
                stop("wrong length for tstart")
            newdata$tstart <- newdata$tstop <- tstart
        }
    }
    else {
        indx <- match(tname, names(data1))
        if (any(is.na(indx))) 
            stop("data1 does not match its own tname attribute")
    }

    # The ... arguments, but don't evaluate them yet
    unused <- Call[is.na(match(names(Call), 
                               c("data1", "data2", "id", "tstart")))]
    if (length(unused) ==1) {
        # An initial call that added nothing.  Very odd.
        attr(data, "tname") <- tname
        attr(data, "tcount") <- NULL  #remove tcount if it exists
        return(data)
    }

    # Now for the actual work, adding a variable into the mix
    # Each of the newvars should be a time variable which fits into the
    # time scale of the starter data set
    unused[[1]] <- as.name("list")  # The as-yet unused arguments
    if (missing(data2)) args <- eval(unused, envir=new)
    else  args <- eval(unused, data2, enclos=new)

    argclass <- sapply(args, function(x) (class(x))[1])
    argname <- names(args)
    if (any(argname== "")) stop("all additional argments must have a name")
       
    check <- match(argclass, c("count", "cumcount", "event", "cumevent"))
    if (any(is.na(check)))
        stop(paste("argument(s)", argname[is.na(check)], 
                   "not a recognized type"))

    indx <- match(id, newdata[[tname["id"]]])
    if (any(is.na(indx))) stop("new data has subjects not in the base data set")

    # The tcount matrix is useful for debugging
    tcount <- matrix(0L, length(argname), 7)
    dimnames(tcount) <- list(argname, c("early","late", "gap", "within", 
                                        "tied edge", "front edge", "back edge"))
    tevent <- attr(data, "tevent") # event type variables
    
    for (i in 1:length(args)) {
        baseid <- newdata[[tname["id"]]]
        dstart <- newdata[[tname["tstart"]]]
        dstop  <- newdata[[tname["tstop"]]]

        # if an event time is missing then skip that obs
        etime <- args[[i]][[1]]
        keep <- !is.na(etime)
        etime <- etime[keep]
        id <- id[keep]

        # indx1 points to the closest start time in the baseline data (data1)
        #  that is <= etime.  indx2 to the closest end time that is >=etime.
        # If etime falls into a (tstart, tstop) interval, these will bracket
        #  it.
        indx1 <- neardate(id, baseid, etime, dstart, best="prior")
        indx2 <- neardate(id, baseid, etime, dstop , best="after")

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
        #   Common: e.g. the subject has time intervals of
        #      (a,b] and (b,c] with a new count at b.
        #  1: intervals of (a,b] (b,d], new count at b
        #  2: intervals of (a,b] (c,d] with c>b, new count at c
        #  3: intervals of (a,b] (c,d] with c>b, new count at b
        #
        subtype <- ifelse(itype!=5, 0, 
                          ifelse(indx1 == indx2+1, 1,
                                 ifelse(etime==dstart[indx1], 2, 3)))
        tcount[i,] <- table(factor(itype+subtype, levels=c(1:4, 6:8)))

        # Look to see if this term has one or two arguments.  If one arg
        #  then the increment is 1, else it is the second arg.  The myfun()
        #  function returns the totals by unique subject/time pair
        #
        if (length(args[[i]]) >1) {
            if (length(args[[i]]) > 2) 
                stop("too many variables in an", argclass[i], "call")
            if (diff(sapply(args[[i]], length)) !=0)
                stop("different lengths in an", argclass[i], "call")
            yinc <- args[[i]][[2]]
            if (!is.numeric(yinc)) 
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
 
        # Now fold it in
        # The "increment" variable contains the "jump" at each time.
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
            if (any(itype==2)) { #extend the window of observation
                count2 <- myfun(indx1, itype==2)
                itemp <- as.numeric(names(count2))
                increment[itemp] <- increment[itemp] + c(count2)
                newdata[itype==2, tname[["tstop"]]] <- etime[itype==2]
            }
            #subtype 2, type 1 and type 3 events are ignored
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
            newdata[rindx,   tname["tstop"]] <- newtimes
            newdata[rindx+1, tname["tstart"]] <- newtimes
            for (i in tevent)
                newdata[rindx, i] <- 0
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

        if (eflag) tevent <- c(tevent, argname[i]) #this is an event type var
    }
    
    #
    # Clean up rows that are redundant, usually a (0,0) row followed by
    #   a (0,t) row, which happens with the first call.
    #
browser()
    time1 <- newdata[[tname[2]]]
    time2 <- newdata[[tname[3]]]
    ties <- (time1==time2)
    if (any(ties)) {
        id <- newdata[[tname[1]]]
        n <- nrow(newdata)
        touch <- ((id[-1] == id[-n]) & time1[-1] == time2[-n])
        toss <- c(ties[-n] & touch, FALSE)
        if (any(toss)) newdata <- newdata[!toss,]
    }
    
    attr(newdata, "tcount") <- tcount
    if (length(tevent)) attr(newdata, "tevent") <- tevent
    row.names(newdata) <- NULL
    newdata
}

                         
    
