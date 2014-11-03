tmerge <- function(data, id, start="start", end, status,  ...) {
    Call <- match.call()
    # The function wants to recognize some special keywords in the
    #  arguments, so define a set of functions which will be used to
    #  mark objects
    new <- new.env(parent=parent.frame())
    assign("count", function(x) {class(x) <- c("count", class(x)); x}, 
           envir=new)
    assign("cumcount", function(x) {class(x) <- c("cumcount", class(x)); x}
           , envir=new)
    assign("event",  function(x) {class(x) <- c("event", class(x)); x}, 
           envir=new)
    assign("cumevent", function(x) {class(x) <- c("cumevent", class(x)); x}
           , envir=new)
    assign("variable", 
           function(...) {x <- list(...); class(x) <- "variable"; x},
           envir= new)
 
    if (missing(data)) stop("a data argument is required")
    dname <- names(data)
    if (missing(id)) stop("the id argument is required")

    if (is.null(attr(data, "tname"))) {
        # The first call has to set up the identities of the 
        #  variables
        # Id variable first
        tname <- rep("", 4)
        names(tname) <- c("id", "start", "end", "status")

        if (length(id)==1 && is.character(id)) {
            j <- match(id, names(data))
            if (missing(j)) stop("id variable name not found in data")
            tname[1] <- id
            id <- data[[id]]
        }
        else stop("on the initial call 'id' should be a variable name")

        if (missing(end)) stop("end argument not supplied")
        if (length(end)==1 && is.character(end)) {
            j <- match(end, names(data))
            if (missing(j)) stop("end variable name not found in data")
            tname[3] <- end
        }
        else stop("on the initial call 'end' should be a variable name")

        if (missing(status)) stop("status argument not supplied")
        if (length(status)==1 && is.character(status)) {
            j <- match(status, names(data))
            if (missing(j)) stop("status variable name not found in data")
            tname[4] <- status
        }
        else stop("on the initial call 'status' should be a variable name")
        temp <- data[[status]]
        if (is.logical(temp)) {
            temp <- ifelse(temp, 0,1)
            data[[status]] <- temp
        }
        if (!is.numeric(temp)) 
            stop("status variable must be numeric or logical")
        if (any(temp!=0 & temp!=1)) 
            stop("numeric status values must be 0 or 1")
        
        # The start variable is allowed to be missing, in which case we 
        #   add one
        if (length(start)==1 && is.character(start)) {
            j <- match(start, names(data))
            if (is.na(j)) data <- cbind(data, start=0)
            tname[2] <- start
        }
        else stop("on the initial call 'start' should be a variable name")

        unused <- Call[is.na(match(names(Call), c("data", "id", "start",
                                                  "end", "status")))]
    }
    else {
        tname <- attr(data, 'tname') # data is a prior tmerge object
        unused <- Call[is.na(match(names(Call), c("data", "id")))]
    }

    if (length(unused) <=1) {
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
    if (any(argname== "")) stop("all argments must have a name")
       
    check <- match(argclass, c("count", "cumcount", "event", 
                               "cumevent", "variable"))
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
        if (argclass[i] == "variable")  etime <- args[[i]][[1]]
        else etime <- args[[i]]
        keep <- !is.na(etime)
        etime <- etime[keep]
        class(etime) <- class(etime)[-1] #throw away my fake class
        id <- id[keep]

        indx1 <- neardate(etime, dstart, id, baseid, best="prior")
        indx2 <- neardate(etime, dstop,  id, baseid, best="after")
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

        if (argclass[i] == "variable") {
            stop("'variable' code not yet finished")
        }
        else {
            increment <- rep(0, nrow(newdata))
            eflag <- (argclass[i] %in% c("event", "cumevent")) # 'event' type
            if (eflag) {
                if (any(subtype==1)) { #subtype 1 events to the earlier interval
                    count1 <- table(indx2[subtype==1])
                    itemp <- as.numeric(names(count1))
                    increment[itemp] <- increment[itemp] + c(count1)
                }
                if (any(subtype==3)) { # go to the matched interval
                    count3 <- table(indx2[subtype==3])
                    itemp <- as.numeric(names(count3))
                    increment[itemp] <- increment[itemp] + c(count3)
                }
                #subtype 2 and type 3 events are ignored
            }
            else {
                if (any(itype==1)) {
                    count <- table(id[itype==1])  #there might be multiples
                    itemp <- match(names(count), baseid)
                    increment[itemp] <- increment[itemp] + c(count)
                }

                if (any(subtype==1)) { #subtype 1 events to the later interval
                    count1 <- table(indx1[subtype==1])
                    itemp <- as.numeric(names(count1))
                    increment[itemp] <- increment[itemp] + c(count1)
                }
  
                if (any(subtype==2)) { # subtype 2 go to the matched interval
                    count2 <- table(indx1[subtype==2])
                    itemp <- as.numeric(names(count2))
                    increment[itemp] <- increment[itemp] + c(count2)
                }
                # subtype 3 are ignored (as are type 3 events)
            } 


            # Type 4 forces us to split rows of the data set
            #  A single subject might have muliple jumps on a single day,
            #   or a single day a jump for multiple subjects, so we need
            #   to count up by day and subject
            if (any(itype==4)) {
                indx4 <- which(itype==4)
                n4 <- length(indx4)
                # first will be true for first of each unique (id, etime) set
                #  etime is the set of new times, indx1=the obs it went into
                #  so each of these represents a new value to insert
                firstid <- c(T, diff(indx1[indx4]) !=0)
                first <- (firstid | c(TRUE, diff(etime[indx4])!=0))
                irows<- (indx1[indx4])[first] #which rows in base to expand
                newrows <- c(table(irows))  #single interval may have >1 insert
                irows <- unique(irows)
                newcount <- diff(c(which(first), 1+n4)) #multiple events, 1 day
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
    }
    attr(newdata, "tcount") <- tcount
    row.names(newdata) <- NULL
    newdata
}

                         
    
