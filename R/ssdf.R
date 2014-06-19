#
# Helper routine to create start-stop data frames
#
ssdf <- function(d1, d2, id, tstart, time, event, zero.one,
                 startname="tstart", stopname="tstop") {
    if (!inherits(d1, "data.frame")) 
        stop ("d1 must be a data.frame or ssdf object")

    # This action is done many times to match a variable in a data frame
    ssmatch(df, miss, arg, defaultname, exclude, label, nullok=FALSE) {
        # if the arg is missing, try the default name first, then the
        #  first variable in the data frame df that is not on the exclude list
        dnames <- names(df)
        fallback <- (dnames[!(dnames %in% exclude)])[1]
        vname <- defaultname
        if (miss) {
            if (defaultname %in% dnames) {
                value <- df[[defaultname]]
                exclude <- c(exclude, defaultname)
            }
            else if (length(fallback) ==0) {
                if (nullok) value <- NULL
                else stop(paste("no value found for", label))
            }
            else {
                value <- df[[fallback]]
                exclude <- c(exclude, fallback)
                vname <- fallback
            }
        }
        else {
            # if it is a variable name. use it
            if (is.character(arg) && length(arg)==1) {
                value <- df[[arg]]
                if (is.null(value)) stop(paste("variable", arg, "not found"))
                else exclude <- c(exclude, arg)
                vname <- arg
            }
            # or allow for a direct value
            else {
                value <- arg
                if (length(arg) != nrow(df))
                stop(paste("wrong length for", label))
            }
        }
        
        list(value=value, exclude=exclude, vname=vname)
    }

    # If argument d1 is an ssdf object, then we are building onto it.
    #  Otherwise we need to find out the names of the key variables and
    #  refactor it using those names

    if (inherits(d1, "ssdf"))  rname <- attr(d1, "rname")
    else {
        rname <- vector("character", 4)
        temp <- ssmatch(d1, missing(tstart), tstart, startname, 
                        names(d1), "tstart", nullok=TRUE)
        start <- temp$value
        if (is.null(start)) start <- rep(0., nrow(d1))
        rname[2] <- temp$vname        
        if (!is.numeric(tstart)) stop("tstart must be numeric")
        
        temp <- ssmatch(d1, missing(id), id, "", NULL, "id")
        id <- temp$value
        if (is.null(id)) stop("no value found for id")
        rname[1] <- temp$vname
 
        temp <- ssmatch(d1, missing(time), time, stopname, temp$exclude, 
                        "time")
        tstop <- temp$value
        rname[3] <- temp$vname
        if (!is.numeric(tstop)) stop("tstop must be numeric")
 
        temp <- ssmatch(d1, missing(event), event, "", temp$exclude, "event")
        event <- temp$event
        rname[4] <- temp$vname
        if (rname[4]=="") stop("no name available for the event variable")
        if (is.logical(event)) event <- ifelse(event, 1L, 0L)
        if (!is.numeric(event) || any(event!=1 & event!=2))
            stop("invalid event code")
        
        remain <- names(d1)
        remain <- remain[!(remain %in% temp$exclude)]  #variables remaining in d1
        
        d1 <- data.frame(id, tstart, tstop, event, d1[,remain])
        names(d1) <- c(rname, remain)
        d1 <- d1[order(d1[,1], d1[,2]),]  #order by time within id
        if (any(d1[,2] >= d1[,3])) stop("start times must be < stop times")

        n <- nrow(d1)
        if (any(d1[-1,3] < d1[-n,3] & d1[-1,1]==d1[-n,1]))
            stop("overlapping time intervals for a subject")
    }
        
    # Four cases: d1 is or is not an ssdf object, d2 is or is not present
    if (!inherits(d1, "ssdf")) {
        if (missing(df2)) {
            #an initial call that first sets up an ssdf object
            attr(d1, "rname") <- rname
            class(d1) <- "ssdf"
            return(d1)
            }
        else stop("d2 is present and d1 is not an ssdf object")
    }

    # Moderately simple case: adding a single "becomes present" variable
    # d2 is missing, and only id, tstop and zero.one are present
    if (missing(d2)) {
        if (missing(zero.one) || missing(tstop) || missing(id) ||
            !missing(event) || !missing(tstart))
            stop(paste("d1 is an ssdf and d2 is missing;"
                      "only id, tstop and zero.one arguments should be present"))
        if (!is.character(zero.one) || length(zero.one) != 1)
            stop("the zero.one argument should specify a variable name")
        if (zero.one %in% names(d1))
            stop("zero.one variable name already present")
        if (any(duplicated(id)))
            stop("id variable has duplicates")
       
        # Check that all id values in d2 are found in d1, if not print a sample
        #  of the mismatches
        idlist <- unique(df1[[1]])
        indx <- match(id, idlist)
        if (any(is.na(indx))) {
            bad <- which(is.na(indx))
            if (length(bad) > 5) {
                bad <- bad[1:5]
                stop(paste("id values not matched in d1:",
                           paste(bad, collapse=" "), "..."))
            }
            else stop(paste("id values not matched in d1:",
                             paste(bad, collapse=" ")))
        }
        
        # 
        sindex <- .Call("ssdf1", match(df1[[1]], idlist), indx,
                        df1[[2]],  df1[[3]], tstop)
        newdf <- df1[sindex$index,]
        for (i in rlist[-(1:3)]) { #for all event variables
            newdf[i, sindex$newend] <- 0
        }
        newdf[[zero.one]] <- sindex$new
        attr(newdf, "rname") <- rname
        class(newdf) <- "ssdf"
        return(newdf)
    }

    #
    # Now for the most general case: d1 is an ssdf object and
    #  d2 is a new data frame
    # First match variable names
    if (!inherits(d2)) stop("d2 must be a data frame")
    if (nrow(d2) ===0) stop("d2 has 0 rows")
    d2name <- names(d2)
    temp <- ssmatch(d2, missing(id), id, rname[1], exclude=d2name, "id")
    id <- temp$value

    temp <- ssmatch(d2, missing(tstart), tstart, rname[2], exclude=d2name,
                      "tstart", nullok=TRUE)
    tstart <- temp$value
    if (!is.null(tstart) && !is.numeric(tstart)) stop("tstart must be numeric")
    
    temp <- ssmatch(d2, missing(time), time, rname[3], exclude=d2name,
                     "time")
    tstop <- temp$value
    if (!is.numeric(tstop)) stop("tstop must be numeric")

    temp <- ssmatch(d2, missing(event), event, "", exclude=d2name, 
                     "event", nullok=TRUE)
    event <- temp$value
    if (!is.null(event)) {
        if (is.logical(event)) event <- ifelse(event, 1, 0)
        ename <- temp$vname
    }
    
    # Check that all id values in d2 are found in d1, if not print a sample
    #  of the mismatches
    idlist <- unique(df1[[1]])
    indx <- match(id, idlist)
    if (any(is.na(indx))) {
        bad <- which(is.na(indx))
        if (length(bad) > 5) {
            bad <- bad[1:5]
            stop(paste("id values not matched in d1:",
                       paste(bad, collapse=" "), "..."))
        }
        else stop(paste("id values not matched in d1:",
                        paste(bad, collapse=" ")))
    }

    # if tstart was missing, create proper start, stop pairs
    #  time is assumed to be when values in d2 were measured, so is
    #  the start of the interval
    n <- nrow(d2)
    if (is.null(tstart)) {
        tstart <- tstop
        tstop <- NULL
    }

