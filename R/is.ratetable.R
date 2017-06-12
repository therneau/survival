is.ratetable <- function(x, verbose=FALSE) {
    att <- attributes(x)
    dlist <- c("dim", "dimnames", "cutpoints")
    if (!verbose) {
	if (!inherits(x, 'ratetable')) return(FALSE)
	if (any(is.na(match(dlist, names(att))))) return(FALSE)
        if (is.null(att$dimid)) att$dimid <- names(att$dimnames)
	nd <- length(att$dim)
	if (length(x) != prod(att$dim)) return(FALSE)
        if (any(att$dimid == "")) return(FALSE)
        
	if (!(is.list(att$dimnames) && is.list(att$cutpoints)))
		 return(FALSE)
	if (length(att$dimnames)!=nd ||
			 length(att$cutpoints)!=nd) return(FALSE)
        # One of 'factor' (old style table) or 'type' (new style) should exist
        if (!is.null(att$factor)) {
            fac <- as.numeric(att$factor)
            if (any(is.na(fac))) return(FALSE)
            if (any(fac <0)) return(FALSE)
            if (length(att$factor)!=nd ) return(FALSE)
            # if any dates use Date cutpoints, all need to
            isDate <- sapply(att$cutpoints, function(x) inherits(x, "Date"))
            if (any(isDate) && any((att$type >2) != isDate)) return(FALSE)
            }
        else if (!is.null(att$type)) {
            if (any(is.na(match(att$type, 1:4)))) return(FALSE)
            fac <- 1*(att$type==1)
            if (length(fac) != nd) return(FALSE)
            }
        else return(FALSE)

        if (length(att$dimid) != nd) return(FALSE)
	for (i in 1:nd) {
	    n <- att$dim[i]
	    if (length(att$dimnames[[i]]) !=n) return(FALSE)
	    if (fac[i]!=1 && length(att$cutpoints[[i]])!=n) return(FALSE)
	    if (fac[i]!=1 && any(order(att$cutpoints[[i]])!= 1:n)) return(FALSE)
	    if (fac[i]==1 && !is.null(att$cutpoints[[i]]))  return(FALSE)
	    if (fac[i]>1 && i<nd) return(FALSE)
	    }
	return(TRUE)
	}

    # verbose return messages, useful for debugging
    msg <- NULL
    if (!inherits(x, 'ratetable')) msg <- c(msg, "wrong class")

    temp <- is.na(match(dlist, names(att)))
    if (any(temp)) 
        msg <- c(msg, paste("missing attribute:", dlist[temp]))
    if (is.null(att$dimid)) att$dimid <- names(att$dimnames)

    # This next error is unlikely, since R itself squawks when you
    #   try to set a wrong dimension.  Ditto with dimnames issues.
    nd <- length(att$dim)
    if (length(x) != prod(att$dim)) 
        msg <- c(msg, 'length of the data does not match prod(dim)')

    if (!is.list(att$dimnames))
	     msg <- c(msg, 'dimnames is not a list')
    if (!is.list(att$cutpoints))
	     msg <- c(msg, 'cutpoints is not a list')

    if (length(att$dimnames)!=nd)
        msg <- c(msg, 'wrong length for dimnames')
    if (length(att$dimid)!=nd)
        msg <- c(msg, "wrong length for dimid, or dimnames do not have names")
    if (any(att$dimid==""))
        msg <- c(msg, "one of the dimnames identifiers is blank")

    if (length(att$cutpoints)!=nd) 
        msg <- c(msg, 'wrong length for cutpoints')

    if (!is.null(att$factor)) {
        fac <- as.numeric(att$factor)
        if (any(is.na(fac))) msg <- c(msg, "illegal 'factor' level of NA")
        if (any(fac <0)) msg <- c(msg, "illegal 'factor' attribute of <0")
        if (length(att$factor)!=nd)
            msg <- c(msg, 'wrong length for factor')
        type <- 1*(fac==1) + 2*(fac==0) + 4*(fac>1)
        }
    else if (!is.null(att$type)) {
        if (any(is.na(match(att$type, 1:4))))
            msg <- c(msg, 'type attribute must be 1, 2, 3, or 4')
        type <- att$type
        if (length(type)!=nd)
            msg <- c(msg, 'wrong length for type attribute')
        # if any dates use Date cutpoints, all need to
        isDate <- sapply(att$cutpoints, function(x) inherits(x, "Date"))
        if (any(isDate) && any((att$type >2) != isDate))
            msg <- c(msg, "all dates must have cutpoints of type Date, or none")
        }
    else msg <- c(msg, "missing the 'type' attribute")

    for (i in 1:nd) {
	n <- att$dim[i]
	if (length(att$dimnames[[i]]) !=n) 
		msg <- c(msg, paste('dimname', i, 'is the wrong length'))

	if (type[i] >1) { #continuous variable
            if (length(att$cutpoints[[i]]) != n)
                msg <- c(msg, paste('wrong length for cutpoints', i))
            else if (any(order(att$cutpoints[[i]])!= 1:n)) 
		msg <- c(msg, paste('unsorted cutpoints for dimension',i))
                }

	if (type[i]==1 && !is.null(att$cutpoints[[i]]))  
		msg <- c(msg, paste('type[', i, 
                                    '] is 1; cutpoint should be null'))
        # This message only applies to old style rate table
	if (!is.null(att$fac) && type[i]==4 && i<nd) 
		msg <- c(msg, 'only the last dimension can be interpolated')
	}
    if (length(msg)==0) TRUE
    else msg
    }
