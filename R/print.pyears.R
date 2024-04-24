# Automatically generated from the noweb directory
print.pyears <- function(x, ...) {
    if (!is.null(cl<- x$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
        }

    if (is.null(x$data)) {
        if (!is.null(x$event))
            cat("Total number of events:", format(sum(x$event)), "\n")
        cat (   "Total number of person-years tabulated:", 
             format(sum(x$pyears)),
             "\nTotal number of person-years off table:",
             format(x$offtable), "\n")
        }
    else {
        if (!is.null(x$data$event))
            cat("Total number of events:", format(sum(x$data$event)), "\n")
        cat (   "Total number of person-years tabulated:", 
             format(sum(x$data$pyears)),
             "\nTotal number of person-years off table:",
             format(x$offtable), "\n")
        }
    if (!is.null(x$summary)) {
        cat("Matches to the chosen rate table:\n  ", 
            x$summary)
        }
    cat("Observations in the data set:", x$observations, "\n")
    if (!is.null(x$na.action))
      cat("  (", naprint(x$na.action), ")\n", sep='')
    cat("\n")
    invisible(x)
}
summary.pyears <- function(object, header=TRUE, call=header,
                           n= TRUE, event=TRUE, pyears=TRUE,
                           expected = TRUE, rate = FALSE, rr = expected,
                           ci.r = FALSE, ci.rr = FALSE, totals=FALSE,
                           legend=TRUE, vline = FALSE, vertical = TRUE,
                           nastring=".", conf.level=0.95, 
                           scale= 1, ...) {
    # Usual checks
    if (!inherits(object, "pyears")) 
        stop("input must be a pyears object")
    temp <- c(is.logical(header), is.logical(call), is.logical(n),
              is.logical(event) , is.logical(pyears), is.logical(expected),
              is.logical(rate), is.logical(ci.r), is.logical(rr),
              is.logical(ci.rr), is.logical(vline), is.logical(vertical),
              is.logical(legend), is.logical(totals))
    tname <- c("header", "call", "n", "event", "pyears", "expected",
               "rate", "ci.r", "rr", "ci.rr", "vline", "vertical", 
               "legend", "totals")
    if (any(!temp) || length(temp) != 14 || any(is.na(temp))) {
        stop("the ", paste(tname[!temp], collapse=", "), 
             "argument(s) must be single logical values")
    }
    if (!is.numeric(conf.level) || conf.level <=0 || conf.level >=1 |
        length(conf.level) > 1 || is.na(conf.level) > 1)
        stop("conf.level must be a single numeric between 0 and 1")
    if (is.na(scale) || !is.numeric(scale) || length(scale) !=1 || scale <=0)
        stop("scale must be a value > 0")
    
    vname <- attr(terms(object), "term.labels")  #variable names

    if (!is.null(object$data)) {
        # Extra work: restore the tables which had been unpacked into a df
        #  All of the categories are factors in this case
        tdata <- object$data[vname]  # the conditioning variables
        dname <- lapply(tdata, function(x) {
            if (is.factor(x)) levels(x) else sort(unique(x))}) # dimnames
        dd  <-   sapply(dname, length)                # dim of arrays
        index <- tapply(tdata[,1], tdata) 
        restore <- c('n', 'event', 'pyears', 'expected') #do these, if present
        restore <- restore[restore %in% names(object$data)] 
        new   <- lapply(object$data[restore],
                        function(x) {
                            temp <- array(0L, dim=dd, dimnames=dname)
                            temp[index] <- x
                            temp} )
        object <- c(object, new)
    }

    if (is.null(object$expected)) {
        expected <- FALSE
        rr <- FALSE
        ci.rr <- FALSE
    }
    if (is.null(object$event)) {
        event <- FALSE
        rate <- FALSE
        ci.r <- FALSE
        rr <- FALSE
        ci.rr <- FALSE
    }
        
    # print out the front matter
    if (call && !is.null(object$call)) {
        cat("Call: ") 
        dput(object$call) 
        cat("\n")
    }
    if (header) {
        cat("number of observations =", object$observations)
        if (length(object$omit))
            cat("  (", naprint(object$omit), ")\n", sep="")
        else cat("\n")
        if (object$offtable > 0)
            cat(" Total time lost (off table)", format(object$offtable), "\n")
        cat("\n")
    }
    
    # Add in totals if requested
    if (totals) {
        # if the pyear object was based on any time dependent cuts, then
        #  the "n" component cannot be totaled up.
        tcut <- if (is.null(object$tcut)) TRUE else object$tcut
        object$n <- pytot(object$n, na=tcut)
        object$pyears <- pytot(object$pyears)
        if (event) object$event <- pytot(object$event)
        if (expected) object$expected <- pytot(object$expected)
    }
        
    dd <- dim(object$n)
    vname <- attr(terms(object), "term.labels")  #variable names
    # Put the elements to be printed onto a list
    pname <- (tname[3:6])[c(n, event, pyears, expected)]
    plist <- object[pname]

    if (rate) {
        pname <- c(pname, "rate")
        plist$r <- scale* object$event/object$pyears
    }
    if (ci.r) {
        pname <- c(pname, "ci.r")
        plist$ci.r <- cipoisson(object$event, object$pyears, p=conf.level) *scale
    }
    if (rr) {
        pname <- c(pname, "rr")
        plist$rr <- object$event/object$expected
    }
    if (ci.rr) {
        pname <- c(pname, "ci.rr")
        plist$ci.rr <-  cipoisson(object$event, object$expected, p=conf.level)
    }

    rname <- c(n = "N", event="Events",
               pyears= "Time", expected= "Expected events",
               rate = "Event rate", ci.r = "CI (rate)",
               rr= "Obs/Exp",   ci.rr= "CI (O/E)")
    rname <- rname[pname]           
    if (length(dd) ==1) {
        # 1 dimensional table
        cname <- names(object$n)  #category names

        if (vertical) {
            # The person-years objects list across the top, categories up and down
            # This makes columns line up in a standard "R" way
            # The first column label is the variable name, content is the categories
            plist <- lapply(plist, pformat, nastring, ...) # make it character
            pcol  <- sapply(plist, function(x) nchar(x[1])) #width of each one
            colwidth <- pmax(pcol, nchar(rname)) +2
            for (i in 1:length(plist)) 
                plist[[i]] <- strpad(plist[[i]], colwidth[i])

            colwidth <- c(max(nchar(vname), nchar(cname)) +2, colwidth)
            leftcol <- list(strpad(cname, colwidth[1]))
            header  <- strpad(c(vname, rname), colwidth)
        }
        else {
            # in this case each column will have different types of objects in it
            #  alignment is the nuisance
            newmat <- pybox(plist, length(plist[[1]]), nastring, ...)
            colwidth <- pmax(nchar(cname), apply(nchar(newmat), 1, max)) +2
            # turn the list sideways
            plist <- split(newmat, row(newmat))
            for (i in 1:length(plist))
                plist[[i]] <- strpad(plist[[i]], colwidth[i])

            colwidth <- c(max(nchar(vname), nchar(rname)) +2, colwidth)
            leftcol <- list(strpad(rname, colwidth[1]))
            header  <- strpad(c(vname, cname), colwidth)
         }

        # Now print it
        if (vline) { # use a pipe table
            cat(paste(header, collapse = "|"), "\n")
            cat(paste(strpad("-", colwidth, "-"), collapse="|"), "\n")

            temp <- do.call("paste", c(leftcol, plist, list(sep ="|")))
            cat(temp, sep= '\n')
        }                      
        else {
            cat(paste(header, collapse = " "), "\n")
            cat(paste(strpad("-", colwidth, "-"), collapse=" "), "\n")
            temp <- do.call("paste", c(leftcol, plist, list(sep =" ")))
            cat(temp, sep='\n')
        }
    } else {
        # more than 1 dimension
        if (header) {
            # the header is itself a table
            width <- max(nchar(rname))
            if (vline) {
                cat('+', strpad('-', width, '-'), "+\n", sep="")
                cat(paste0('|',strpad(rname, width), '|'), sep='\n')
                cat('+', strpad('-', width, '-'), "+\n\n", sep="")
            } else {
                cat(strpad('-', width, '-'), "\n")
                cat(strpad(rname, width), sep='\n')
                cat(strpad('-', width, '-'), "\n\n")
            }
        }
        tname <- vname[1:2]  #names for the row and col
        rowname  <- dimnames(object$n)[[1]]
        colname  <- dimnames(object$n)[[2]]
        if (length(dd) > 2) 
            newmat <- pybox(plist, c(dd[1],dd[2], prod(dd[-(1:2)])), 
                            nastring, ...)
        else  newmat <- pybox(plist, dd,  nastring, ...)

        if (length(dd) > 2) {
            newmat <- pybox(plist, c(dd[1],dd[2], prod(dd[-(1:2)])), 
                            nastring, ...)
            outer.label <- do.call("expand.grid", dimnames(object$n)[-(1:2)])
            temp <- names(outer.label)
            for (i in 1:nrow(outer.label)) {
                # first the caption, then data
                cat(paste(":", paste(temp, outer.label[i,], sep="=")), '\n')
                pyshow(newmat[,,i,], tname, rowname, colname, vline)
            }
        }
        else {
            newmat <- pybox(plist, dd,  nastring, ...)
            pyshow(newmat, tname, rowname, colname, vline)
        }
    }
    invisible(object)
}

strpad <- function(x, width, pad=' ') {
    # x = the string(s) to be padded out
    # width = width of desired string. 
    nc <- nchar(x)
    added <- width - nc

    left  <- pmax(0, floor(added/2))       # can't add negative space
    right <- pmax(0, width - (nc + left))  # right will be >= left

    if (all(right <=0)) {
        if (length(x) >= length(width)) x  # nothing needs to be done
        else rep(x, length.out=length(width))
    }
    else {
        # Each pad could be a different length.
        # Make a long string from which we can take a portion
        longpad <- paste(rep(pad, max(right)), collapse='') 
        paste0(substring(longpad, 1, left), x, substring(longpad,1, right))
    }
}

pformat <- function(x, nastring, ...) {
    # This is only called for single index tables, in vertical mode
    # Any matrix will be a confidence interval
    if (is.matrix(x)) 
        ret <- paste(ifelse(is.na(x[,1]), nastring,
                            format(x[,1],  ...)), "-", 
                     ifelse(is.na(x[,2]), nastring, 
                            format(x[,2],  ...)))
    else ret <- ifelse(is.na(x), nastring, format(x,  ...))
}
pybox <- function(plist, dd, nastring, ...) {
    ci <- (substring(names(plist), 1,3) == "ci.")  # the CI components
    int <- sapply(plist, function(x) all(x == floor(x) | is.na(x)))
    int <- (!ci & int)
    real<- (!ci & !int)
    nc <- prod(dd)
    final <- matrix("", nrow=nc, ncol=length(ci))
    
    if (any(int)) { # integers
        if (any(sapply(plist[int], length) != nc))
            stop("programming length error, notify package author")
        temp <- unlist(plist[int])
        final[,int] <- ifelse(is.na(temp), nastring, format(temp))
    }
    if (any(real)) { # floating point
        if (any(sapply(plist[real], length) != nc))
            stop("programming length error, notify package author")
        temp <- unlist(plist[real])
        final[,real] <- ifelse(is.na(temp), nastring, 
                               format(temp,  ...))
    }
    
    if (any(ci)) {
        if (any(sapply(plist[ci], length) != nc*2))
            stop("programming length error, notify package author")
        temp <- unlist(plist[ci])    
        temp <- array(ifelse(is.na(temp), nastring,
                             format(temp,  ...)),
                      dim=c(nc, 2, sum(ci)))
        final[,ci] <- paste(temp[,1,], temp[,2,], sep='-')
    }
    array(final, dim=c(dd, length(ci)))
}
pyshow <- function(dmat, labels, rowname, colname, vline) {
    # Every column is the same width, except the first
    colwidth <- c(max(nchar(rowname), nchar(labels[1])),
                  rep(max(nchar(dmat[1,1,]), nchar(colname)), length(colname)))
    colwidth[2] <- max(colwidth[2], nchar(labels[2]))
    ncol <- length(colwidth)

    dd <- dim(dmat)  # vector of length 3, third dim is the statistics
    rline <- ceiling(dd[3]/2)  #which line to put the row label on.
    if (vline) { # use a grid table
        cat("+", paste(strpad('-', colwidth, pad='-'), collapse='+'), "+\n",
            sep='')
        temp <- rep(' ', ncol); temp[2] <- labels[2]
        cat("|", paste(strpad(temp, colwidth), collapse="|"), "|\n",
            sep='')
        cat("|", paste(strpad(c(labels[1], colname), colwidth), collapse="|"),
            "|\n", sep='')
        cat("+", paste(strpad('=', colwidth, pad='='), collapse="+"), "+\n",
            sep='')
        for (i in 1:dd[1]) {
            for (j in 1:dd[3]) { #one printout line per stat
                if (j==rline) temp <- c(rowname[i], dmat[i,,j])
                else temp <- c("", dmat[i,,j])
                cat("|", paste(strpad(temp, colwidth), collapse='|'), "|\n",
                    sep='')
            }
            cat("+", paste(strpad('-', colwidth, '-'), collapse='+'), "+\n",
                sep='')
        }
    }
    else { # use a multiline table
        cat(paste(strpad('-', colwidth, '-'), collapse='-'), "\n")
        temp <- rep(' ', ncol); temp[2] <- labels[2]
        cat(paste(strpad(temp, colwidth), collapse=" "), "\n")
        cat(paste(strpad(c(labels[1], colname), colwidth), collapse=" "),
            "\n")
        cat(paste(strpad('-', colwidth, pad='-'), collapse=" "), "\n")
        for (i in 1:dd[1]) {
            for (j in 1:dd[3]) { #one printout line per stat
                if (j==rline) temp <- c(rowname[i], dmat[i,,j])
                else temp <- c("", dmat[i,,j])
                cat(paste(strpad(temp, colwidth), collapse=' '), "\n")
            }
            if (i< dd[1]) cat(" \n") #blank line
        }
        cat(paste(strpad('-', colwidth, '-'), collapse='-'), "\n")
    }
}
pytot <- function(x, na=FALSE) {
    dd <- dim(x)
    if (length(dd) ==1) {
        if (na) array(c(x, NA), dim= length(x) +1,
                              dimnames=list(c(dimnames(x)[[1]], "Total")))
        else array(c(x, sum(x)), dim= length(x) +1,
                              dimnames=list(c(dimnames(x)[[1]], "Total")))
    }
    else if (length(dd) ==2) {
        if (na) new <- rbind(cbind(x, NA), NA)
        else {
            new <- rbind(x, colSums(x))
            new <- cbind(new, rowSums(new))
            }
        array(new, dim=dim(x) + c(1,1), 
              dimnames=list(c(dimnames(x)[[1]], "Total"),
                            c(dimnames(x)[[2]], "Total")))
    }
    else {
        # The general case
        index <- 1:length(dd)
        if (na) sum1 <- sum2 <- sum3 <- NA
        else {
            sum1 <- apply(x, index[-1], sum)    # row sums
            sum2 <- apply(x, index[-2], sum)    # col sums
            sum3 <- apply(x, index[-(1:2)], sum) # total sums
            }
        
        # create a new matrix and then fill it in
        d2 <- dd
        d2[1:2] <- dd[1:2] +1
        dname <- dimnames(x)
        dname[[1]] <- c(dname[[1]], "Total")
        dname[[2]] <- c(dname[[2]], "Total")
        new <- array(x[1], dim=d2, dimnames=dname)

        # say dim(x) =(5,8,4); we want new[6,-9,] <- sum1; new[-6,9,] <- sum2
        #  and new[6,9,] <- sum3
        # if dim is longer, we need to add more commas
        commas <- rep(',', length(dd) -2)
        eval(parse(text=paste("new[1:dd[1], 1:dd[2]", commas, "] <- x")))
        eval(parse(text=paste("new[ d2[1],-d2[2]", commas, "] <- sum1")))
        eval(parse(text=paste("new[-d2[1], d2[2]", commas, "] <- sum2")))
        eval(parse(text=paste("new[ d2[1], d2[2]", commas, "] <- sum3")))
        new
    }
}
