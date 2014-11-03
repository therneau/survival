#
# Create a "nearest date" index
#  date1: the trial date
#  date2: target to match to
#
# result: an index vector for data set 1, which shows the row in data set
#  2 that has the same id, and the best date.
#
# best = "after"  The closest date in #2 that is on or after the date in #1
#        "prior"  The closest date in #2 that is on or before the date in #1
# 
neardate <- function(date1, date2, id1, id2, best=c("after", "prior"),
                     nomatch=NA_integer_) {
    if (!missing(id1)) {
        if (length(id1) != length(date1))
            stop("id1 and date1 have different lengths")
        if (missing(id2))
            stop("either both or neither of id1 and id2 must be supplied")
        if (length(id2) != length(date2))
            stop("id2 and date2 have different lengths")
    }
    else if (!missing(id2))
        stop("either both or neither of id1 and id2 must be supplied")
    best <- match.arg(best)
        
    # This check could be more sophisticated (though I don't see how to do it)
    #  We want to make sure that the "alldate" line below it make sense for the 
    #  data types that the user passed in.
    if (is.factor(date1) || is.factor(date2))
        stop("date1 and date2 must be sortable")
    if (inherits(date1, 'POSIXt')) 
        if (!inherits(date2, 'POSIXt')) date2 <- as(date2, class(date1))
    else if (inherits(date2, 'POSIXt')) date1 <- as(date1, class(date2))

    alldate <- sort(unique(c(date1, date2)))
    date1 <- match(date1, alldate)
    date2 <- match(date2, alldate)

    # Throw out missing dates in the second arg, but remember which ones
    rowid <- 1:length(date2)
    if (any(is.na(date2))) {
        toss <- is.na(date2)
        date2 <- date2[!toss]
        if (!missing(id2)) id2 <- id2[!toss]
        rowid <- rowid[!toss]
    }
    n2 <- length(date2)
    if (n2 ==0) stop("No valid entries in data set 2")

    # Simpler case: no id variables
    rowid <- 1:length(date2)
    if (missing(id1)) {
        if (best=="prior")
            indx2 <- approx(date2, 1:n2, date1, method="constant", yleft=NA,
                            yright=n2, rule=2, f=0)$y
        else 
            indx2 <- approx(date2, 1:n2, date1, method="constant", yleft=1,
                            yright=NA, rule=2, f=1)$y
        return(rowid[indx2])
    }

    # match id values as well
    #   First toss out any rows in id2 that are not possible targets for id1
    #   (id2 is usually the larger data set, thinning speeds it up)
    indx1 <- match(id2, id1)
    toss <- is.na(indx1) 
    if (any(toss)) {
        id2 <- id2[!toss]
        date2 <- date2[!toss]
        indx1 <- indx1[!toss]
        rowid <- rowid[!toss]
    }
    
    n2 <- length(date2)
    if (n2 ==0) stop("No valid entries in data set 2")

    # We need to create a merging id.  A minimal amount of
    #  spread for the dates keeps numeric overflow at bay
    delta <- 1.0 + length(alldate)  #numeric, not integer, on purpose
    hash1 <- match(id1, id1)*delta + date1
    hash2 <- indx1*delta + date2 

    if (best=="prior")
        indx2 <- approx(hash2, 1:n2, hash1, method="constant", yleft=NA,
                        yright=n2, rule=2, f=0)$y
    else 
        indx2 <- approx(hash2, 1:n2, hash1, method="constant", yleft=1,
                        yright=NA, rule=2, f=1)$y
    result <- rowid[ifelse(id1== id2[indx2], indx2, 0)]
    ifelse(result==0, nomatch, result)
}
