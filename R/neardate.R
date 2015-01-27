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
neardate <- function(id1, id2, y1, y2, best=c("after", "prior"),
                     nomatch=NA_integer_) {
    if (missing(id1)) stop("id1 argument is required")
    if (missing(id2)) stop("id2 argument is required")
    if (missing(y1))  stop("y1 argument is required")
    if (missing(y2))  stop("y2 argument is required")
    if (length(id1) != length(y1))
            stop("id1 and y1 have different lengths")
    if (length(id2) != length(y2))
            stop("id2 and y2 have different lengths")

    best <- match.arg(best)
        
    # This check could be more sophisticated (though I don't see how to do it)
    #  We want to make sure that the "alldate" line below makes sense for the 
    #  data types that the user passed in.
    if (is.factor(y1) || is.factor(y2))
        stop("y1 and y2 must be sortable")
    if (inherits(y1, 'POSIXt')) 
        if (!inherits(y2, 'POSIXt')) y2 <- as(y2, class(y1))
    else if (inherits(y2, 'POSIXt')) y1 <- as(y1, class(y2))

    alldate <- sort(unique(c(y1, y2)))
    y1 <- match(y1, alldate)
    y2 <- match(y2, alldate)

    # Throw out lines with missing y2, but remember which ones
    rowid <- 1:length(y2)
    if (any(is.na(y2))) {
        toss <- is.na(y2)
        y2 <- y2[!toss]
        if (!missing(id2)) id2 <- id2[!toss]
        rowid <- rowid[!toss]
    }
    n2 <- length(y2)
    if (n2 ==0) stop("No valid entries in data set 2")

    #   Toss out any rows in id2 that are not possible targets for id1
    #   (id2 is usually the larger data set, thinning speeds it up)
    indx1 <- match(id2, id1)
    toss <- is.na(indx1) 
    if (any(toss)) {
        id2 <- id2[!toss]
        y2 <- y2[!toss]
        indx1 <- indx1[!toss]
        rowid <- rowid[!toss]
    }
    
    n2 <- length(y2)
    if (n2 ==0) stop("No valid entries in data set 2")

    # We need to create a merging id.  A minimal amount of
    #  spread for the dates keeps numeric overflow at bay
    delta <- 1.0 + length(alldate)  #numeric, not integer, on purpose
    hash1 <- match(id1, id1)*delta + y1
    hash2 <- indx1*delta + y2 

    if (best=="prior")
        indx2 <- approx(hash2, 1:n2, hash1, method="constant", yleft=NA,
                        yright=n2, rule=2, f=0)$y
    else 
        indx2 <- approx(hash2, 1:n2, hash1, method="constant", yleft=1,
                        yright=NA, rule=2, f=1)$y

    rowid[ifelse(id1== id2[indx2], indx2, NA)]
}
