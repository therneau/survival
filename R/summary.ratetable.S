# $Id: summary.ratetable.S 11437 2010-10-28 02:21:16Z therneau $
#
# Print out information about a rate table: it's dimensions and keywords
#
summary.ratetable <- function(object, ...) {
    rtable<-object
    if (!inherits(rtable, 'ratetable')) stop("Argument is not a rate table")

    att <- attributes(rtable)
    ncat <- length(dim(rtable))
    cat (" Rate table with", ncat, "dimensions:\n")
    for (i in 1:ncat) {
        # One of 'factor' (old style table) or "type" (new style) should exist
        if (!is.null(att$factor)) {
            if (att$factor[i]==0) {
                cat("\t", att$dimid[i], " ranges from ", 
                    format(min(att$cutpoints[[i]])), " to ", 
                    format(max(att$cutpoints[[i]])), "; with ", att$dim[i],
                    " categories\n", sep='')
                }
            else if(att$factor[i]==1) {
                cat("\t", att$dimid[i], " has levels of: ",
                    paste(att$dimnames[[i]], collapse=' '), "\n", sep='')
                }
            else {
                cat("\t", att$dimid[i], " ranges from " , 
                    format(min(att$cutpoints[[i]])), " to ", 
                    format(max(att$cutpoints[[i]])), "; with ", att$dim[i],
                    " categories,\n\t\tlinearly interpolated in ",
                    att$factor[i], " steps per division\n", sep='')
                }
            }
        else {
            if (att$type[i]==1) {
                cat("\t", att$dimid[i], " has levels of: ",
                    paste(att$dimnames[[i]], collapse=' '), "\n", sep='')
                }
            else if (att$type[i]>2) { #date
                cat("\t", att$dimid[i], " ranges from " , 
                 format(as.Date(min(att$cutpoints[[i]]), origin='1960/01/01')),
                    " to ", 
                 format(as.Date(max(att$cutpoints[[i]]), origin='1960/01/01')),
                    "; with ", att$dim[i],
                    " categories\n", sep='')
                }

            else {
                cat("\t", att$dimid[i], " ranges from ", 
                    format(min(att$cutpoints[[i]])), " to ", 
                    format(max(att$cutpoints[[i]])), "; with ", att$dim[i],
                    " categories\n", sep='')
                }
            }
        }
            
    invisible(att)
    }

