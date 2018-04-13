# Methods for actual rate tables, like
#  the table of US survival rates by age and sex (survexp.us).  Rate tables
#  have the "ratetable" class.  However, since each one is rather unique, 
#  there is no function to create a rate table.  Each consists of a multi-way
#  array of event rates along with a set of attributes. 
#
"[.ratetable" <- function(x, ..., drop=TRUE) {
    aa <- attributes(x)
    attributes(x) <- aa[c("dim", "dimnames")]
    y <- NextMethod("[", drop=FALSE)
    newdim <- attr(y, 'dim')
    if (is.null(newdim)) return(y)  #when the subscript was a single vector
    dropped <- (newdim==1)
    if (drop)  change <- (newdim!=aa$dim & !dropped)
    else       change <- (newdim!=aa$dim)

    if (any(change)) {  #dims that got smaller, but not dropped
	newcut <- aa$cutpoints
	for (i in (1:length(change))[change])
	    if (!is.null(newcut[[i]])) newcut[[i]] <-
		(newcut[[i]])[match(dimnames(y)[[i]], aa$dimnames[[i]])]
	aa$cutpoints <- newcut
	}
    if (drop && any(dropped)){
	if (all(dropped)) as.numeric(y)   #single element
	else {
	    #Note that we have to drop the summary function
	    attributes(y) <- list( dim = dim(y)[!dropped],
				   dimnames = dimnames(y)[!dropped],
				   dimid = aa$dimid[!dropped],
				   factor = aa$factor[!dropped],
				   cutpoints =aa$cutpoints[!dropped],
                                   type = aa$type[!dropped])
	    class(y) <- 'ratetable'
            y
	    }
	}
    else {
	aa$dim <- aa$dimnames <- NULL
	attributes(y) <- c(attributes(y), aa)
	y
	}
    }

is.na.ratetable  <- function(x)
    structure(is.na(as.vector(x)), dim=dim(x), dimnames=dimnames(x))

Math.ratetable <- function(x, ...) {
    attributes(x) <- attributes(x)[c("dim", "dimnames")]
    NextMethod(.Generic)
    }

Ops.ratetable <- function(e1, e2) {
    #just treat it as an array
    if (nchar(.Method[1])) attributes(e1) <- attributes(e1)[c("dim","dimnames")]
    if (nchar(.Method[2])) attributes(e2) <- attributes(e2)[c("dim","dimnames")]
    NextMethod(.Generic)
    }

as.matrix.ratetable <- function(x, ...) {
    attributes(x) <- attributes(x)[c("dim", "dimnames")]
    x
    }
