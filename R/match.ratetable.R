# Do a set of error checks on whether any categorical vars match the
#   level set of the actual ratetable.  If so they are mapped to the levels
#   found in the ratetable.  Dates need to match dates, and others are set
#   to simple numerics with unclass().  A matrix is returned.
# This is called by pyears and survexp, but not by users
#
# The categoricals are turned into integer subscripts
#
match.ratetable <- function(R, ratetable) {
    datecheck <- function(x) 
        inherits(x, c("Date", "POSIXt", "date", "chron")

    if (!is.ratetable(ratetable)) stop("Invalid rate table")
    dimid <- names(dimnames(ratetable))
    if (is.null(dimid)) attr(ratetable, 'dimid')  # older ratetable
    datecut <- sapply(attr(ratetable, "cutpoints"), datecheck)

    rtype  <- attr(ratetable, 'type') # 1= class, 2=cont, 3=date, 4=US yr
    if (is.null(rtype)) { #old style ratetable, be backwards compatable
        temp <- attr(ratetable, 'factor')
        rtype <- 1*(temp==1) + ifelse(datecut, 3,2)*(temp==0) + 4*(temp >1)
    }
    # is.ratetable has ensured that rtype agrees with datecut

    if (is.matrix(R)) {  
        # not possible any more, but just in case
        stop("internal error in match.ratetable, contact author")
    }

    Rnames <- names(R)
    ord <- match(dimid, Rnames)
    # This should have already been checked in pyears or survexp
    if (any(is.na(ord)))
       stop(paste("Argument '", dimid[is.na(ord)],
	    "' needed by the ratetable was not found in the data", sep=''))
    # Neither should this -- two argments matched one of the dimids -- since
    #  I demand an exact match
    if (any(duplicated(ord)))
        stop("A ratetable argument appears twice in the data")
    R <- R[,ord,drop=FALSE]  #put cols in same order as the ratetable
    levlist<- lapply(R, levels)
    isDate <- sapply(R, datecheck)

    # Now, go through the dimensions of the ratetable 1 by 1, and
    #  verify that the user's variable is compatable
    #  with the rate table's dimensions
    #
    dtemp <-dimnames(ratetable)
    if (any((rtype<3) & isDate)) {
        indx <- which(rtype<3 & isDate)
        stop(paste("Data has a date type variable, but the reference",
                   "ratetable is not a date variable", dimid[indx]))
        }
    if (any((rtype>2) & !isDate)) {
        indx <- which(rtype>2 & !isDate)
        stop(paste("the reference ratetable expects a date for variable",
                    dimid[indx]))
        }
    for (i in (1:ncol(R))) {
	if (length(levlist[[i]]) >0) {  #factor or character variable
	    if (rtype[i]!=1) stop(paste("for this ratetable,", dimid[i],
				     "must be a continuous variable"))
	    temp <- charmatch(casefold(levlist[[i]]), casefold(dtemp[[i]]))
	    if (any(is.na(temp)))
		stop(paste("Levels do not match for ratetable() variable",
			    dimid[i]))
            if (any(temp==0)) 
                stop(paste("Non-unique ratetable match for variable",
                               dimid[i]))
	    R[,i] <- temp[as.numeric(R[,i])]
	    }

	else {   # user's data isn't a factor or date
            R[,i] <- unclass(R[,i])  # get rid of difftimes & other such
	    if (rtype[i]==1) {   #ratetable is a factor: ok if data is integer
		temp <- R[,i]
		if (any(floor(temp)!=temp) || any(temp<=0) ||
			    max(temp) > length(dtemp[[i]]))
		stop(paste("The variable", dimid[i], "is out of range"))
		}
	    }
	}
    R <- as.matrix(R)

    summ <- attr(ratetable, 'summary')
    if (is.null(summ))
	 list(R= R)
    else list(R= R, summ=summ(R))
    }
