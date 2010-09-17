# $Id: match.ratetable.S 11289 2009-06-12 18:10:51Z therneau $
# Do a set of error checks on whether any categorical tetable() vars match the
#   level set of the actual ratetable.  Continuous variables are left alone.
# This is called by pyears and survexp, but not by users
#
# The categoricals are turned into integer subscripts
#
match.ratetable <- function(R, ratetable) {
    attR <- attributes(R)
    attributes(R) <- attR['dim']     #other attrs get in the way later
    if (!is.ratetable(ratetable)) stop("Invalid rate table")
    dimid <- attr(ratetable, 'dimid')

    ord <- match(attR$dimnames[[2]], dimid)

    if (any(is.na(ord)))
       stop(paste("Argument '", (attR$dimnames[[2]])[is.na(ord)],
	    "' in ratetable()",
	    " does not match the given table of event rates", sep=''))
    nd <- length(ord)
    if (nd != length(dimid))
	stop("The ratetable() call has the wrong number of arguments")
    ord[ord] <- 1:nd   #reverse the index, so "ord" can be on the right-hand
    R <- R[,ord,drop=FALSE]  #put cols in same order as the ratetable

    # Check out the dimensions of R --
    isDate <- attR[["isDate"]][ord]
    levlist <- attR[['levlist']][ord]
    dtemp <-dimnames(ratetable)
    rtype  <- attr(ratetable, 'type') # 1= class, 2=cont, 3=date, 4=US yr
    if (is.null(rtype)) { #old style ratetable, be backwards compatable
        temp <- attr(ratetable, 'factor')
        # we map 'old continuous' to 'new date'; since it might be a date
        rtype <- 1*(temp==1) + 3*(temp==0) + 4*(temp >1)  
        }

    # Now, go through the dimensions of the ratetable 1 by 1, and
    #  verify that the user's variable is compatable
    #  with the rate table's dimensions
    #
    if (any(rtype<3 & isDate)) {
        indx <- which(rtype<1 & isDate)
        stop(paste("Data has a date type variable, but the reference",
                   "ratetable is not a date for variable", dimid[indx]))
        }
    for (i in (1:nd)) {
	if (length(levlist[[i]]) >0) {  #factor or character variable
	    if (rtype[i]!=1) stop(paste("In ratetable(),", dimid[i],
				     "must be a continuous variable"))
	    temp <- charmatch(casefold(levlist[[i]]), casefold(dtemp[[i]]))
	    if (any(is.na(temp)))
		stop(paste("Levels do not match for ratetable() variable",
			    dimid[i]))
            if (any(temp==0)) 
                stop(paste("Non-unique ratetable match for variable",
                               dimid[i]))
	    R[,i] <- temp[R[,i]]
	    }

	else {   # ratetable() thinks it is a continuous variable
	    if (rtype[i]==1) {   #but it's not-- make sure it is an integer
		temp <- R[,i]
		if (any(floor(temp)!=temp) || any(temp<=0) ||
			    max(temp) > length(dtemp[[i]]))
		stop(paste("In ratetable(),",dimid[i],"is out of range"))
		}
	    }
	}

    summ <- attr(ratetable, 'summary')
    if (is.null(summ))
	 list(R= R)
    else list(R= R, summ=summ(R))
    }
