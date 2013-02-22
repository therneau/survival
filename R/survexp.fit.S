#  Actually compute the expected survival for one or more cohorts
#    of subjects.  If each subject is his/her own group, it gives individual
#    survival
#  group = groups (one curve per group)
#  x matrix contains the rate
#    table indices = starting point for each obs in the rate table.
#  y is the number of follow-up days for each subject
#  times = the time points at which survival is desired
#  death = T if we want the conditional estimate
survexp.fit <- function(group, x, y, times, death, ratetable) {
   if (!is.matrix(x)) stop("x must be a matrix")
    if (ncol(x) != length(dim(ratetable)))
	stop("x matrix does not match the rate table")
    atts <- attributes(ratetable)
    ngrp <- max(group)
    times <- sort(unique(times))
    if (any(times <0)) stop("Negative time point requested")
    if (missing(y))  y <- rep(max(times), nrow(x))
    ntime <- length(times)
    if (!is.logical(death)) stop("Invalid value for death indicator")

    cuts <- atts$cutpoints

    if (is.null(atts$type)) {
        # old style rate table
        rfac <- atts$factor
        us.special <- (rfac >1)
        }
    else {
        rfac <- 1*(atts$type ==1)
        us.special <- (atts$type==4)
        }
    if (any(us.special)) {  #special handling for US pop tables
	if (sum(us.special) >1)
	    stop("Two columns marked for special handling as a US rate table")
	# Now, the 'entry' date on a US rate table is the number of days 
        #  since 1/1/1960, and the user data has been aligned to the
        #  same system by match.ratetable and marked as "year".
        # US rate tables are odd: the entry for age (year=1970, age=55)
        #  contains the daily rate for anyone who turns 55 in that year,
        #  from their birthday forward for 365 days.  So if your birthday
        #  is on Oct 2, the 1970 table applies from 2Oct 1970 to 1Oct 1971.
        # The underlying C code wants to make the 1970 rate table apply
        #  from 1Jan 1970 to 31Dec 1970.  The easiest way to finess this is
        #  to fudge everyone's enter-the-study date.  If you were born
        #  in March but entered in April, make it look like you entered in
        #  Febuary; that way you get the first 11 months at the entry 
        #  year's rates, etc.  This is the same as being born on Jan 1.
        # The birth date is entry date - age in days (based on 1/1/1960).
        #
        cols <- match(c("age", "year"), atts$dimid)
        if (any(is.na(cols))) 
            stop("Ratetable does not have expected shape")
        if (exists("as.Date")) {  # true for modern version of R
            bdate <- as.Date('1960/1/1') + (x[,cols[2]] - x[,cols[1]])
            byear <- format(bdate, "%Y") # year of birth
            offset <- as.numeric(bdate - 
                                 as.Date(paste(byear, '01/01', sep='/')))
            }
        # The lines below were commented out to stop spurious warning
        #   messages from "CMD check".  They are very unlikely to ever
        #   be needed, so no big loss.
        #else if (exists('month.day.year')) { # Splus, usually
        #    bdate <- x[,cols[2]] - x[,cols[1]]
        #    byear <- month.day.year(bdate)$year
        #    offset <- bdate - julian(1,1,byear)
        #    }
        #else if (exists('date.mdy')) { # the TMT date class is available
        #    bdate <- as.date(x[,cols[2]] - x[,cols[1]])
        #    byear <- date.mdy(bdate)$year
        #    offset <- bdate - mdy.date(1,1,byear)
        #    }
        else stop("Can't find an appropriate date class\n") 
        x[,cols[2]] <- x[,cols[2]] - offset

	# Doctor up "cutpoints" - only needed for old style rate tables
        #  for which the C code does interpolation on the fly
        if (any(rfac >1)) {
            temp <-  which(us.special)
            nyear <- length(cuts[[temp]])
            nint <- rfac[temp]       #intervals to interpolate over
            cuts[[temp]] <- round(approx(nint*(1:nyear), cuts[[temp]],
					nint:(nint*nyear))$y - .0001)
            }
	}

   storage.mode(x) <- storage.mode(y) <- "double"
   storage.mode(times) <- "double"
   temp <- .Call(Cpyears3b,
                 as.integer(death),
                 as.integer(rfac),
                 as.integer(atts$dim),
                 as.double(unlist(cuts)),
                 ratetable,
                 as.integer(group),
                 x, y, times,
                 as.integer(ngrp))

   if (ntime==1) list(surv=temp$surv, n=temp$n)
   else if (ngrp >1)
	 list(surv=apply(matrix(temp$surv, ntime, ngrp),2,cumprod),
		 n=   matrix(temp$n, ntime, ngrp))
    else list(surv=cumprod(temp$surv), n=temp$n)
    }
