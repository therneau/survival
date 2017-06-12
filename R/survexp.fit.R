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
        # Someone born in June of 1945, say, gets the 1945 US rate until their
        #  next birthday.  But the underlying logic of the code would change
        #  them to the 1946 rate on 1/1/1946, which is the cutpoint in the
        #  rate table.  We fudge by faking their birthdate back to January 1.
        # The cutpoint for year can be simple numeric (older) or a Date object,
        #  and the older style ones use days since 1/1/1960.  (Date objects in
        #  R didn't exist when rate tables were conceived.)  So watch out for
        #  that as well.
        if (is.null(atts$dimid)) dimid <- names(atts$dimnames)
        else dimid <- atts$dimid
        cols <- match(c("age", "year"), dimid)
        if (any(is.na(cols))) 
            stop("Ratetable does not have expected shape")

        if (inherits(atts$cutpoints[[which(us.special)]], "Date"))
            origin <- "1970/01/01"  else origin <- "1960/01/01"
        bdate <- as.Date(origin) + (x[,cols[2]] - x[,cols[1]])
        byear <- format(bdate, "%Y")
        offset <- bdate - as.Date(paste(byear, "01/01", sep='/'), 
                                      origin=origin)
        x[,cols[2]] <- x[,cols[2]] - offset

	# Doctor up "cutpoints" - only needed for (very) old style rate tables
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
