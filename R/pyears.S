# $Id: pyears.S 11430 2010-07-29 22:34:23Z therneau $
pyears <- function(formula, data,
	weights, subset, na.action,
	ratetable=survexp.us, scale=365.25,  expect=c('event', 'pyears'),
	model=FALSE, x=FALSE, y=FALSE, data.frame=FALSE) {

    expect <- match.arg(expect)
    call <- match.call()
    m <- match.call(expand.dots=FALSE)
    m <- m[c(1, match(c('formula', 'data', 'weights', 'subset', 'na.action'),
		      names(m), nomatch=0))]
    
    Terms <- if(missing(data)) terms(formula, 'ratetable')
	     else              terms(formula, 'ratetable',data=data)

    rate <- attr(Terms, "specials")$ratetable
    if (length(rate) >1 )
	stop ("Can have only 1 ratetable() call in a formula")
    else if (length(rate) == 0 && !missing(ratetable)) {
	# add a 'ratetable' call to the internal formula
        # The dummy function stops an annoying warning message "Looking for
        #  'formula' of mode function, ignored one of mode ..."
	xx <- function(x) formula(x)
    
	if(is.ratetable(ratetable))   varlist <- attr(ratetable, "dimid")
	else stop("Invalid rate table")

	ftemp <- deparse(substitute(formula))
	formula <- xx( paste( ftemp, "+ ratetable(",
			  paste( varlist, "=", varlist, collapse = ","), ")"))
	Terms <- if (missing(data)) terms(formula, "ratetable")
	         else               terms(formula, "ratetable", data = data)
	rate <- attr(Terms, "specials")$ratetable
	}

    if (any(attr(Terms, 'order') >1))
	    stop("Pyears cannot have interaction terms")
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())

    Y <- model.extract(m, 'response')
    if (is.null(Y)) stop ("Follow-up time must appear in the formula")
    if (!is.Surv(Y)){
	if (any(Y <0)) stop ("Negative follow up time")
	Y <- as.matrix(Y)
	if (ncol(Y) >2) stop("Y has too many columns")
	if (ncol(Y)==2 && any(Y[,2] <= Y[,1]))
	    stop("stop time must be > start time")
	}
    else {
        stype <- attr(Y, 'type')
        if (stype == 'right') {
            if (any(Y[,1] <0)) stop("Negative survival time")
            nzero <- sum(Y[,1]==0 & Y[,2] ==1)
            if (nzero >0) 
                warning(paste(nzero, 
                         "observations with an event and 0 follow-up time,",
                       "any rate calculations are statistically questionable"))
            }
        else if (stype != 'counting')
            stop("Only right-censored and counting process survival types are supported")
        }
    
    n <- nrow(Y)
    if (is.null(n) || n==0) stop("Data set has 0 observations")

    weights <- model.extract(m, 'weights')
    if (is.null(weights)) weights <- rep(1.0, n)

    if (length(rate)==1) {
	ovars <- (dimnames(attr(Terms, 'factors'))[[1]])[-c(1, rate)]
	rtemp <- match.ratetable(m[,rate], ratetable)
	R <- rtemp$R
	if (!is.null(rtemp$call)) {#need to drop some dimensions from ratetable
	    ratetable <- eval(parse(text=rtemp$call))
	    }
	}
    else ovars <- (dimnames(attr(Terms, 'factors'))[[1]])[-1]

    # Now process the other (non-ratetable) variables
    if (length(ovars)==0)  {
	# no categories!
	X <- rep(1,n)
	ofac <- odim <- odims <- ocut <- 1
	}
    else {
	odim <- length(ovars)
	ocut <- NULL
	odims <- ofac <- double(odim)
	X <- matrix(0, n, odim)
	outdname <- vector("list", odim)
	for (i in 1:odim) {
	    temp <- m[[ovars[i]]]
	    ctemp <- oldClass(temp)
	    if (!is.null(ctemp) && ctemp=='tcut') {
		X[,i] <- temp
		temp2 <- attr(temp, 'cutpoints')
		odims[i] <- length(temp2) -1
		ocut <- c(ocut, temp2)
		ofac[i] <- 0
		outdname[[i]] <- attr(temp, 'labels')
		}
	    else {
		temp2 <- as.factor(temp)
		X[,i] <- temp2
		temp3 <- levels(temp2)
		odims[i] <- length(temp3)
		ofac[i] <- 1
		outdname[[i]] <- temp3
		}
	    }
	}

    # Now do the computations
    ocut <-c(ocut,0)   #just in case it were of length 0
    osize <- prod(odims)
    if (length(rate)) {  #include expected
	atts <- attributes(ratetable)
	cuts <- atts$cutpoints
        if (is.null(atts$type)) {
            #old stlye table
            rfac <- atts$factor
            us.special <- (rfac >1)
            }
        else {
            rfac <- 1*(atts$type ==1)
            us.special <- (atts$type==4)
            }
	if (any(us.special)) {  #special handling for US pop tables
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
            #  year's rates, etc. 
            # The birth date is entry date - age in days (based on 1/1/1960).
	    # I don't much care which date functions I use to do the arithmetic
            #  below.  Unfortunately R and Splus don't share one.  My "date"
            #  class is simple, but is also one of the earlier date class
            #  attempts, has less features than others, and will one day fade
            #  away; so I don't want to depend on it alone.
            #
	    cols <- match(c("age", "year"), atts$dimid)
      	    if (any(is.na(cols))) 
                 stop("Ratetable does not have expected shape")
            if (exists("as.Date")) {  # true for modern version of R
                bdate <- as.Date('1960/1/1') + (R[,cols[2]] - R[,cols[1]])
                byear <- format(bdate, "%Y")
                offset <- bdate - as.Date(paste(byear, "01/01", sep='/'), 
                                          origin="1960/01/01")
                }
            else if (exists('month.day.year')) { # Splus, usually
                bdate <- R[,cols[2]] - R[,cols[1]]
                byear <- month.day.year(bdate)$year
                offset <- bdate - julian(1,1,byear)
                }
            else if (exists('date.mdy')) { # Therneau's date class is available
                bdate <- as.date(R[,cols[2]] - R[,cols[1]])
                byear <- date.mdy(bdate)$year
                offset <- bdate - mdy.date(1,1,byear)
                }
            else stop("Can't find an appropriate date class\n") 
            R[,cols[2]] <- R[,cols[2]] - offset

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
	temp <- .C("pyears1",
			as.integer(n),
			as.integer(ncol(Y)),
			as.integer(is.Surv(Y)),
			as.double(Y),
		        as.double(weights),
			as.integer(length(atts$dim)),
			as.integer(rfac),
			as.integer(atts$dim),
			as.double(unlist(cuts)),
			as.double(ratetable),
			as.double(R),
			as.integer(odim),
			as.integer(ofac),
			as.integer(odims),
			as.double(ocut),
			as.integer(expect=='event'),
			as.double(X),
			pyears=double(osize),
			pn    =double(osize),
			pcount=double(if(is.Surv(Y)) osize else 1),
			pexpect=double(osize),
			offtable=double(1))[18:22]
	}
    else {
	temp <- .C('pyears2',
			as.integer(n),
			as.integer(ncol(Y)),
			as.integer(is.Surv(Y)),
			as.double(Y),
		        as.double(weights),
			as.integer(odim),
			as.integer(ofac),
			as.integer(odims),
			as.double(ocut),
			as.double(X),
			pyears=double(osize),
			pn    =double(osize),
			pcount=double(if(is.Surv(Y)) osize else 1),
			offtable=double(1)) [11:14]
	}

    if (data.frame) {
        # Create a data frame as the output, rather than a set of
        #  rate tables
        keep <- (temp$pyears >0)  # what rows to keep in the output
        names(outdname) <- ovars
        if (length(outdname) ==1) {
            # if there is only one variable, the call to "do.call" loses
            #  the variable name, since expand.grid returns a factor
            df <- data.frame((outdname[[1]])[keep], 
                             pyears= temp$pyears[keep]/scale,
                             n = temp$pn[keep])
            names(df) <- c(names(outdname), 'pyears', 'n')
            }
        else {
            df <- cbind(do.call("expand.grid", outdname)[keep,],
                             pyears= temp$pyears[keep]/scale,
                             n = temp$pn[keep])
            }
        row.names(df) <- 1:nrow(df)
        if (length(rate)) df$expected <- temp$pexpect[keep]
        if (expect=='pyears') df$expected <- df$expected/scale
        if (is.Surv(Y)) df$event <- temp$pcount[keep]

        out <- list(call=call,
                    data= df, offtable=temp$offtable/scale)  
        if (length(rate) && !is.null(rtemp$summ))
            out$summary <- rtemp$summ
        }

    else if (prod(odims) ==1) {  #don't make it an array
	out <- list(call=call, pyears=temp$pyears/scale, n=temp$pn,
		    offtable=temp$offtable/scale)
	if (length(rate)) {
	    out$expected <- temp$pexpect
	    if (expect=='pyears') out$expected <- out$expected/scale
	    if (!is.null(rtemp$summ)) out$summary <- rtemp$summ
	    }
	if (is.Surv(Y)) out$event <- temp$pcount
	}
    else {
	out <- list(call = call,
		pyears= array(temp$pyears/scale, dim=odims, dimnames=outdname),
		n     = array(temp$pn,     dim=odims, dimnames=outdname),
		offtable = temp$offtable/scale)
	if (length(rate)) {
	    out$expected <- array(temp$pexpect, dim=odims, dimnames=outdname)
	    if (expect=='pyears') out$expected <- out$expected/scale
	    if (!is.null(rtemp$summ)) out$summary <- rtemp$summ
	    }
	if (is.Surv(Y))
		out$event <- array(temp$pcount, dim=odims, dimnames=outdname)
	}
    na.action <- attr(m, "na.action")
    if (length(na.action))  out$na.action <- na.action
    if (model) out$model <- m
    else {
	if (x) out$x <- cbind(X, R)
	if (y) out$y <- Y
	}
    oldClass(out) <- 'pyears'
    out
    }

