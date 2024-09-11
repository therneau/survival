# Automatically generated from the noweb directory
pyears <- function(formula, data,
        weights, subset, na.action, rmap,
        ratetable, scale=365.25,  expect=c('event', 'pyears'),
        model=FALSE, x=FALSE, y=FALSE, data.frame=FALSE) {

    expect <- match.arg(expect)
    Call <- match.call()
        
    # create a call to model.frame() that contains the formula (required)
    #  and any other of the relevant optional arguments
    # then evaluate it in the proper frame
    indx <- match(c("formula", "data", "weights", "subset", "na.action"),
                      names(Call), nomatch=0) 
    if (indx[1] ==0) stop("A formula argument is required")
    tform <- Call[c(1,indx)]  # only keep the arguments we wanted
    tform[[1L]] <- quote(stats::model.frame)  # change the function called

    Terms <- if(missing(data)) terms(formula)
             else              terms(formula, data=data)
    if (any(attr(Terms, 'order') >1))
            stop("Pyears cannot have interaction terms")

    if (!missing(rmap) || !missing(ratetable)) {
        has.ratetable <- TRUE
        if (missing(ratetable)) stop("No rate table specified")
        if (!missing(rmap)) {
            rcall <- substitute(rmap)
            if (!is.call(rcall) || rcall[[1]] != as.name('list'))
                stop("Invalid rcall argument")
            }
        else rcall <- NULL   # A ratetable, but no rcall argument

        # Check that there are no illegal names in rcall, then expand it
        #  to include all the names in the ratetable
        if (is.ratetable(ratetable))   {
            varlist <- names(dimnames(ratetable))
            if (is.null(varlist)) varlist <- attr(ratetable, "dimid") # older style
        }
        else if(inherits(ratetable, "coxph") && !inherits(ratetable, "coxphms")) {
            ## Remove "log" and such things, to get just the list of
            #   variable names
            varlist <- all.vars(delete.response(ratetable$terms))
            }
        else stop("Invalid rate table")

        temp <- match(names(rcall)[-1], varlist) # 2,3,... are the argument names
        if (any(is.na(temp)))
            stop("Variable not found in the ratetable:", (names(rcall))[is.na(temp)])
            
        if (any(!(varlist %in% names(rcall)))) {
            to.add <- varlist[!(varlist %in% names(rcall))]
            temp1 <- paste(text=paste(to.add, to.add, sep='='), collapse=',')
            if (is.null(rcall)) rcall <- parse(text=paste("list(", temp1, ")"))[[1]]
            else {
                temp2 <- deparse(rcall)
                rcall <- parse(text=paste("c(", temp2, ",list(", temp1, "))"))[[1]]
                }
            }
        # Create a temporary formula, used only in the call to model.frame
        newvar <- all.vars(rcall)
        if (length(newvar) > 0) {
            temp <- paste(paste(deparse(Terms), collapse=""),  
                           paste(newvar, collapse='+'), sep='+')
            tform$formula <- as.formula(temp, environment(Terms))
            }
        }
    else has.ratetable <- FALSE

    mf <- eval(tform, parent.frame())

    Y <- model.extract(mf, 'response')
    if (is.null(Y)) stop("Follow-up time must appear in the formula")
    if (!is.Surv(Y)){
        if (any(Y <0)) stop("Negative follow up time")
        Y <- as.matrix(Y)
        if (ncol(Y) >2) stop("Y has too many columns")
        }
    else {
        stype <- attr(Y, 'type')
        if (stype == 'right') {
            if (any(Y[,1] <0)) stop("Negative survival time")
            nzero <- sum(Y[,1]==0 & Y[,2] ==1)
            if (nzero >0) 
                warning(sprintf(ngettext(nzero,
                         "%d observation with an event and 0 follow-up time, any rate calculations are statistically questionable",
                         "%d observations with an event and 0 follow-up time, any rate calculations are statistically questionable",
                          domain = "R-survival"), nzero), domain = NA)
            }
        else if (stype != 'counting')
            stop("Only right-censored and counting process survival types are supported")
        }

    n <- nrow(Y)
    if (is.null(n) || n==0) stop("Data set has 0 observations")

    weights <- model.extract(mf, 'weights')
    if (is.null(weights)) weights <- rep(1.0, n)
    # rdata contains the variables matching the ratetable
    if (has.ratetable) {
        rdata <- data.frame(eval(rcall, mf), stringsAsFactors=TRUE)  
        if (is.ratetable(ratetable)) {
            israte <- TRUE
            rtemp <- match.ratetable(rdata, ratetable)
            R <- rtemp$R
            }
        else if (inherits(ratetable, 'coxph') && !inherits(ratetable, "coxphms")) {
            israte <- FALSE
            Terms <- ratetable$terms
            if (!is.null(attr(Terms, 'offset')))
                stop("Cannot deal with models that contain an offset")
            strats <- attr(Terms, "specials")$strata
            if (length(strats))
                stop("pyears cannot handle stratified Cox models")

            R <- model.matrix.coxph(ratetable, data=rdata)
            }
        else stop("Invalid ratetable")
        }
    ovars <- attr(Terms, 'term.labels')
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
        names(outdname) <- attr(Terms, 'term.labels')
        for (i in 1:odim) {
            temp <- mf[[ovars[i]]]
            if (inherits(temp, 'tcut')) {
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
    ocut <-c(ocut,0)   #just in case it were of length 0
    osize <- prod(odims)
    if (has.ratetable) {  #include expected
        atts <- attributes(ratetable)
        datecheck <- function(x) 
            inherits(x, c("Date", "POSIXt", "date", "chron"))
        cuts <- lapply(attr(ratetable, "cutpoints"), function(x)
            if (!is.null(x) & datecheck(x)) ratetableDate(x) else x)

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
            if (sum(us.special) > 1) stop("more than one type=4 in a rate table")
            # Someone born in June of 1945, say, gets the 1945 US rate until their
            #  next birthday.  But the underlying logic of the code would change
            #  them to the 1946 rate on 1/1/1946, which is the cutpoint in the
            #  rate table.  We fudge by faking their enrollment date back to their
            #  birth date.
            #
            # The cutpoint for year has been converted to days since 1/1/1970 by
            #  the ratetableDate function.  (Date objects in R didn't exist when 
            #  rate tables were conceived.) 
            if (is.null(atts$dimid)) dimid <- names(atts$dimnames)
            else dimid <- atts$dimid
            cols <- match(c("age", "year"), dimid)
            if (any(is.na(cols))) 
                stop("ratetable does not have expected shape")

            # The format command works for Dates, use it to get an offset
            bdate <- as.Date("1970-01-01") + (R[,cols[2]] - R[,cols[1]])
            byear <- format(bdate, "%Y")
            offset <- as.numeric(bdate - as.Date(paste0(byear, "-01-01")))
            R[,cols[2]] <- R[,cols[2]] - offset
       
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
        docount <- is.Surv(Y)
        temp <- .C(Cpyears1,
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
                        pcount=double(if(docount) osize else 1),
                        pexpect=double(osize),
                        offtable=double(1))[18:22]
        }
    else {   #no expected
        docount <- as.integer(ncol(Y) >1)
        temp <- .C(Cpyears2,
                        as.integer(n),
                        as.integer(ncol(Y)),
                        as.integer(docount),
                        as.double(Y),
                        as.double(weights),
                        as.integer(odim),
                        as.integer(ofac),
                        as.integer(odims),
                        as.double(ocut),
                        as.double(X),
                        pyears=double(osize),
                        pn    =double(osize),
                        pcount=double(if (docount) osize else 1),
                        offtable=double(1)) [11:14]
        }
    has.tcut <- any(sapply(mf, function(x) inherits(x, 'tcut')))
    if (data.frame) {
        # Create a data frame as the output, rather than a set of
        #  rate tables
        if (length(ovars) ==0) {  # no variables on the right hand side
            keep <- TRUE
            df <- data.frame(pyears= temp$pyears/scale,
                             n = temp$n)
        }
        else {
            keep <- (temp$pyears >0)  # what rows to keep in the output
            # grab prototype rows from the model frame, this preserves class
            #  (unless it is a tcut variable, then we know what to do)
            tdata <- lapply(1:length(ovars), function(i) {
                temp <- mf[[ovars[i]]]
                if (inherits(temp, "tcut")) { #if levels are numeric, return numeric
                    if (is.numeric(outdname[[i]])) outdname[[i]]
                    else  factor(outdname[[i]], outdname[[i]]) # else factor
                }
                else temp[match(outdname[[i]], temp)]
            })
            tdata$stringsAsFactors <- FALSE  # argument for expand.grid
            df <- do.call("expand.grid", tdata)[keep,,drop=FALSE]
            names(df) <- ovars
            df$pyears <- temp$pyears[keep]/scale
            df$n <- temp$pn[keep]
        }
        row.names(df) <- NULL   # toss useless 'creation history'
        if (has.ratetable) df$expected <- temp$pexpect[keep]
        if (expect=='pyears') df$expected <- df$expected/scale
        if (docount) df$event <- temp$pcount[keep]
        # if any of the predictors were factors, make them factors in the output
        for (i in 1:length(ovars)){
            if (is.factor( mf[[ovars[i]]]))
                df[[ovars[i]]] <- factor(df[[ovars[i]]], levels( mf[[ovars[i]]]))
        }

        out <- list(call=Call,
                    data= df, offtable=temp$offtable/scale,
                    tcut=has.tcut)
        if (has.ratetable && !is.null(rtemp$summ))
            out$summary <- rtemp$summ
    }

    else if (prod(odims) ==1) {  #don't make it an array
        out <- list(call=Call, pyears=temp$pyears/scale, n=temp$pn,
                    offtable=temp$offtable/scale, tcut = has.tcut)
        if (has.ratetable) {
            out$expected <- temp$pexpect
            if (expect=='pyears') out$expected <- out$expected/scale
            if (!is.null(rtemp$summ)) out$summary <- rtemp$summ
        }
        if (docount) out$event <- temp$pcount
    }
    else {
        out <- list(call = Call,
                pyears= array(temp$pyears/scale, dim=odims, dimnames=outdname),
                n     = array(temp$pn,     dim=odims, dimnames=outdname),
                offtable = temp$offtable/scale, tcut=has.tcut)
        if (has.ratetable) {
            out$expected <- array(temp$pexpect, dim=odims, dimnames=outdname)
            if (expect=='pyears') out$expected <- out$expected/scale
            if (!is.null(rtemp$summ)) out$summary <- rtemp$summ
        }
        if (docount)
                out$event <- array(temp$pcount, dim=odims, dimnames=outdname)
    }
    out$observations <- nrow(mf)
    out$terms <- Terms
    na.action <- attr(mf, "na.action")
    if (length(na.action))  out$na.action <- na.action
    if (model) out$model <- mf
    else {
        if (x) out$x <- X
        if (y) out$y <- Y
    }
    class(out) <- 'pyears'
    out
    }
