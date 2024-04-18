# residuals for a survfit object
residuals.survfit <- function(object, times, type= "pstate",
                              collapse=FALSE, weighted=collapse, 
                              data.frame= FALSE, ...){
    if (!inherits(object, "survfit"))
        stop("argument must be a survfit object")
    if (object$type=="interval") {
        stop("residuals for interval-censored data are not available")
        }
    survfitms <- inherits(object, "survfitms")  # if multi-state
    coxsurv   <- inherits(object, "survfitcox")
    if (coxsurv) stop("residuals method for coxph survival curve not found")
    
    if (!is.logical(collapse)) stop("collapse must be TRUE/FALSE")
    if (!is.logical(weighted)) stop("weighted must be TRUE/FALSE")
    if (!is.logical(data.frame)) stop("data.frame must be TRUE/FALSE")

    if (!is.null(object$oldstates)) 
        stop("residuals not available for a subscripted survfit object")
 
    # allow a set of alias for 'type'
    temp <- c("pstate", "cumhaz", "sojourn", "survival",
                              "chaz", "rmst", "rmts", "auc")
    type <- match.arg(casefold(type), temp)
    itemp <-  c(1,2,3,1,2,3,3,3)[match(type, temp)]
    type <- c("pstate", "cumhaz", "auc")[itemp]

    # Skip roundoff correction for the times only if the survfit call did so
    timefix <- (is.null(object$timefix) || object$timefix)
    if (missing(times)) 
        stop ("the times argument is required")
    else {
        if (!is.numeric(times)) stop("times must be a numeric vector")
        times <- sort(unique(times))
        if (timefix) times <- aeqSurv(Surv(times))[,1]
    }
    timelab <- signif(times, 4)  # used for dimnames
    if (any(duplicated(timelab)) timelab <- NULL  # give up on nice values

    if (!is.logical(collapse)) stop("collapse must be TRUE/FALSE")
 
    # We need all of time, status, id, cluster, and weight, so grab the model
    # frame
    mf <- model.frame(object)
    n <- nrow(mf)
    if (n==0) stop("data set has no non-missing observations")
    Call <- object$call
    Terms <- terms(mf)

    if (is.null(object$y)) Y <- model.response(mf)
    else Y <- object$y
    if (inherits(Y, "Surv2")) stop("Surv2 objects not supported")
    ny <- ncol(Y)

    id <- model.extract(mf, "id")
    cluster <- model.extract(mf, "cluster")
    if (is.null(cluster)) cluster <- id
    if (is.null(cluster)) collapse <- FALSE
    if (collapse && !weighted) 
        stop("invalid combination of options: collapse=TRUE and weighted=FALSE")
     
    if (is.null(object$weight)) casewt <- model.weights(mf)
    if (is.null(casewt)) {
        weighted <- FALSE  # no weights available
        casewt <- rep(1.0, n) # we still need them for C calls
    }


    if (is.null(object$x)) {
        ll <- attr(Terms, 'term.labels')
        if (length(ll) == 0) X <- factor(rep(1,n))  # ~1 on the right
        else X <- strata(mf[ll])
    } else X <- object$x
    xlev <- levels(X)

    if (collapse && !is.null(id) && any(X != X[1])) {
        # If the same id shows up in multiple curves, we just can't deal
        #  with it.
        temp <- unlist(lapply(split(id, X), unique))
        if (any(duplicated(temp)))
            stop("same id appears in multiple curves, cannot collapse")
    }

    start.time <- object$start.time
    if (is.null(start.time)) start.time <- min(c(0, object$time))

    # remember the name of the id variable, if present, to use for a label
    #  but we don't try to parse it:  id= mydata$clinic becomes NULL
    idname <- Call$id
    if (is.name(idname)) idname <- as.character(idname)
    else idname <- "(id)"  
   
    # if missing id, create one that is 'row number in original data'
    if (length(id) ==0) {
        id <- seq(n + length(object$na.action))
        if (!is.null(object$na.action)) id <- id[-object$na.action]  
    }

   # What type of curve?  # to do, support the old 'type' arg
    stype <- Call$stype
    if (is.null(stype)) stype <- 1
    ctype <- Call$ctype
    if (is.null(ctype)) ctype <- 1
    
    
    rowid <- as.integer(X) # which rows of X go with each curve
    ncurve <- max(rowid)
    ntime <- length(times)
    if (!survfitms) {
        resid <- matrix(0, n, ntime)
        for (i in 1:ncurve) {
            i1 <- which(rowid ==i)
            resid[i1,] <- rsurvpart1(Y[i1,], times, type, stype, ctype,
                                object[i])
        }
 
        if (collapse) {
            resid <- rowsum(resid*casewt, id, reorder=FALSE)
            dimnames(resid) <- list(id= unique(id), times=timelab)
            curve <- (as.integer(X))[!duplicated(id)] #which curve for each row
        } else {
            if (weighted && any(casewt !=1)) resid <- resid*casewt
            dimnames(resid) <- list(id=id, times=timelab)
            curve <- as.integer(X)
        }
    }
    else {  # multi-state
        if (length(id) ==0) { # only possible if ny=2, competing risks
            states <- c("(s0)", attr(Y, "states"))
            istate <- factor(rep("(s0)", nrow(Y)), states)
        } else {
            istate  <- model.extract(mf, "istate")
            if (is.null(istate)) mcheck <- survcheck2(Y, id)  
            else mcheck <- survcheck2(Y, id, istate)
            istate <- mcheck$istate  # Normalize the states
            states <- mcheck$states
        }
        if (!identical(states, object$states))
            stop("error in residuals.survfit, non-matching states")
        if (!identical(states, attr(Y, "states"))) { # new states are a superset
            smap <- c(0, match(attr(Y, "states"), states))
            Y[,ncol(Y)] <- smap[Y[,ncol(Y)] +1]  
            attr(Y, 'states') <- states
        }       

        nstate <- length(states)
        nhaz <- ncol(object$cumhaz)
        if (type== "cumhaz") {
            sname <- colnames(object$cumhaz)
            resid <- array(0, dim=c(n, nhaz, ntime))
        } else {
            resid <- array(0., dim=c(n, nstate, ntime))
            sname <- object$states
        } 

        if (ncurve==1) 
            resid <- rsurvpart2(Y, casewt, istate,
                                     times, type, object)
        else for (i in 1:ncurve) {
            i1 <- which(rowid ==i)
            resid[i1,,] <- rsurvpart2(Y[i1,], casewt[i1], istate[i1],
                                     times, type, object[i,])
        }

        if (collapse) {
             dd <- dim(resid)
             if (weighted) 
                 resid <- rowsum(casewt*matrix(resid, nrow=dd[1]),
                                   cluster, reorder=FALSE)
             else resid <- rowsum(matrix(resid, nrow=dd[1]), cluster, 
                                   reorder=FALSE)
             dim(resid) <- c(nrow(resid), dd[-1])
             if (length(times) >1)
                 dimnames(resid) <- list(id= unique(id), sname, times=timelab)
             else dimnames(resid) <- list(id= unique(id), sname)
             curve <- (as.integer(X))[!duplicated(id)] #which curve for each row
         } else {
             if (weighted && any(casewt != 1)) resid <- resid*casewt
             if (length(times) ==1) dimnames(resid) <- list(id=id, sname)
             else dimnames(resid) <- list(id=id, sname, times=timelab)
             curve <- as.integer(X)
        }       
    }       

    names(dimnames(resid))[1] <- idname
    if (ncurve==1) curve <- NULL
    
    # deal with na.action
    if (!is.null(object$na.action) && !collapse && !data.frame) {
        test <- seq(dim(resid)[1])
        # if naresid does nothing, i.e., the default na.omit, do nothing
        if (!identical(test, naresid(object$na.action, test))) {
            if (length(dim(resid)) > 2) {
               r2 <- naresid(object$na.action, matrix(resid, nrow=dim(resid)[1]))
               d2 <- dim(resid)[-1]
               resid <- array(r2, dim= c(length(r2)/prod(d2), d2))
            } else resid <- naresid(object$na.action, resid)
            if (length(id)) id <- naresid(object$na.action, id)
            if (length(curve)) curve <- naresid(object$na.action,curve)
        }
    }
                          
   if (!data.frame) resid
   else {
       rname <- dimnames(resid)
       rd <- dim(resid)
       if (length(rd) < 2) {
           # single time point, simple survival
           rdat <- data.frame(id=id, time=times, resid=resid)
           if (length(curve)>0) rd$curve <- curve
       } else {
           id <- rep(id, prod(rd[-1]))
           if (!survfitms) # simple surv, multiple times
                rdat <- data.frame(id=id, time=times[col(resid)],resid=c(resid))
           else { #multistate
               if (length(times) ==1)
                   rdat <- data.frame(id=id, 
                                    state= rname[[2]][col(resid)], 
                                    time= times, resid=c(resid))
               else rdat <- data.frame(id=id, 
                                  state= rep(rep(rname[[2]], each=rd[1]),rd[3]),
                                  time= rep(times, each= rd[1]*rd[2]),
                                  resid= c(resid))
               if (type=="cumhaz") names(rdat)[2] <- "transition"
           }
           if (length(curve) >0) rdat$curve <- rep(curve, prod(rd[-1]))
       }
       names(rdat)[1] <- idname        
       rdat
   }        
}


# The working code for single endpoint survival
rsurvpart1 <- function(Y, times, type, stype, ctype, fit) {
    # Y = data for a single survival curve
    # times, the desired output times (ordered from first to last)
    # type = type of residual desired
    # stype, ctype type of survival curve, and type of hazard
    # fit = the fitted survival object for this curve
    ntime <- length(times)
    events <- (fit$n.event >0)  # minor speedup by only looking at events
    hazard <- diff(c(0, fit$cumhaz[events]))
    nrisk  <- fit$n.risk[events]
    surv   <- fit$surv[events]
    dtime  <- fit$time[events]
    n <- nrow(Y)
    ny <- ncol(Y)
    status <- Y[,ncol(Y)]

    # Create the index of the reporting times into the survival curve
    #   tindex = largest event time <= reporting time
    #   yindex = largest event time <= per-subject event/censor time
    #   sindex = largest event time <= per-subject entry time
    tindex <- findInterval(times, dtime, left.open=FALSE)
    yindex <- findInterval(Y[, ny-1], dtime, left.open=FALSE)
    if (ny==3) sindex <- findInterval(Y[,1], dtime, left.open=FALSE)
    
    # A common operation is that resid[i,j] is updated using 
    #    xxx[min(yindex[i], tindex[j])] for some vector xxx
    # For the dN term the rule is 
    #   if (death and tindex >= yindex) then yindex  else 0
    # i.e. for any row dN applies to all reporting times at or
    #  after the death time for that observation
    ymin <- outer(yindex, tindex, pmin)
    dmin <- outer(ifelse(status==0L, 0L, yindex), tindex,
                  function(a, b) ifelse(a==0 | a>b, 0L, a)) 
    if (ny==3) smin <- outer(sindex, tindex, pmin)
 
    if (type=="cumhaz" || (type=="pstate" && stype ==2)) {
        # hazard is the primary thing
        hsum <- cumsum(c(0, hazard/nrisk))
        term1 <- c(0, 1/nrisk)[1+ dmin]   # the dN part
        term2 <- hsum[1+ ymin]            # the d\lambda part
        if (ny ==2) resid <- matrix(term1- term2, nrow =n)
        else {
            # events happen at the end of an interval, so no dN at the start
            term3 <- hsum[1 + smin]
            resid <- matrix(term1 + term3 - term2, nrow =n)
        }

        if (ctype==2) { 
            warning("code for ctype=2 not yet completed, result is approximate")
            # If there are d tied deaths at some time, we need to think of
            #  the increment as d separate steps, and the deriv as a sum over
            #  those steps.  
            # We can't pull that off the survfit object, so it's going to be
            #  real work using the raw data.  It is a low priority task and
            #  may never be filled in.
        }
            
        if (type=="pstate") {
            # survival is exp(-cumhaz), deriv is -S(t)* deriv(cumhaz)
            resid <- -resid * c(1, surv)[1+ tindex[col(resid)]]
        }
    } else if (type=="pstate") {
        # this is likely the most used branch
        # avoid a 0/0 issue which will arise when S(t) =0 and hazard =1
        temp <- ifelse(hazard==1, 1, 1-hazard)
        hsum <- cumsum(c(0, hazard/(nrisk* temp)))
        term1 <- c(0, 1/(temp*nrisk))[1+ dmin]  # the dN portion
        term2 <- hsum[1+ ymin]                  # the d\lambda portion
        stemp <- c(1, surv)[1+ rep(tindex, each=n)]
        if (ny==2) resid <- matrix(stemp *(term2 - term1), nrow= n)
        else {
            term3 <- hsum[1+smin]
            resid <- matrix(stemp*(term2 - (term3 + term1)), nrow=n)
        }
    } else if (type== "auc") {
        # The AUC is the area under the survival curve
        # see survfit:AUC in the methods document
        dtime <- fit$time[events]
        delta <- diff(dtime)
        j <- length(dtime) #temp index
        aucd <- cumsum(c(0,surv[-j] * delta))   #AUC from dtime[1] to dtime[k]
        if (max(times) > max(dtime)) {
            auctau <- approx(c(dtime, max(times)),
                             c(aucd, aucd[j] + surv[j]*(max(times)-dtime[j])),
                             times, yleft=0, rule=2)$y
        } else auctau <- approx(dtime, aucd, times, yleft=0)$y #dtime[1] to times
         
        if (stype==2) dd <- 1/nrisk   # the denominator for S= exp(-cumhaz)
        else          dd <- 1/(nrisk* (1-hazard)) # for KM
        dd <- ifelse(is.finite(dd), dd, 0)  # past the last death
        wtmat <- outer(auctau, -aucd, '+') # row i is AUC from dtime[i] to tau
        # Each column of resid has a different weight vector
        resid <- matrix(0, n, ntime)

        for (i in 1:ntime ){
            hsum <- cumsum(c(0, wtmat[i,]*hazard*dd))
            term1 <- c(0, wtmat[i,]*dd)[1+ dmin[,i]]   # the dN part
            term2 <- hsum[1+ ymin[,i]]                # the d\lambda part
            if (ny ==2) resid[,i] <- term2-term1
            else {
                term3 <- hsum[1 + smin[,i]]
                resid[,i] <- term2 - (term1 + term3) 
            }
        }
     } else stop("unknown type")

    resid
}
