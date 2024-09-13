# Automatically generated from the noweb directory
finegray <- function(formula, data, weights, subset, na.action= na.pass,
                     etype, prefix="fg", count="", id, timefix=TRUE) {
    Call <- match.call()
    indx <- match(c("formula", "data", "weights", "subset", "id"),
              names(Call), nomatch=0) 
    if (indx[1] ==0) stop(gettextf("'%s' argument is required", "formula"))
    temp <- Call[c(1,indx)]  # only keep the arguments we wanted
    temp$na.action <- na.action
    temp[[1L]] <- quote(stats::model.frame)  # change the function called

    special <- c("strata", "cluster")
    temp$formula <- if(missing(data)) terms(formula, special)
    else              terms(formula, special, data=data)

    mf <- eval(temp, parent.frame())
    if (nrow(mf) ==0) stop("No (non-missing) observations")
    Terms <- terms(mf)

    Y <- model.extract(mf, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")
    type <- attr(Y, "type")
    if (type!='mright' && type!='mcounting')
        stop("Fine-Gray model requires a multi-state survival")
    nY <- ncol(Y)
    states <- attr(Y, "states")
    if (timefix) Y <- aeqSurv(Y)

    strats <- attr(Terms, "specials")$strata
    if (length(strats)) {
        stemp <- untangle.specials(Terms, 'strata', 1)
        if (length(stemp$vars)==1) strata <- mf[[stemp$vars]]
        else strata <- survival::strata(mf[,stemp$vars], shortlabel=TRUE)
        istrat <- as.numeric(strata)
        mf[stemp$vars] <- NULL
        }
    else istrat <- rep(1, nrow(mf))
    
    id <- model.extract(mf, "id")
    if (!is.null(id)) mf["(id)"] <- NULL  # don't leave it in result
    user.weights <- model.weights(mf)
    if (is.null(user.weights)) user.weights <- rep(1.0, nrow(mf))

    cluster<- attr(Terms, "specials")$cluster
    if (length(cluster)) {
        stop("a cluster() term is not valid")
    }
    
    # If there is start-stop data, then there needs to be an id
    #  also check that this is indeed a competing risks form of data.
    # Mark the first and last obs of each subject, as we need it later.
    #  Observations may not be in time order within a subject
    delay <- FALSE  # is there delayed entry?
    if (type=="mcounting") {
        if (is.null(id)) stop("(start, stop] data requires a subject id")
        else {
            index <- order(id, Y[,2]) # by time within id
            sorty <- Y[index,]
            first <- which(!duplicated(id[index]))
            last  <- c(first[-1] -1, length(id))
            if (any(sorty[-last, 3] != 0))
                stop("a subject has a transition before their last time point")
            delta <- c(sorty[-1,1], 0) - sorty[,2]
            if (any(delta[-last] !=0)) 
                stop("a subject has gaps in time")
            if (any(Y[first,1] > min(Y[,2]))) delay <- TRUE
            temp1 <- temp2 <- rep(FALSE, nrow(mf))
            temp1[index[first]] <- TRUE
            temp2[index[last]]  <- TRUE
            first <- temp1  #used later
            last <-  temp2
         }
    } else last <- rep(TRUE, nrow(mf))  

    if (missing(etype)) enum <- 1  #generate a data set for which endpoint?
    else {
        index <- match(etype, states)
        if (any(is.na(index)))
            stop("etype argument has a state that is not in the data")
        enum <- index[1]
        if (length(index) > 1) warning("only the first endpoint was used")
    }
    
    # make sure count, if present is syntactically valid
    if (!missing(count)) count <- make.names(count) else count <- NULL
    oname <- paste0(prefix, c("start", "stop", "status", "wt"))
        
    if (ncol(Y) ==2) {
        temp <- min(Y[,1], na.rm=TRUE)
        if (temp >0) zero <- 0
        else zero <- 2*temp -1  # a value less than any observed y
        Y <- cbind(zero, Y)  # add a start column
    }

    utime <- sort(unique(c(Y[,1:2])))  # all the unique times
    newtime <- matrix(findInterval(Y[,1:2], utime), ncol=2) 
    status <- Y[,3]

    newtime[status !=0, 2] <- newtime[status !=0,2] - .2
    Gsurv <- survfit(Surv(newtime[,1], newtime[,2], last & status==0) ~ istrat, 
                     se.fit=FALSE)
    if (delay) 
        Hsurv <- survfit(Surv(-newtime[,2], -newtime[,1], first) ~ istrat, 
                         se.fit =FALSE)
    status <- Y[, 3]

    # Do computations separately for each stratum
    stratfun <- function(i) {
        keep <- (istrat ==i)
        times <- sort(unique(Y[keep & status == enum, 2])) #unique event times 
        if (length(times)==0) return(NULL)  #no events in this stratum
        tdata <- mf[keep, -1, drop=FALSE]
        maxtime <- max(Y[keep, 2])

        Gtemp <- Gsurv[i]
        if (delay) {
            Htemp <- Hsurv[i]
            dtime <- rev(-Htemp$time[Htemp$n.event > 0])
            dprob <- c(rev(Htemp$surv[Htemp$n.event > 0])[-1], 1)
            ctime <- Gtemp$time[Gtemp$n.event > 0]
            cprob <- c(1, Gtemp$surv[Gtemp$n.event > 0]) 
            temp <- sort(unique(c(dtime, ctime))) # these will all be integers
            index1 <- findInterval(temp, dtime)
            index2 <- findInterval(temp, ctime)
            ctime <- utime[temp]
            cprob <- dprob[index1] * cprob[index2+1]  # G(t)H(t), eq 11 Geskus
        }
        else {
            ctime <- utime[Gtemp$time[Gtemp$n.event > 0]]
            cprob <- Gtemp$surv[Gtemp$n.event > 0]
        }
        
        ct2 <- c(ctime, maxtime)
        cp2 <- c(1.0, cprob)
        index <- findInterval(times, ct2, left.open=TRUE)
        index <- sort(unique(index))  # the intervals that were actually seen
        # times before the first ctime get index 0, those between 1 and 2 get 1
        ckeep <- rep(FALSE, length(ct2))
        ckeep[index] <- TRUE
        expand <- (Y[keep, 3] !=0 & Y[keep,3] != enum & last[keep]) #which rows to expand
        split <- .Call(Cfinegray, Y[keep,1], Y[keep,2], ct2, cp2, expand, 
                       c(TRUE, ckeep)) 
        tdata <- tdata[split$row,,drop=FALSE]
        tstat <- ifelse((status[keep])[split$row]== enum, 1, 0)


        tdata[[oname[1]]] <- split$start
        tdata[[oname[2]]] <- split$end
        tdata[[oname[3]]] <- tstat
        tdata[[oname[4]]] <- split$wt * user.weights[split$row]
        if (!is.null(count)) tdata[[count]] <- split$add
        tdata
    }

    if (max(istrat) ==1) result <- stratfun(1)
    else {
        tlist <- lapply(1:max(istrat), stratfun)
        result <- do.call("rbind", tlist)
    }

    rownames(result) <- NULL   #remove all the odd labels that R adds
    attr(result, "event") <- states[enum]
    result
}  
