# Automatically generated from the noweb directory
# residuals for a survfit object
residuals.survfit <- function(object, times, 
                              type=c("surv", "cumsurv", "cumhaz")){
    if (!inherits(object, "survfit"))
        stop("argument must be a survfit object")
    survfitms <- inherits(object, "survfitms")
    coxsurv <- inherits(object, "survfitcox")
    timefix <- (is.null(object$timefix) || object$timefix)
    type <- match.arg(type)
    type.int <- match(type, c("surv", "cumsurv", "cumhaz"))
    
    start.time <- object$start.time
    if (is.null(start.time)) start.time <- min(c(0, object$time))

    # check input arguments
    if (missing(times)) 
        stop ("the times argument is required")
    else {
        if (!is.numeric(times)) stop("times must be a numeric vector")
        times <- sort(unique(times))
        if (timefix) times <- aeqSurv(Surv(times))[,1]
    }
    
    # get the data
    # I always need the model frame
    if (coxsurv) {
        mf <- model.frame(object)
        if (is.null(object$y)) Y <- model.response(mf)
        else Y <- object$y
    }
    else {
        Call <- object$call
        formula <- formula(object)

        # the chunk below is shared with survfit.formula 
        # create a copy of the call that has only the arguments we want,
        #  and use it to call model.frame()
        indx <- match(c('formula', 'data', 'weights', 'subset','na.action',
                        'istate', 'id', 'cluster', "etype"), names(Call), nomatch=0)
        #It's very hard to get the next error message other than malice
        #  eg survfit(wt=Surv(time, status) ~1) 
        if (indx[1]==0) stop("a formula argument is required")
        temp <- Call[c(1, indx)]
        temp[[1L]] <- quote(stats::model.frame)
        mf <- eval.parent(temp)

        Terms <- terms(formula, c("strata", "cluster"))
        ord <- attr(Terms, 'order')
        if (length(ord) & any(ord !=1))
                stop("Interaction terms are not valid for this function")

        n <- nrow(mf)
        Y <- model.response(mf)
        if (inherits(Y, "Surv2")) {
            # this is Surv2 style data
            # if there are any obs removed due to missing, remake the model frame
            if (length(attr(mf, "na.action"))) {
                temp$na.action <- na.pass
                mf <- eval.parent(temp)
            }
            if (!is.null(attr(Terms, "specials")$cluster))
                stop("cluster() cannot appear in the model statement")
            new <- surv2data(mf)
            mf <- new$mf
            istate <- new$istate
            id <- new$id
            Y <- new$y
            if (anyNA(mf[-1])) { #ignore the response still found there
                if (missing(na.action)) temp <- get(getOption("na.action"))(mf[-1])
                else temp <- na.action(mf[-1])
                omit <- attr(temp, "na.action")
                mf <- mf[-omit,]
                Y <- Y[-omit]
                id <- id[-omit]
                istate <- istate[-omit]
            }                      
            n <- nrow(mf)
        }       
        else {
            if (!is.Surv(Y)) stop("Response must be a survival object")
            id <- model.extract(mf, "id")
            istate <- model.extract(mf, "istate")
        }
        if (n==0) stop("data set has no non-missing observations")

        casewt <- model.extract(mf, "weights")
        if (is.null(casewt)) casewt <- rep(1.0, n)
        else {
            if (!is.numeric(casewt)) stop("weights must be numeric")
            if (any(!is.finite(casewt))) stop("weights must be finite") 
            if (any(casewt <0)) stop("weights must be non-negative")
            casewt <- as.numeric(casewt)  # transform integer to numeric
        }

        if (!is.null(attr(Terms, 'offset'))) warning("Offset term ignored")

        cluster <- model.extract(mf, "cluster")
        temp <- untangle.specials(Terms, "cluster")
        if (length(temp$vars)>0) {
            if (length(cluster) >0) stop("cluster appears as both an argument and a model term")
            if (length(temp$vars) > 1) stop("can not have two cluster terms")
            cluster <- mf[[temp$vars]]
            Terms <- Terms[-temp$terms]
        }

        ll <- attr(Terms, 'term.labels')
        if (length(ll) == 0) X <- factor(rep(1,n))  # ~1 on the right
        else X <- strata(mf[ll])

        # Backwards support for the now-depreciated etype argument
        etype <- model.extract(mf, "etype")
        if (!is.null(etype)) {
            if (attr(Y, "type") == "mcounting" ||
                attr(Y, "type") == "mright")
                stop("cannot use both the etype argument and mstate survival type")
            if (length(istate)) 
                stop("cannot use both the etype and istate arguments")
            status <- Y[,ncol(Y)]
            etype <- as.factor(etype)
            temp <- table(etype, status==0)

            if (all(rowSums(temp==0) ==1)) {
                # The user had a unique level of etype for the censors
                newlev <- levels(etype)[order(-temp[,2])] #censors first
            }
            else newlev <- c(" ", levels(etype)[temp[,1] >0])
            status <- factor(ifelse(status==0,0, as.numeric(etype)),
                                 labels=newlev)

            if (attr(Y, 'type') == "right")
                Y <- Surv(Y[,1], status, type="mstate")
            else if (attr(Y, "type") == "counting")
                Y <- Surv(Y[,1], Y[,2], status, type="mstate")
            else stop("etype argument incompatable with survival type")
        }
        # end of shared code 
    }

    xlev <- levels(X)

    # Deal with ties
    if (is.null(Call$timefix) || Call$timefix) newY <- aeqSurv(Y) else newY <- Y

    # Find the clustering, if any
    if (!is.null(cluster)) {
        if (is.factor(cluster)) {
            clname <- levels(cluster)
            cluster <- as.integer(cluster)
        } else {
            clname  <- sort(unique(cluster))
            cluster <- match(cluster, clname)
        }
        ncluster <- length(clname)
    } 
    else if (!is.null(id)) {
        # treat the id as both identifier and clustering
        clname <- levels(id)
        cluster <- as.integer(id)
        ncluster <- length(clname)
    }
    else {
        # create our own clustering
        n <- nrow(Y)
        cluster <- 1:n
        ncluster <- n
        clname <- NULL
    }   

    ny <- ncol(newY)
    # order the data by group, which is the order we want for the output
    sort1 <- order(cluster)
    newY <- newY[sort1,]
    X <- X[sort1]
    casewt <- casewt[sort1]
    cluster <- cluster[sort1]

    # What survival curves was used?
    if (!coxsurv) {
        stype <- Call$stype
        if (is.null(stype)) stype <- 1
        ctype <- Call$ctype
        if (is.null(ctype)) ctype <- 1
        
        if (!survfitms) {
            resid <- rsurvpart1(newY, X, casewt, times,
                                type.int, stype, ctype, object)
            dimnames(resid) <- list(clname, times)
        }
        else stop("not finished")
    }
    resid
}
 rsurvpart1 <- function(Y, X, casewt, times,
         type, stype, ctype, fit) {
     
     ntime <- length(times)
     # We need an index of where the times vector matches the
     #  fit$time vector.  If the latter were (1, 24, 40) and
     #  times= (0, 5, 24, 50) the result would be 0, 1, 2, 3)
     # This has to be done separately for each curve though.
     # Also, we only want to match death times, not others
     if (is.null(fit$strata)) fitrow <- list(seq(along=fit$time))
     else {
         temp1 <- cumsum(fit$strata)
         temp2 <- c(1, temp1+1)
         fitrow <- lappy(1:length(j), function(i) seq(temp2[i], temp1[i]))
     }
     matchfun <- function(x, fit, index) {
         tt <- fit$time[index]  # subset to this curve
         deaths <- tt[fit$n.event[index] >0]
         i1 <- match(deaths, tt)
         i2 <- findInterval(x, deaths, left.open=FALSE) 
         # why pmax?  indices of 0 drop elements, so i1[i2] will be the
         #  wrong length
         ifelse(i2==0, 0, i1[pmax(1,i2)] + index[1] -1)  # index in fit
     }
         
     tindex <- matrix(0L, nrow(Y), length(times))
     for (i in 1:length(fitrow)) {
         yrow <- (as.integer(X) ==i)
         temp <- matchfun(times, fit, fitrow[[i]])
         tindex[yrow, ] <- rep(temp, each= length(yrow))
     }

     # repeat the indexing for Y onto fit$time
     ny <- ncol(Y)
     yindex <- matrix(0L, nrow(Y), length(times))
     event <- Y[,ny]
     if (ny==3) startindex <- yindex
     for (i in 1:length(fitrow)) {
         temp <- matchfun(Y[,ny-1], fit, fitrow[[i]])
         yrow <- (as.integer(X) ==i)
         yindex[yrow,] <- rep(temp, ncol(yindex))
         if (ny==3) {
             temp <- matchfun(Y[,1], fit, fitrow[[i]])
             startindex[yrow,] <- rep(temp, ncol(yindex))
         }
     }                    

     if (type==3 || stype==2) {
         if (ctype==1) {
             add1 <- (yindex <= tindex & rep(event, ntime))
             lsum <- unlist(lapply(fitrow, function(i) 
                          cumsum(fit$n.event[i]/fit$n.risk[i]^2)))
                 
             term1 <- c(0, 1/fit$n.risk)[ifelse(add1, 1+yindex, 1)]
             term2 <- c(0, lsum)[1+pmin(yindex, tindex)]
             if (ny==3) term3 <- c(0, lsum)[1 + startindex]

             if (ny==2) D <- matrix(term1 -  term2, ncol=ntime)
             else       D <- matrix(term1 + term3 - term2, ncol=ntime)

             # survival is exp(-H) so the derivative is a simple transform of D
             if (type==1) D <- -D* c(1,fit$surv)[1+ tindex]

             if (type==2) { # type 2 = AUC
                 auc <- unlist(lapply(fitrow, function(i) {
                     temp <- c(1, fit$surv[i])
                     cumsum(temp[-length(temp)] * diff(c(0, fit$time[i])))
                 }))

                 overtime <- (times[col(D)]- c(0,fit$time)[1+tindex]) # t - last event time
                 auc2 <- c(0, auc)[1 + tindex]  + overtime* c(1,fit$surv)[1+tindex] # A(0, t)
                 aterm1 <- -D * auc2
                 
                 lsum2 <- unlist(lapply(fitrow, function(i) 
                          cumsum(auc[i]*fit$n.event[i]/fit$n.risk[i]^2)))
                 aterm2 <- c(0, lsum2)[1 + pmin(yindex, tindex)]
                 aterm3 <- c(0, auc/fit$n.risk)[ifelse(add1, 1+yindex, 1)]

                 D <- matrix(aterm1 + aterm3 - aterm2, ncol=ntime)
             }
         } else {
             nevent <- fit$n.event/casewt[1]
             if (any(casewt != casewt[1])) {
                 # Have to reconstruct the number of obs with an event
                 temp <- lappy(seq(along=levels(X)), function(i) {
                     keep <- which(as.numeric(X) ==i)
                     dtime <- Y[keep & status==1, ny-1]
                     count <- table(dtime)
                     nindx <- (fitrow[[i]])[fit$n.event[fitrow[[i]]] > 0]
                     nevent[nindx] <- as.vector(count)
                     })
             }

             risk2 <- fit$n.risk
             ltemp <- fit$n.event/fit$n.risk^2
             for (i in which(nevent>1)) {
                 denom <- risk2[i] - fit$n.event*(0:(nevent[i]-1))/nevent[i]
                 risk2[i] <- mean(1/denom) # multiplier for the event
                 ltemp[i] <- fit$n.event* mean(1/denom^2)
             }

             add1 <- (yindex >= tindex & rep(event, ntime))
             lsum <- unlist(lapply(fitrow, function(i) cumsum(ltemp[i])))
             term1 <- c(0, 1/risk2)[ifelse(add1, 1+yindex, 1)]
             term2 <- c(0, lsum)[1+pmin(yindex, tindex)]
             if (ny==3) term3 <- c(0, lsum)[1 + startindex]

             if (ny==2) D <- matrix(term1 - term2, ncol=ntime)
             else D <- matrix(term1 + term3 - term2, ncol=ntime)

             if (type==1) D <- -D* c(0,fit$surv)[1+ tindex]
             else if (type==2) { #RMST
                 auc <- unlist(lapply(fitrow, function(i) {
                     temp <- c(1, fit$surv[i])
                     cumsum(temp[-length(temp)] * diff(c(0, fit$time[i])))
                 }))

                 overtime <- (times[col(D)]- c(0,fit$time)[1+tindex]) # t - last event time
                 auc2 <- c(0, auc)[1 + tindex]  + overtime* c(1,fit$surv)[1+tindex] # A(0, t)
                 aterm1 <- -D * auc2
                 
                 lsum2 <- unlist(lapply(fitrow, function(i) 
                          cumsum(auc[i]*ltemp)))
                 aterm2 <- c(0, lsum2)[1 + pmin(yindex, tindex)]
                 aterm3 <- c(0, auc/risk2)[ifelse(add1, 1+yindex, 1)]

                 D <- matrix(aterm1 + aterm3 - aterm2, ncol=ntime)
             }
         }
     } else {
         add1 <- (yindex <= tindex & rep(event, ntime))
         # dtemp avoids 1/0
         dtemp <- ifelse(fit$n.risk==fit$n.event, 0, 1/(fit$n.risk- fit$n.event))
         hsum <- unlist(lapply(fitrow, function(i) 
                      cumsum(dtemp[i]*fit$n.event[i]/fit$n.risk[i])))
             
         term1 <- c(0, dtemp)[ifelse(add1, 1+yindex, 1)]
         term2 <- c(0, hsum)[1+pmin(yindex, tindex)]
         if (ny==3) term3 <- c(0, hsum)[1 + startindex]

         if (ny==2) D <- matrix(term1 -  term2, ncol=ntime)
         else       D <- matrix(term1 + term3 - term2, ncol=ntime)

         browser()
         if (type==1) D <- -D* c(1,fit$surv)[1+ tindex]
         else {
             stop("not done")
         }
     }
     D
}
