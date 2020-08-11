# Automatically generated from the noweb directory
# residuals for a survfit object
residuals.survfit <- function(object, times, 
                              type=c("surv", "RMTS", "cumhaz"),
                              collapse=TRUE, weighted=FALSE){
    if (!inherits(object, "survfit"))
        stop("argument must be a survfit object")
    survfitms <- inherits(object, "survfitms")
    coxsurv <- inherits(object, "survfitcox")
    timefix <- (is.null(object$timefix) || object$timefix)
    type <- match.arg(type)
    type.int <- match(type, c("surv", "RMTS", "cumhaz"))
    
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
            dimnames(resid) <- list(NULL, times)
        }
        else {
            if (is.null(istate)) istate <-survcheck2(newY, id)$istate
            resid <- rsurvpart2(newY, X, casewt, istate, times, cluster,
                                type.int, object)
        }
    }
    else stop("coxph survival curves not yet available")
    
    if (weighted && any(casewt !=1)) resid <- resid*casewt
    if (collapse) {
        dd <- dim(resid)
        if (length(dd) ==3) {
            resid <- rowsum(matrix(resid, nrow=dd[1]), cluster, reorder=FALSE)
            dim(resid) <- c(length(resid)/(dd[2]*dd[3]), dd[2:3])
            dimnames(resid)[[1]] <- unique(cluster)
        }
        else {
            resid <- rowsum(resid, cluster, reorder=FALSE)
            rownames(resid) <- unique(cluster)
        }
    }
    resid
}
 rsurvpart1 <- function(Y, X, casewt, times,
         type, stype, ctype, fit) {
     
     ntime <- length(times)
     etime <- (fit$n.event >0)
     event <- (Y[,ny] >0)
     # 
     #  Create a list whose first element contains the location of
     #   the death times in curve 1, second element for curve 2, etc.
     #  
     if (is.null(fit$strata)) fitrow <- list(which(etime))
     else {
         temp1 <- cumsum(fit$strata)
         temp2 <- c(1, temp1+1)
         fitrow <- lappy(1:length(fit$strata), function(i) {
             indx <- seq(temp2[i], temp1[i])
             indx[etime[indx]] # keep the death times
         }) 
     }
 
     # for each time x, the index of the last death time which is <=x.
     #  0 if x is before the first death time
     matchfun <- function(x, fit, index) {
         dtime <- fit$time[index]  # subset to this curve
         i2 <- findInterval(x, dtime, left.open=FALSE)
         c(0, index)[i2 +1]
     }
     
     # output matrix D will have one row per observation, one col for each
     #  reporting time. tindex and yindex have the same dimension as D.
     # tindex points to the last death time in fit which
     #  is <= the reporting time.  (If there is only 1 curve, each col of
     #  tindex will be a repeat of the same value.)
     tindex <- matrix(0L, nrow(Y), length(times))
     for (i in 1:length(fitrow)) {
         yrow <- (as.integer(X) ==i)
         temp <- matchfun(times, fit, fitrow[[i]])
         tindex[yrow, ] <- rep(temp, each= length(yrow))
     }

     # repeat the indexing for Y onto fit$time.  Each row of yindex points
     #  to the last row of fit with death time <= Y[,ny]
     ny <- ncol(Y)
     yindex <- matrix(0L, nrow(Y), length(times))
     event <- (Y[,ny] >0)
     if (ny==3) startindex <- yindex
     for (i in 1:length(fitrow)) {
         yrow <- (as.integer(X) ==i)  # rows of Y for this curve
         temp <- matchfun(Y[yrow,ny-1], fit, fitrow[[i]])
         yindex[yrow,] <- rep(temp, ncol(yindex))
         if (ny==3) {
             temp <- matchfun(Y[yrow,1], fit, fitrow[[i]])
             startindex[yrow,] <- rep(temp, ncol(yindex))
         }
     }                    

     # Now do the work
     if (type==3 || stype==2) {
         if (ctype==1) {
             add1 <- (yindex <= tindex & rep(event, ntime))
             lsum <- unlist(lapply(fitrow, function(i) 
                          cumsum(fit$n.event[i]/fit$n.risk[i]^2)))
                 
             term1 <- c(0, 1/fit$n.risk)[ifelse(add1, 1+yindex, 1)]
             term2 <- c(0, lsum)[1+pmin(yindex, tindex)]
             if (ny==3) term3 <- c(0, lsum)[1 + pmin(startindex, tindex)]

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
         # dtemp avoids 1/0.  (When this occurs the influence is 0, since
         #  the curve has dropped to zero; and this avoids Inf-Inf in term1-term2).
         dtemp <- ifelse(fit$n.risk==fit$n.event, 0, 1/(fit$n.risk- fit$n.event))
         hsum <- unlist(lapply(fitrow, function(i) 
                      cumsum(dtemp[i]*fit$n.event[i]/fit$n.risk[i])))
             
         term1 <- c(0, dtemp)[ifelse(add1, 1+yindex, 1)]
         term2 <- c(0, hsum)[1+pmin(yindex, tindex)]
         if (ny==3) term3 <- c(0, hsum)[1 + startindex]

         if (ny==2) D <- matrix(term1 -  term2, ncol=ntime)
         else       D <- matrix(term1 + term3 - term2, ncol=ntime)

         if (type==1) D <- -D* c(1,fit$surv)[1+ tindex]
         else if (type==2){
             auc <- unlist(lapply(fitrow, function(i) {
                 temp <- c(1, fit$surv[i])
                 cumsum(temp[-length(temp)] * diff(c(0, fit$time[i])))
             }))

             overtime <- (times[col(D)]- c(0,fit$time)[1+tindex]) # t - last event time
             auc2 <- c(0, auc)[1 + tindex]  + overtime* c(1,fit$surv)[1+tindex] # A(0, t)
             aterm1 <- -D * auc2
             
             lsum2 <- unlist(lapply(fitrow, function(i) 
                      cumsum(auc[i]*(dtemp[i]*fit$n.event[i]/fit$n.risk[i]))))
             aterm2 <- c(0, lsum2)[1 + pmin(yindex, tindex)]
             aterm3 <- c(0, auc/fit$n.risk)[ifelse(add1, 1+yindex, 1)]

             D <- matrix(aterm1 + aterm3 - aterm2, ncol=ntime)
         }
     }
     D
}
rsurvpart2 <- function(Y, X, casewt, istate, times, cluster, type, fit) {
    ny <- ncol(Y)
    ntime <- length(times)
    nstate <- length(fit$states)
    
    # ensure that Y, istate, and fit all use the same set of states
    states <- fit$states
    if (!identical(attr(Y, "states"), fit$states)) {
        map <- match(attr(Y, "states"), fit$states)
        Y[,ny] <- c(0, map)[1+ Y[,ny]]    # 0 = censored
    }
    if (is.character(istate)) istate <- factor(istate)
    if (is.factor(istate)) {
        if (!identical(levels(istate), fit$states)) {
            map <- match(levels(istate), fit$states)
            if (any(is.na(map))) stop ("invalid levels in istate")
            istate <- map[istate]
        }
    } # istate is numeric, we take what we get and hope it's right

    # Compute the initial leverage
    inf0 <- NULL
    if (is.null(fit$call$p0)) { #p0 was not supplied by the user
        inf0 <- matrix(0., nrow=nrow(Y), ncol=nstate)
        i0fun <- function(i, fit, inf0) {
            # reprise algorithm in survfitCI
            p0 <- fit$p0
            t0 <- fit$time[1]
            if (ny==2) at.zero <- which(as.numeric(X) ==i)
            else       
                at.zero <- which(as.numeric(X) ==i &
                          (Y[,1] < t0 & Y[,2] >= t0))
            for (j in 1:nstate) {
                inf0[at.zero, j] <- (ifelse(istate[at.zero]==states[j], 1, 0) -
                                     p0[j])/sum(casewt[at.zero])
            }
            inf0
        }

        if (is.null(fit$strata)) inf0 <- i0fun(1, fit, inf0)
        else for (i in 1:length(levels(X)))
            inf0 <- ifun0(i, fit[i], inf0)
    }
    
    # collapse redundant rows in Y, for efficiency
    if (ny==3 && any(duplicated(cluster))) {
        ord <- order(cluster, X, istate, Y[,1])
        cfit <- .Call(Ccollapse, Y, X, istate, cluster, casewt, ord -1L) 
        if (nrow(cfit) < .8*length(X))  {
            # shrinking the data by 20% is worth it
            temp <- Y[ord,]
            Y <- cbind(temp[cfit[,1], 1], temp[cfit[2], 2:3])
            X <- X[cfit[,1]]
            istate <- istate[cfit[1,]]
            cluster <- cluster[cfit[1,]]
        }       
    }

    fit <- survfit0(fit)  # package the initial state into the picture
  
    # This next block is identical to the one in rsurvpart1, more comments are
    #  there
    etime <- (rowSums(fit$n.event) >0)
    event <- (Y[,ny] >0)
    # 
    #  Create a list whose first element contains the location of
    #   the death times in curve 1, second element for curve 2, etc.
    #  
    if (is.null(fit$strata)) fitrow <- list(which(etime))
    else {
        temp1 <- cumsum(fit$strata)
        temp2 <- c(1, temp1+1)
        fitrow <- lappy(1:length(fit$strata), function(i) {
            indx <- seq(temp2[i], temp1[i])
            indx[etime[indx]] # keep the death times
        }) 
    }
 
     # for each time x, the index of the last death time which is <=x.
     #  0 if x is before the first death time
     matchfun <- function(x, fit, index) {
         dtime <- fit$time[index]  # subset to this curve
         i2 <- findInterval(x, dtime, left.open=FALSE)
         c(0, index)[i2 +1]
     }
     

     if (type==3) {
         # output matrix D will have one row per observation, one col for each
         #  reporting time. tindex and yindex have the same dimension as D.
         # tindex points to the last death time in fit which
         #  is <= the reporting time.  (If there is only 1 curve, each col of
         #  tindex will be a repeat of the same value.)
         tindex <- matrix(0L, nrow(Y), length(times))
         for (i in 1:length(fitrow)) {
             yrow <- (as.integer(X) ==i)
             temp <- matchfun(times, fit, fitrow[[i]])
             tindex[yrow, ] <- rep(temp, each= length(yrow))
         }

         # repeat the indexing for Y onto fit$time.  Each row of yindex points
         #  to the last row of fit with death time <= Y[,ny]
         ny <- ncol(Y)
         yindex <- matrix(0L, nrow(Y), length(times))
         event <- (Y[,ny] >0)
         if (ny==3) startindex <- yindex
         for (i in 1:length(fitrow)) {
             yrow <- (as.integer(X) ==i)  # rows of Y for this curve
             temp <- matchfun(Y[yrow,ny-1], fit, fitrow[[i]])
             yindex[yrow,] <- rep(temp, ncol(yindex))
             if (ny==3) {
                 temp <- matchfun(Y[yrow,1], fit, fitrow[[i]])
                 startindex[yrow,] <- rep(temp, ncol(yindex))
             }
         }                    

         dstate <- Y[,ncol(Y)]
         istate <- as.numeric(istate)
         ntrans <- ncol(fit$cumhaz)  # the number of possible transitions
         D <- array(0, dim=c(nrow(Y), ntime, ntrans))

         scount <- table(istate[dstate!=0], dstate[dstate!=0]) # observed transitions
         state1 <- row(scount)[scount>0]
         state2 <- col(scount)[scount>0]
         temp <- paste(state1, state2, sep='.')
         if (!identical(temp, colnames(fit$cumhaz))) stop("setup error")

         for (k in length(state1)) {
             e2 <- Y[,ny] == state2[k]
             add1 <- (yindex <= tindex & rep(e2, ntime))
             lsum <- unlist(lapply(fitrow, function(i) 
                      cumsum(fit$n.event[i,k]/fit$n.risk[i,k]^2)))
             
             term1 <- c(0, 1/fit$n.risk[,k])[ifelse(add1, 1+yindex, 1)]
             term2 <- c(0, lsum)[1+pmin(yindex, tindex)]
             if (ny==3) term3 <- c(0, lsum)[1 + startindex]

             if (ny==2) D[,,k] <- matrix(term1 -  term2, ncol=ntime)
             else       D[,,k] <- matrix(term1 + term3 - term2, ncol=ntime)
         }
     } else {
         ndeath <- sum(fit$time <= max(times) & etime)  # number of unique event times
         # Expand Y
         if (ny==2) split <- .Call(rep(0., nrow(Y)), Y[,1], times)
         else split <- .Call(Csurvsplit, Y[,1], Y[,2], times)
         X <- X[split$row]
         casewt <- casewt[split$row]
         istate <- istate[split$row]
         Y <- cbind(split$start, split$end, 
                     ifelse(split$censor, 0, Y[split$row,ny]))
         ny <- 3

         # Create a vector containing the index of each end time into the fit object
         yindex <- ystart <- double(nrow(Y))
         for (i in 1:length(fitrow)) {
             yrow <- (as.integer(X) ==i)  # rows of Y for this curve
             yindex[yrow] <- matchfun(Y[yrow, 2], fit, fitrow[[i]])
             ystart[yrow] <- matchfun(Y[yrow, 1], fit, fitrow[[i]])
         }
         # And one indexing the reporting times into fit
         tindex <- matrix(0L, nrow=length(fitrow), ncol=ntime)
         for (i in 1:length(fitrow)) {
             tindex[i,] <- matchfun(times, fit, fitrow[[i]])
         }

         # Create the array of C matrices
         cmat <- array(0, dim=c(nstate, nstate, ndeath)) # max(i2) = ndeath, by design
         Hmat <- cmat

         # We only care about observations that had a transition; any transitions
         #  after the last reporting time are not relevant
         transition <- (Y[,ny] !=0 & Y[,ny] != istate &
                        Y[,ny-1] <= max(times)) # obs that had a transition
         i2 <- match(yindex, sort(unique(yindex)))  # which C matrix this obs goes to
         i2 <- i2[transition]
         from <- as.numeric(istate[transition])  # from this state
         to   <- Y[transition, ny]   # to this state
         nrisk <- fit$n.risk[cbind(yindex[transition], from)]  # number at risk
         wt <- casewt[transition]
         for (i in seq(along=from)) {
             j <- c(from[i], to[i])
             haz <- wt[i]/nrisk[i]
             cmat[from[i], j, i2[i]] <- cmat[from[i], j, i2[i]] + c(-haz, haz)
         }
         for (i in 1:ndeath) Hmat[,,i] <- cmat[,,i] + diag(nstate)

         # The transformation matrix H(t) at time t  is cmat[,,t] + I
         # Create the set of W matrices.
         # 
         dindex <- which(etime & fit$time <= max(times))
         Wmat <- array(0, dim=c(nstate, nstate, ndeath))
         for (i in ndeath:1) {
             j <- match(dindex[i], tindex, nomatch=0) 
             if (j > 0) {
                 # this death matches one of the reporting times
                 Wmat[,,i] <- diag(nstate)
             } 
             else Wmat[,,i] <- Hmat[,,i+1] %*% Wmat[,,i+1]
         }
         iterm <- array(0, dim=c(nstate, nstate, ndeath)) # term in equation
         itemp <- matrix(0, nstate, nstate)  # cumulative sum, temporary
         isum  <- isum2 <- iterm  # cumulative sum
         for (i in 1:ndeath) {
             j <- dindex[i]
             n0 <- ifelse(fit$n.risk[j,] ==0, 1, fit$n.risk[j,]) # avoid 0/0
             iterm[,,i] <- ((fit$pstate[j-1,]/n0) * cmat[,,i]) %*% Wmat[,,i]
             itemp <- itemp + iterm[,,i]
             isum[,,i] <- itemp
             j <- match(dindex[i], tindex, nomatch=0)
             if (j>0) itemp <- matrix(0, nstate, nstate)  # reset
             isum2[,,i] <- itemp
         }

         # We want to add isum[state,, entry time] - isum[state,, exit time] for
         #  each subject, and for those with an a:b transition there will be an 
         #  additional vector with -1, 1 in the a and b position.
         i1 <- match(ystart, sort(unique(yindex)), nomatch=0) # start at 0 gives 0
         i2 <- match(yindex, sort(unique(yindex)))
         D <- matrix(0., nrow(Y), nstate)
         keep <- (Y[,2] <= max(times))  # any intervals after the last reporting time
                                         # will have 0 influence
         for (i in which(keep)) {
             if (Y[i,3] !=0 && istate[i] != Y[i,3]) {
                 z <- fit$pstate[yindex[i]-1, istate[i]]/fit$n.risk[yindex[i], istate[i]]
                 temp <- double(nstate)
                 temp[istate[i]] = -z
                 temp[Y[i,3]]    =  z
                 temp <- temp %*% Wmat[,,i2[i]] - isum[istate[i],,i2[i]]
                 if (i1[i] >0) temp <- temp + isum2[istate[i],, i1[i]]
                 D[i,] <- temp
             }
             else {
                 if (i1[i] >0) D[i,] = isum2[istate[i],,i1[i]] - isum[istate[i],, i2[i]]
                 else  D[i,] =  -isum[istate[i],, i2[i]]
             }
         }
         Dsave <- D
         if (!is.null(inf0)) {
             # add in the initial influence, to the first row of each obs
             #   (inf0 was created on unsplit data)
             j <- which(!duplicated(split$row))
             D[j,] <- D[j,] + (inf0%*% Hmat[,,1] %*% Wmat[,,1])
         }
         if (ntime > 1) {
             interval <- findInterval(yindex, tindex, left.open=TRUE)
             D2 <- array(0., dim=c(dim(D), ntime))
             D2[interval==0,,1] <- D[interval==0,]
             for (i in 1:(ntime-1)) {
                 D2[interval==i,,i+1] = D[interval==i,]
                 browser()
                 j <- tindex[i]
                 D2[,,i+1] = D2[,,i+1] + D2[,,i] %*% (Hmat[,,j] %*% Wmat[,,j])
             } 
         browser()
             D <- D2
         }
         else browser()

         # undo any artificial split
         if (any(duplicated(split$row))) {
             if (ntime==1) D <- rowsum(D, split$row)
             else {
                 # rowsums has to be fooled
                 temp <- rowsum(matrix(D, ncol=(nstate*ntime)), split$row)
                 # then undo it
                 D <- array(temp, dim=c(nrow(temp), nstate, ntime))
             }
         }
         D
     }
     D
}
