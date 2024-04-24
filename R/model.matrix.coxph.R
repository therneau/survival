# Automatically generated from the noweb directory
# In internal use "data" will often be an already derived model frame.
#  We detect this via it having a terms attribute.
model.matrix.coxph <- function(object, data=NULL, 
                               contrast.arg=object$contrasts, ...) {
    # 
    # If the object has an "x" component, return it, unless a new
    #   data set is given
    if (is.null(data) && !is.null(object[['x']])) 
        return(object[['x']]) #don't match "xlevels"

    Terms <- delete.response(object$terms)
    if (is.null(data)) mf <- stats::model.frame(object)
    else {
        if (is.null(attr(data, "terms")))
            mf <- stats::model.frame(Terms, data, xlev=object$xlevels)
        else mf <- data  #assume "data" is already a model frame
    }

    cluster <- attr(Terms, "specials")$cluster
    if (length(cluster)) {
        temp <- untangle.specials(Terms, "cluster")
        dropterms <- temp$terms
    }
    else dropterms <- NULL
    
    strats <- attr(Terms, "specials")$strata
    hasinteractions <- FALSE
    if (length(strats)) {
        stemp <- untangle.specials(Terms, 'strata', 1)
        if (length(stemp$vars)==1) strata.keep <- mf[[stemp$vars]]
        else strata.keep <- strata(mf[,stemp$vars], shortlabel=TRUE)
        istrat <- as.integer(strata.keep)

        for (i in stemp$vars) {  #multiple strata terms are allowed
            # The factors attr has one row for each variable in the frame, one
            #   col for each term in the model.  Pick rows for each strata
            #   var, and find if it participates in any interactions.
            if (any(attr(Terms, 'order')[attr(Terms, "factors")[i,] >0] >1))
                hasinteractions <- TRUE  
        }
        if (!hasinteractions) dropterms <- c(dropterms, stemp$terms) 
    } else istrat <- NULL


    if (length(dropterms)) {
        Terms2 <- Terms[-dropterms]
        X <- model.matrix(Terms2, mf, constrasts.arg=contrast.arg)
        # we want to number the terms wrt the original model matrix
        temp <- attr(X, "assign")
        shift <- sort(dropterms)
        for (i in seq(along.with=shift))
            temp <- temp + 1*(shift[i] <= temp)
        attr(X, "assign") <- temp 
    }
    else X <- model.matrix(Terms, mf, contrasts.arg=contrast.arg)

    # drop the intercept after the fact, and also drop strata if necessary
    Xatt <- attributes(X) 
    if (hasinteractions) adrop <- c(0, untangle.specials(Terms, "strata")$terms)
    else adrop <- 0
    xdrop <- Xatt$assign %in% adrop  #columns to drop (always the intercept)
    X <- X[, !xdrop, drop=FALSE]
    attr(X, "assign") <- Xatt$assign[!xdrop]
    attr(X, "contrasts") <- Xatt$contrasts
    X
}
model.frame.coxph <- function(formula, ...) {
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset", "weights",
                          "id", "cluster", "istate"), 
                        names(dots), 0)] 
    # If nothing has changed and the coxph object had a model component,
    #   simply return it.
    if (length(nargs) ==0  && !is.null(formula$model)) return(formula$model)
    else {
        # Rebuild the original call to model.frame
        Terms <- terms(formula)
        fcall <- formula$call
        indx <- match(c("formula", "data", "weights", "subset", "na.action",
                        "cluster", "id", "istate"),
                  names(fcall), nomatch=0) 
        if (indx[1] ==0) stop("The coxph call is missing a formula!")
   
        temp <- fcall[c(1,indx)]  # only keep the arguments we wanted
        temp[[1]] <- quote(stats::model.frame)  # change the function called
        temp$xlev <- formula$xlevels  # this will turn strings to factors
        temp$formula <- Terms   #keep the predvars attribute
        # Now, any arguments that were on this call overtake the ones that
        #  were in the original call.  
        if (length(nargs) >0)
            temp[names(nargs)] <- nargs

        # Make "tt" visible for coxph formulas, 
        if (!is.null(attr(temp$formula, "specials")$tt)) {
            coxenv <- new.env(parent= environment(temp$formula))
            assign("tt", function(x) x, envir=coxenv)
            environment(temp$formula) <- coxenv
        }

        # The documentation for model.frame implies that the environment arg
        #  to eval will be ignored, but if we omit it there is a problem.
        if (is.null(environment(formula$terms))) 
            mf <- eval(temp, parent.frame())
        else mf <- eval(temp, environment(formula$terms), parent.frame())

        if (!is.null(attr(formula$terms, "dataClasses")))
            .checkMFClasses(attr(formula$terms, "dataClasses"), mf)
       
        if (is.null(attr(Terms, "specials")$tt)) return(mf)
        else {
            # Do time transform
            tt <- eval(formula$call$tt)
            Y <- aeqSurv(model.response(mf))
            strats <- attr(Terms, "specials")$strata
            if (length(strats)) {
                stemp <- untangle.specials(Terms, 'strata', 1)
                if (length(stemp$vars)==1) strata.keep <- mf[[stemp$vars]]
                else strata.keep <- strata(mf[,stemp$vars], shortlabel=TRUE)
                istrat <- as.numeric(strata.keep)
            }
            
            timetrans <- untangle.specials(Terms, 'tt')
            ntrans <- length(timetrans$terms)

            if (is.null(tt)) {
                tt <- function(x, time, riskset, weights){ #default to O'Brien's logit rank
                    obrien <- function(x) {
                        r <- rank(x)
                        (r-.5)/(.5+length(r)-r)
                    }
                    unlist(tapply(x, riskset, obrien))
                }
            }
            if (is.function(tt)) tt <- list(tt)  #single function becomes a list
                
            if (is.list(tt)) {
                if (any(!sapply(tt, is.function))) 
                    stop("The tt argument must contain function or list of functions")
                if (length(tt) != ntrans) {
                    if (length(tt) ==1) {
                        temp <- vector("list", ntrans)
                        for (i in 1:ntrans) temp[[i]] <- tt[[1]]
                        tt <- temp
                    }
                    else stop("Wrong length for tt argument")
                }
            }
            else stop("The tt argument must contain a function or list of functions")

            if (ncol(Y)==2) {
                if (length(strats)==0) {
                    sorted <- order(-Y[,1], Y[,2])
                    newstrat <- rep.int(0L, nrow(Y))
                    newstrat[1] <- 1L
                    }
                else {
                    sorted <- order(istrat, -Y[,1], Y[,2])
                    #newstrat marks the first obs of each strata
                    newstrat <-  as.integer(c(1, 1*(diff(istrat[sorted])!=0))) 
                    }
                if (storage.mode(Y) != "double") storage.mode(Y) <- "double"
                counts <- .Call(Ccoxcount1, Y[sorted,], 
                                as.integer(newstrat))
                tindex <- sorted[counts$index]
            }
            else {
                if (length(strats)==0) {
                    sort.end  <- order(-Y[,2], Y[,3])
                    sort.start<- order(-Y[,1])
                    newstrat  <- c(1L, rep(0, nrow(Y) -1))
                }
                else {
                    sort.end  <- order(istrat, -Y[,2], Y[,3])
                    sort.start<- order(istrat, -Y[,1])
                    newstrat  <- c(1L, as.integer(diff(istrat[sort.end])!=0))
                }
                if (storage.mode(Y) != "double") storage.mode(Y) <- "double"
                counts <- .Call(Ccoxcount2, Y, 
                                as.integer(sort.start -1L),
                                as.integer(sort.end -1L), 
                                as.integer(newstrat))
                tindex <- counts$index
            }
            Y <- Surv(rep(counts$time, counts$nrisk), counts$status)
            type <- 'right'  # new Y is right censored, even if the old was (start, stop]

            mf <- mf[tindex,]
            istrat <- rep(1:length(counts$nrisk), counts$nrisk)
            weights <- model.weights(mf)
            if (!is.null(weights) && any(!is.finite(weights)))
                stop("weights must be finite") 
            id <- model.extract(mf, "id")   # update the id and/or cluster, if present
            cluster <- model.extract(mf, "cluster")

            tcall <- attr(Terms, 'variables')[timetrans$terms+2]
            pvars <- attr(Terms, 'predvars')
            pmethod <- sub("makepredictcall.", "", as.vector(methods("makepredictcall")))
            for (i in 1:ntrans) {
                newtt <- (tt[[i]])(mf[[timetrans$vars[i]]], Y[,1], istrat, weights)
                mf[[timetrans$vars[i]]] <- newtt
                nclass <- class(newtt)
                if (any(nclass %in% pmethod)) { # It has a makepredictcall method
                    dummy <- as.call(list(as.name(class(newtt)[1]), tcall[[i]][[2]]))
                    ptemp <- makepredictcall(newtt, dummy)
                    pvars[[timetrans$terms[i]+2]] <- ptemp
                }
            }
            attr(Terms, "predvars") <- pvars
            mf[[".strata."]] <- istrat
            return(mf)
        }
    }
}
