# Automatically generated from all.nw using noweb
# In internal use "data" will often be an already derived model frame.
#  We detect this via it having a terms attribute.
model.matrix.coxph <- function(object, data=NULL, 
                               contrast.arg=object$contrasts, ...) {
    # 
    # If the object has an "x" component, return it, unless a new
    #   data set is given
    if (is.null(data) && !is.null(object[['x']])) 
        return(object[['x']]) #object$x might fetch "xlevels" instead

    Terms <- delete.response(object$terms)
    if (is.null(data)) mf <- model.frame(object)
    else {
        if (is.null(attr(data, "terms")))
            mf <- model.frame(Terms, data, xlev=object$xlevels)
        else mf <- data  #assume "data" is already a model frame
    }

    cluster <- attr(Terms, "specials")$cluster
    if (length(cluster)) {
        temp <- untangle.specials(Terms, "cluster")
        dropterms <- temp$terms
    }
    else dropterms <- NULL
    
    attr(Terms, "intercept") <- TRUE
    adrop <- 0  #levels of "assign" to be dropped; 0= intercept
    stemp <- untangle.specials(Terms, 'strata', 1)
    if (length(stemp$vars) > 0) {  #if there is a strata statement
        hasinteractions <- FALSE
        for (i in stemp$vars) {  #multiple strata terms are allowed
            # The factors att has one row for each variable in the frame, one
            #   col for each term in the model.  Pick rows for each strata
            #   var, and find if it participates in any interactions.
            if (any(attr(Terms, 'order')[attr(Terms, "factors")[i,] >0] >1))
                hasinteractions <- TRUE  
            }
        if (!hasinteractions) 
            dropterms <- c(dropterms, stemp$terms)
        else adrop <- c(0, match(stemp$var, colnames(attr(Terms, 'factors'))))
    }

    if (length(dropterms)) {
        temppred <- attr(terms, "predvars")
        Terms2 <- Terms[ -dropterms]
        if (!is.null(temppred)) {
            # subscripting a Terms object currently drops predvars, in error
            attr(Terms2, "predvars") <- temppred[-(1+dropterms)] # "Call" object
        }
        X <- model.matrix(Terms2, mf, constrasts=contrast.arg)
        # we want to number the terms wrt the original model matrix
        # Do not forget the intercept, which will be a zero
        renumber <- match(colnames(attr(Terms2, "factors")), 
                          colnames(attr(Terms,  "factors")))
        attr(X, "assign") <- c(0, renumber)[1+attr(X, "assign")]
    }
    else X <- model.matrix(Terms, mf, contrasts=contrast.arg)

    # drop the intercept after the fact, and also drop strata if necessary
    Xatt <- attributes(X) 
    xdrop <- Xatt$assign %in% adrop  #columns to drop (always the intercept)
    X <- X[, !xdrop, drop=FALSE]
    attr(X, "assign") <- Xatt$assign[!xdrop]
    #if (any(adrop>0)) attr(X, "contrasts") <- Xatt$contrasts[-adrop]
    #else attr(X, "contrasts") <- Xatt$contrasts
    attr(X, "contrasts") <- Xatt$contrasts
    X
}
model.frame.coxph <- function(formula, ...) {
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset", "weights"), 
                        names(dots), 0)] 
    # If nothing has changed and the coxph object had a model component,
    #   simply return it.
    if (length(nargs) ==0  && !is.null(formula$model)) return(formula$model)
    else {
        # Rebuild the original call to model.frame
        Terms <- terms(formula)
        fcall <- formula$call
        indx <- match(c("formula", "data", "weights", "subset", "na.action"),
                  names(fcall), nomatch=0) 
        if (indx[1] ==0) stop("The coxph call is missing a formula!")
   
        temp <- fcall[c(1,indx)]  # only keep the arguments we wanted
        temp[[1]] <- as.name('model.frame')  # change the function called
        temp$xlev <- formula$xlevels
        temp$formula <- Terms   #keep the predvars attribute
        # Now, any arguments that were on this call overtake the ones that
        #  were in the original call.  
        if (length(nargs) >0)
            temp[names(nargs)] <- nargs

        # The documentation for model.frame implies that the environment arg
        #  to eval will be ignored, but if we omit it there is a problem.
        if (is.null(environment(formula$terms))) 
            mf <- eval(temp, parent.frame())
        else mf <- eval(temp, environment(formula$terms), parent.frame())

        if (!is.null(attr(formula$terms, "dataClasses")))
            .checkMFClasses(attr(formula$terms, "dataClasses"), mf)
       
       if (!is.null(attr(Terms, "specials")$tt)) {
          # Do time transform
          tt <- eval(formula$call$tt)
          Y <- model.response(mf)
          strats <- attr(Terms, "specials")$strata
          if (length(strats)) {
              stemp <- untangle.specials(Terms, 'strata', 1)
              if (length(stemp$vars)==1) strata.keep <- mf[[stemp$vars]]
              else strata.keep <- strata(mf[,stemp$vars], shortlabel=TRUE)
              strats <- as.numeric(strata.keep)
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
                  sorted <- order(strats, -Y[,1], Y[,2])
                  #newstrat marks the first obs of each strata
                  newstrat <-  as.integer(c(1, 1*(diff(strats[sorted])!=0))) 
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
                  sort.end  <- order(strats, -Y[,2], Y[,3])
                  sort.start<- order(strats, -Y[,1])
                  newstrat  <- c(1L, as.integer(diff(strats[sort.end])!=0))
              }
              if (storage.mode(Y) != "double") storage.mode(Y) <- "double"
              counts <- .Call(Ccoxcount2, Y, 
                              as.integer(sort.start -1L),
                              as.integer(sort.end -1L), 
                              as.integer(newstrat))
              tindex <- counts$index
          }
          mf <- mf[tindex,]
          Y <- Surv(rep(counts$time, counts$nrisk), counts$status)
          type <- 'right'  # new Y is right censored, even if the old was (start, stop]
          strats <- rep(1:length(counts$nrisk), counts$nrisk)
          weights <- model.weights(mf)
          for (i in 1:ntrans) 
              mf[[timetrans$var[i]]] <- (tt[[i]])(mf[[timetrans$var[i]]], Y[,1], strats, 
                                                 weights)
          mf[[".strata."]] <- strats
       }
       mf
    }
}
