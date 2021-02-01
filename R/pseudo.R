#
# Function for pseudo-values
#
pseudo <- function(fit, times, type, addNA=TRUE, data.frame=FALSE,
                   minus1=FALSE, ...) {
    if (!inherits(fit, "survfit"))
        stop("fit argument must be a survfit object")
    if (inherits(fit, "survfit.coxph"))
        stop("psuedo values not defined for a coxph survival curve")

    if (missing(type)) type <- "pstate"
    else {
        # allow various aliases for the type
        temp <- c("pstate", "cumhaz", "sojourn", "survival",
                              "chaz", "rmst", "rmts", "auc")
        type <- match.arg(casefold(type), temp)
        itemp <-  c(1,2,3,1,2,3,3,3)[match(type, temp)]
        type <- c("pstate", "cumhaz", "auc")[itemp]
    }

    # all the real work is done by the residuals function
    rfit <- survresid.fit(fit, times=times, type=type, ...)
    resid <- rfit$residuals
    curve <- rfit$curve
    id    <- rfit$id
    idname <- rfit$idname

    ddr <- dim(resid)
    ddf <- dim(fit)

    if (is.null(fit$strata)) ncurve <- 1L  
    else ncurve <- length(fit$strata)  # ncurve will also = max(curve)

    ntime <- length(times)
    nstate <- length(fit$states)

    minus <- ifelse(minus1, 1, 0)  # use n-1 or n to multiply

    # get the estimates
    # retain the dimension of the residuals, which will have subjects first,
    #  timepoints second and (state or cumhaz) last.
    # The residuals have a row for each subject.  To recenter them for a
    #  pseudovalue, however, we need to match each one to the curve from
    #  whence it came.
    # If the data has been collapsed, n is the number of subjects after
    #  the collapse, otherwise the number of observations.
    nn <- as.vector(table(curve)) - minus

    if (type == "pstate") {
        temp <- summary(fit, times=times, extend=TRUE)
        if (!is.null(temp$pstate)) { # multi-state
            # pstate will have times, curve, state,  as the order
            # resid will have subject, times, state
            yhat <-  array(temp$pstate, dim=c(ntime, ncurve, nstate))
            yhat <-  aperm(yhat, c(2,1,3))
            # yhat[curve,,] is the trick: resid has first dimension of n =
            #  number of subjects, yhat has 1 row per curve, yhat[curve,,]
            #  will have n rows, each pulled from the correct curve
            ptemp <- yhat[curve,,] + nn[curve]*resid
        }
        else {
            yhat <- matrix(temp$surv, nrow = ncurve, byrow=TRUE)
            ptemp <- yhat[curve,] + nn[curve]*resid
        }
    }       
    else if (type=="cumhaz") {
        temp <- summary(fit, times=times, extend=TRUE)$cumhaz
        if (is.matrix(temp)) { # multi-state
            nhaz <- ncol(fit$cumhaz)
            yhat <- array(temp, dim=c(ntime, ncurve, nhaz))
            yhat <- aperm(yhat, c(2,1,3)) 
            ptemp <- yhat[curve,,] + nn[curve]*resid
        }
         else {
            yhat <- matrix(temp, nrow=ncurve, byrow=TRUE)
            ptemp <- yhat[curve,] + nn[curve]*resid
        }
    }

    else { #auc
        fit <- survfit0(fit) # add in time 0
        aucfun <- function(x, y, times) {
            # the nuisance is allowing for in-between values  
            # the solution is to expand the list of times, and fill in y for
            #  those insertions.
            alltime <- unique(sort(c(x, times)))
            alltime <- alltime[alltime <= max(times)]
            index <- findInterval(alltime, x, left.open=FALSE)
            delta <- c(diff(alltime), 0)
            i2 <- match(times, alltime)
            if (is.matrix(y)) {
               temp <- apply(y[index,,drop=FALSE], 2, 
                             function(y) cumsum(delta*y))
               temp[i2-1,,drop=FALSE]
            }
            else {
                temp <- cumsum(y[index]*delta)
                temp[i2-1]
            }       
        }                  

        if (inherits(fit, "survfitms")) { # multi-state
            if (is.null(fit$strata)) {
                auc <- aucfun(fit$time, fit$pstate, times)
                yhat <- array(auc, dim=c(1, ntime, nstate))
                ptemp <- yhat[curve,,] + nn[curve]*resid
            } else {
                yhat <- array(0, dim=c(ncurve, ntime, nstate))
                for (i in 1:ncurve) {
                    temp <- fit[i,]
                    yhat[i,,] <- aucfun(temp$time, temp$pstate, times)
                }
                ptemp <- yhat[curve,,] + nn[curve]*resid
            }
         } else { # simple survival
             if (is.null(fit$strata)) {
                 yhat <- aucfun(fit$time, fit$surv, times) # will be a vector
                 ptemp <- yhat[col(resid)] + nn[curve]*resid
             } else {
                 yhat <- matrix(0, ncurve, ntime)
                 for (i in 1:ncurve) {
                     temp <- fit[i]
                     yhat[i,] <- aucfun(temp$time, temp$surv, times)
                 }
             ptemp <- yhat[curve,] + nn[curve]*resid
             } 
         }
    }
        
    if (missing(addNA) && !is.null(id) && 
               (length(id) > length(cluster))) addNA <- FALSE
    if (addNA && length(fit$na.action) >0) {
        if (!is.null(id)) {
            if (length(id) > length(cluster)) {
                # the data set was collapsed over id, we can't do it
                warning("data collapsed on id, addNA option ignored")
            } else {
                ptemp <- pexpand(fit$na.action, ptemp)
                if (data.frame) {
                    omit <- fit$na.action  #observations tossed out
                    n2 <- dim(ptemp)[1]
                    keep <- rep.int(NA, n2)
                    keep[-omit] <- seq(along=id)
                    temp <- id
                    id <- vector(class(id), n2)
                    id[keep] <- temp
                }
            }
        }
        else ptemp <- pexpand(fit$na.action, ptemp)
    }

    if (data.frame) { # return this as a data.frame
        # warning: expand.grid defaults to stringsAsFactors = TRUE
        if (is.null(id)) id <- seq.int(1, dim(ptemp)[1])
        else if (length(id) > length(cluster)) id <- unique(id)
        if (length(fit$states) >0)
             temp <- expand.grid(id=id, times=times, state=fit$states,
                                 stringsAsFactors=FALSE)
        else temp <-  expand.grid(id =id, time=times, stringsAsFactors=FALSE)
        ptemp <- data.frame(temp, pseudo= as.vector(ptemp))
        # if there is no id variable we want to return a column labeled as (Id).
        #  like (Intercept) a convention that made up names not conflict with
        #  standard names.  But anything you do to a dataframe with such a
        #  name will cause it to be replaced with X.Id.  So we have to do this
        #  at the very last.
        if (!is.null(idname)) names(ptemp)[1] <- idname
        else names(ptemp)[1] <- "(Id)"
    }

    ptemp
}        

# this is a copy of the naomit.exclude function
pexpand <- function (omit, x) 
{
    if (is.null(x)) 
        return(x)
    n <- NROW(x)
    keep <- rep.int(NA, n + length(omit))
    keep[-omit] <- 1:n
        
    if (is.matrix(x)) {
        x <- x[keep, , drop = FALSE]
        temp <- rownames(x)
        if (length(temp)) {
            temp[omit] <- names(omit)
            rownames(x) <- temp
        }
    }
    else if (is.array(x) && length(d <- dim(x)) > 2L) {
        x <- x[keep, , , drop = FALSE]
        temp <- (dn <- dimnames(x))[[1L]]
        if (!is.null(temp)) {
            temp[omit] <- names(omit)
            dimnames(x)[[1L]] <- temp
        }
    }
    else {
        x <- x[keep]
        temp <- names(x)
        if (length(temp)) {
            temp[omit] <- names(omit)
            names(x) <- temp
        }
    }
    x
}
