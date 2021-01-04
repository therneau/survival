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
    nstate <- if (inherits(fit, "survfitms")) length(fit$states) else 1L
    nhaz   <- if (is.matrix(fit$cumhaz)) ncol(fit$cumhaz) else 1L

    minus <- ifelse(minus1, 1, 0)  # use n-1 or n to multiply

    # get the estimates
    # retain the dimension of the residuals, which will have subjects first,
    #  timepoints second and (state or cumhaz) last.
    # The residuals have a row for each subject.  If there are multiple groups
    # (strata), the centering constant is a weighted average of the survival
    #   estimates.
    curvewt <- fit$n/sum(fit$n)
    ones <- rep.int(1L, ddr[1])
    if (type== "pstate") {
        temp <- summary(fit, times=times, extend=TRUE)
        if (inherits(fit, "survfitms")) {
            if (ncurve==1) yhat <- array(temp$pstate, dim=c(1,ntime, nstate))
            else {
                yhat <- array(temp$pstate, dim=c(ntime, ncurve, nstate))
                ymean <- apply(yhat, c(1,3), function(x) sum(x* curvewt))
                yhat <- array(ymean, dim=c(1, ncurve, nstate))
            }  
            ptemp <- yhat[ones,,] + (fit$n[curve] -minus)*resid
        } else {
            if (ncurve ==1) yhat <- matrix(temp$surv, nrow=1)
            else {
                ymean <- apply(matrix(temp$surv, ncurve, ntime), 2,
                                 function(x) sum(x*curvewt))
                yhat <- matrix(ymean, nrow=1)
            }
            ptemp <- yhat[ones,] + (fit$n[curve] -minus)*resid
        }
    }
    else if (type == "cumhaz") {
        temp <- summary(fit, times=times, extend=TRUE)$cumhaz
        if (inherits(fit, "survfitms")) {
            if (ncurve==1) yhat <- array(temp$cumhaz, dim=c(1,ntime, nstate))
            else {
                yhat <- array(temp$cumhaz, dim=c(ntime, ncurve, nstate))
                ymean <- apply(yhat, c(1,3), function(x) sum(x* curvewt))
                yhat <- array(ymean, dim=c(1, ncurve, nstate))
            }       
            ptemp <- yhat[ones,,] + (fit$n[curve] -minus)*resid
        } else {
            if (ncurve ==1) yhat <- matrix(temp$cumhaz, nrow=1)
            else {
                ymean <- apply(matrix(temp$cumhaz, ncurve, ntime), 2,
                                 function(x) sum(x*curvewt))
                yhat <- matrix(ymean, nrow=1)
            }
            ptemp <- yhat[ones,] + (fit$n[curve] -minus)*resid
        }
    }
            
    else { #AUC 
        # first compute the AUC at all the timepoints
        fit <- survfit0(fit)  # add the 0 time point
        aucfun <- function(x, y, times) {
            # the nuisance is allowing for in-between values  
            # the solution is to expand the list of times, and fill in y for
            #  those insertions.
            alltime <- sort(c(x, times))
            alltime <- alltime[alltime <= max(times)]
            index <- findInterval(alltime, x, left.open=FALSE)
            delta <- c(diff(alltime), 0)
            i2 <- match(times, alltime)
            if (is.matrix(y)) {
               temp <- apply(y[index,,drop=FALSE], 2, 
                             function(y) cumsum(delta*y))
               temp[i2,,drop=FALSE]
            }
            else {
                temp <- cumsum(y[index]*delta)
                temp[i2]
            }       
        }                  

        if (inherits(fit, "survfitms")) { # multi-state
            if (ncurve==1) {
                auc <- aucfun(fit$time, fit$pstate, times)
                yhat <- array(auc, dim=c(1, ntime, nstate))
            } else {
                yhat <- array(0, dim=c(ntime, nstate, ncurve))
                for (i in 1:ncurve) {
                    temp <- fit[i,]
                    yhat[,,i] <- aucfun(temp$time, temp$pstate, times)
                }
                ymean <- apply(yhat, 1:2, function(x) sum(x*curvewt))
                yhat <- array(ymean, dim=c(1, ntime, nstate))
                }               
            ptemp <- yhat[ones,,] + (fit$n[curve] -minus)*resid

         } else { # simple survival
             if (ncurve==1) {
                 yhat <- aucfun(fit$time, fit$surv, times) # will be a vector
                 ptemp <- yhat[col(resid)] + (fit$n[curve]- minus)*resid
             } else {
                 yhat <- matrix(0, ncurve, ntime)
                 for (i in 1:ncurve) {
                     temp <- fit[i]
                     yhat[i,] <- aucfun(temp$time, temp$surv, times)
                 }
                 yhat <- colMeans(yhat)
             ptemp <- yhat[col(resid)] + (fit$n[curve] -minus)*resid
             } 
         }
    }
        
    if (addNA && length(fit$na.action) >0) {
        ptemp <- pexpand(fit$na.action, ptemp)
        if (!is.null(id) && data.frame) {
            omit <- fit$na.action  #observations tossed out
            n2 <- dim(ptemp)[1]
            keep <- rep.int(NA, n2)
            keep[-omit] <- seq(along=id)
            temp <- id
            id <- vector(class(id), n2)
            id[keep] <- temp
        }
    }
    if (data.frame) { # return this as a data.frame
        # warning: expand.grid defaults to stringsAsFactors = TRUE
        if (is.null(id)) id <- seq.int(1, dim(ptemp)[1])
        if (length(fit$states) >0)
             temp <- expand.grid(id=id, times=times, state=fit$states,
                                 stringsAsFactors=FALSE)
        else temp <-  expand.grid(id =id, time=times, stringsAsFactors=FALSE)
        ptemp <- data.frame(temp, pseudo= as.vector(ptemp))
        # if there is no id variable we want to return a column labeldd as (Id).
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
