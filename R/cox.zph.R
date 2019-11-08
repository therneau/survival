# Automatically generated from the noweb directory
cox.zph <- function(fit, transform='km', terms=TRUE, singledf =FALSE, 
                    global=TRUE) {
    Call <- match.call()
    if (!inherits(fit, "coxph") && !inherits(fit, "coxme")) 
        stop ("argument must be the result of Cox model fit")
    if (inherits(fit, "coxph.null"))
        stop("there are no score residuals for a Null model")
    if (!is.null(attr(terms(fit), "specials")[["tt"]]))
        stop("function not defined for models with tt() terms")

    cget <- coxph.getdata(fit, y=TRUE, x=TRUE, stratax=TRUE, weights=TRUE)
    y <- cget$y
    ny <- ncol(y)
    event <- (y[,ny] ==1)
    if (length(cget$strata)) 
        istrat <- as.integer(cget$strata) - 1L # number from 0 for C
    else istrat <- rep(0L, nrow(y))

    varnames <- names(fit$coefficients)
    nvar <- length(varnames)
    # if terms==FALSE the singledf argument is moot, but FALSE leads to a
    #  simpler path through the code
    if (!terms) singledf <- FALSE 
    
    eta <- fit$linear.predictors
    if (!terms) {
        # create a fake asgn that has one value per coefficient
        asgn <- as.list(1:nvar)
        names(asgn) <- names(fit$coefficients)
    }
    else if (inherits(fit, "coxme")) {
        asgn <- attrassign(cget$x, terms(fit))
        # allow for a spelling inconsistency in coxme, later fixed
        if (is.null(fit$linear.predictors)) 
            eta <- fit$linear.predictor
    }
    else   asgn <- fit$assign
        
    if (!is.list(asgn)) stop ("unexpected assign component")

    frail <- grepl("frailty(", names(asgn), fixed=TRUE)
    if (any(frail)) {
        dcol <- unlist(asgn[frail])    # remove these columns from X
        cget$x <- cget$x[, -dcol, drop=FALSE]
        asgn <- asgn[!frail]
    }
    nterm <- length(asgn)
    termname <- names(asgn)

    if (any(is.na(fit$coefficients))) {
        # fix up assign so as to ignore missing coefs, this should be rare
        mcoef <- which(is.na(fit$coefficients))
        asgn <- lapply(asgn, function(i) i[!(i %in% mcoef)])
        asgn <- asgn[sapply(asgn, length)>0]  # drop any that were lost
        termname <- names(asgn)
        nterm <- length(asgn)   # asgn will be 1, 2,2,2, 3, etc
    }
    times <- y[,ny-1]
    if (is.character(transform)) {
        tname <- transform
        ttimes <- switch(transform,
                         'identity'= times,
                         'rank'    = rank(times),
                         'log'     = log(times),
                         'km' = {
                             temp <- survfitKM(factor(rep(1L, nrow(y))),
                                               y, se.fit=FALSE)
                             # A nuisance to do left cont KM
                             indx <- findInterval(times, temp$time, left.open=TRUE)
                             1.0 - c(1, temp$surv)[indx+1]
                         },
                         stop("Unrecognized transform"))
            }
        else {
            tname <- deparse(substitute(transform))
            if (length(tname) >1) tname <- 'user'
            ttimes <- transform(times)
            }
        gtime <- ttimes - mean(ttimes[event]) 

        # Now get the U, information, and residuals
        if (ny==2) {
            ord <- order(istrat, y[,1]) -1L
            resid <- .Call(Czph1, gtime, y, cget$x, eta,
                            cget$weights, istrat, fit$method=="efron", ord)
        }
        else {
            ord1 <- order(-istrat, -y[,1]) -1L   # reverse time for zph2
            ord  <- order(-istrat, -y[,2]) -1L
            resid <- .Call(Czph2, gtime, y, cget$x, eta,
                            cget$weights, istrat, fit$method=="efron", 
                            ord1, ord)
        }
    test <- double(nterm+1)
    df   <- rep(1L, nterm+1)
    u0 <- rep(0, nvar)
    if (!is.null(fit$coxlist2)) { # there are penalized terms
        pmat <- matrix(0., 2*nvar, 2*nvar) # second derivative penalty
        pmat[1:nvar, 1:nvar] <- fit$coxlist2$second
        pmat[1:nvar + nvar, 1:nvar + nvar] <- fit$coxlist2$second
        imatr <- resid$imat + pmat
    }
    else imatr <- resid$imat

    for (ii in 1:nterm) {
        jj <- asgn[[ii]]
        kk <- c(1:nvar, jj+nvar)
        imat <- imatr[kk, kk]
        u <- c(u0, resid$u[jj+nvar])
        if (singledf && length(jj) >1) {
            vv <- solve(imat)[-(1:nvar), -(1:nvar)]
            t1 <- sum(fit$coef[jj] * resid$u[jj+nvar])
            test[ii] <- t1^2 * (fit$coef[jj] %*% vv %*% fit$coef[jj])
            df[ii] <- 1
        }
        else {
            test[ii] <- drop(solve(imat,u) %*% u)
            if (is.null(fit$df)) df[ii] <- length(jj)
            else df[ii] <- fit$df[ii]
        }
    }

    #Global test
    if (global) {
        u <- c(u0, resid$u[-(1:nvar)])
        test[nterm+1] <- solve(imatr, u) %*% u
        if (is.null(fit$df))  df[nterm+1]   <- nvar
        else df[nterm+1] <- sum(fit$df)

        tbl <- cbind(test, df, pchisq(test, df, lower.tail=FALSE))
        dimnames(tbl) <- list(c(termname, "GLOBAL"), c("chisq", "df", "p"))
    }
    else {
        tbl <- cbind(test, df, pchisq(test, df, lower.tail=FALSE))[1:nterm,, drop=FALSE]
        dimnames(tbl) <- list(termname, c("chisq", "df", "p"))
    }

    # The x, y, residuals part is sorted by time within strata; this is
    #  what the C routine zph1 and zph2 return
    indx <- if (ny==2) ord +1 else rev(ord) +1  # return to 1 based subscripts
    indx <- indx[event[indx]]                   # only keep the death times
    rval <- list(table=tbl, x=unname(ttimes[indx]), time=unname(y[indx, ny-1]))
    if (length(cget$strata)) rval$strata <- cget$strata[indx]
    # Watch out for a particular edge case: there is a factor, and one of the
    #   strata happens to not use one of its levels.  The element of resid$used will
    #   be zero, but it really should not.
    used <-resid$used
    for (i in asgn) {
        if (length(i) > 1 && any(used[,i] ==0)) 
            used[,i] <- apply(used[,i,drop=FALSE], 1, max)
    }
        
    # Make the weight matrix
    wtmat <- matrix(0, nvar, nvar)
    for (i in 1:nrow(used))
        wtmat <- wtmat + outer(used[i,], used[i,], pmin)
    # with strata*covariate interactions (multi-state models for instance) the
    #  imatr matrix will be block diagonal.  Don't divide these off diagonal zeros
    #  by a wtmat value of zero.
    vmean <- imatr[1:nvar, 1:nvar, drop=FALSE]/ifelse(wtmat==0, 1, wtmat)

    sresid <- resid$schoen
    if (terms && any(sapply(asgn, length) > 1)) { # collase multi-column terms
        temp <- matrix(0, ncol(sresid), nterm)
        for (i in 1:nterm) {
            j <- asgn[[i]]
            if (length(j) ==1) temp[j, i] <- 1
            else temp[j, i] <- fit$coefficients[j]
        }

        sresid <- sresid %*% temp
        vmean <- t(temp) %*% vmean %*% temp
        used <- used[, sapply(asgn, function(x) x[1]), drop=FALSE]
    }

    dimnames(sresid) <- list(signif(rval$time, 4), termname)

    # for each stratum, rescale the Schoenfeld residuals in that stratum
    sgrp <- rep(1:nrow(used), apply(used, 1, max))
    for (i in 1:nrow(used)) {
        k <- which(used[i,] > 0)
        if (length(k) >0)  { # there might be no deaths in the stratum
            j <- which(sgrp==i)
            if (length(k) ==1) sresid[j,k] <- sresid[j,k]/vmean[k,k]
            else sresid[j, k] <- t(solve(vmean[k, k], t(sresid[j, k, drop=FALSE])))
            sresid[j, -k] <- NA
        }
    } 

    # Add in beta-hat.  For a term with multiple columns we are testing zph for
    #  the linear predictor X\beta, which always has a coefficient of 1
    for (i in 1:nterm) {
        j <- asgn[[i]]
        if (length(j) ==1) sresid[,i] <- sresid[,i] + fit$coefficients[j]
        else sresid[,i] <- sresid[,i] +1
    }

    rval$y <- sresid
    rval$var <- solve(vmean)  

    rval$transform <- tname
    rval$call <- Call
    class(rval) <- "cox.zph"
    return(rval)
}

print.cox.zph <- function(x, digits = max(options()$digits - 4, 3),
                          signif.stars=FALSE, ...)  {
    invisible(printCoefmat(x$table, digits=digits, signif.stars=signif.stars, 
                           P.values=TRUE, has.Pvalue=TRUE, ...))
}
"[.cox.zph" <- function(x, ..., drop=FALSE) {
    i <- ..1
    if (!is.null(x$strata)) {
        y2 <- x$y[,i,drop=FALSE]
        ymiss <- apply(is.na(y2), 1, all)
        if (any(ymiss)) {
            # some deaths played no role in these coefficients
            #  due to a strata * covariate interaction, drop unneeded rows
            z<- list(table=x$table[i,,drop=FALSE], x=x$x[!ymiss], 
                     time= x$time[!ymiss], 
                     strata = x$strata[!ymiss],
                     y = y2[!ymiss,,drop=FALSE],
                     var=x$var[i,i, drop=FALSE], 
                     transform=x$transform, call=x$call)
            }
        else z<- list(table=x$table[i,,drop=FALSE], x=x$x, time= x$time, 
                      strata = x$strata,
                      y = y2,  var=x$var[i,i, drop=FALSE], 
                      transform=x$transform, call=x$call)
    }
    else
        z<- list(table=x$table[i,,drop=FALSE], x=x$x, time= x$time, 
                 y = x$y[,i,drop=FALSE],
                 var=x$var[i,i, drop=FALSE],
                 transform=x$transform, call=x$call)
    class(z) <- class(x)
    z
}
