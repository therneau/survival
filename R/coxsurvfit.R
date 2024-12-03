# Routine to generate survival curves from a fitted coxph model.
# This is called by survfit.coxph, which does a lot of data and argument
#  checking. It is exposed to users, the rms package for instance calls this
#  interface directly.  If you do so, make sure the arguments are valid!
#
# ctype: cumulative hazard type, 1= Nelson (Breslow) 2= F-H (Efron)
# stype: survival type: 1= product limit, 2= exp(cum hazard)
#   stype will almost always be 2, ctype matches the ties option of coxph
# se.fit: compute the standard error
# varmat: the variance-covariance matrix of the coxph model
# cluster: TRUE/FALSE, do a clustered estimatate
# y, x, wt: data fed into the coxph fit (y a Surv object)
# risk = risk score for each of the fitted subjects
# position: 1*(first time interval for a subject) + 2*(last time interval for
#   a subject); used to distinguish actual censoring
# strata: integer, there will be separate curves per stratum
# oldid: id variable for the original data
# y2, x2, risk2, strata2, id2: variables for the newdata argument
#
# See coxph:survival in the methods document.
#
coxsurv.fit <- function(ctype, stype, se.fit, varmat, cluster, 
                         y, x, wt, risk, position, strata, oldid,
                         y2, x2, risk2, strata2, id2, unlist=TRUE) {

    if (missing(strata) || length(strata)==0) strata <- rep(0L, nrow(y))

    if (is.factor(strata)) ustrata <- levels(strata)
    else                   ustrata <- sort(unique(strata))
    nstrata <- length(ustrata)
    survlist <- vector('list', nstrata)
    names(survlist) <- ustrata
    survtype <- if (stype==1) 1 else ctype+1
    vartype <- survtype
    if (is.null(wt)) wt <- rep(1.0, nrow(y))
    if (is.null(strata)) strata <- rep(1L, nrow(y))
    for (i in 1:nstrata) {
        indx <- which(strata== ustrata[i])
        # agsurv does the actual work
        survlist[[i]] <- agsurv(y[indx,,drop=F], x[indx,,drop=F], 
                                wt[indx], risk[indx],
                                survtype, vartype)
        }

    # Expand out the result; x2 has a row for each row in newdata.
    expand <- function(fit, x2, risk2, varmat, se.fit) {
        if (survtype==1) 
            surv <- cumprod(fit$surv) # Kalbfleisch-Prentice
        else surv <- exp(-fit$cumhaz) # Breslow

        if (is.matrix(x2) && nrow(x2) >1) {  #more than 1 row in newdata
            fit$surv <- outer(surv, risk2, '^')
            dimnames(fit$surv) <- list(NULL, row.names(x2))
            if (se.fit) {
                varh <- matrix(0., nrow=length(fit$varhaz), ncol=nrow(x2))
                for (i in 1:nrow(x2)) {
                    dt <- outer(fit$hazard, x2[i,], '*') - fit$xbar
                    dt <- apply(dt, 2, cumsum) # integrate
                    varh[,i] <- (cumsum(fit$varhaz) + rowSums((dt %*% varmat)* dt))*
                        risk2[i]^2
                    }
                fit$std.err <- sqrt(varh)
                }
            fit$cumhaz <- outer(fit$cumhaz, risk2, '*')
            }
        else {
            fit$surv <- surv^risk2
            if (se.fit) {
                dt <-  outer(fit$hazard, c(x2), '*') - fit$xbar
                dt <-  apply(dt, 2, cumsum) # integrate
                varh <- (cumsum(fit$varhaz) + rowSums((dt %*% varmat)* dt)) * 
                    risk2^2
                fit$std.err <- sqrt(varh)
                }
            fit$cumhaz <- fit$cumhaz * risk2
            }
        fit
        }

    if (missing(id2) || is.null(id2)) 
        result <- lapply(survlist, expand, x2, risk2, varmat, se.fit)
    else {
        # The harder case of time-dependent coefficients or covariates.
        # The output object will have one curve per id, stored in the output
        #  as a "stratum", since each curve might have a different number
        #  of rows.
        # Any given subject will usually spend time in multiple strata; they
        #  can even revisit one of them (odd case but possible.
        # Basic algorithm: the survlist object will have one curve per stratum.
        #  * If subject 1 occupies 3 rows of newdata, then pick off the 
        #  (time1, time2) interval of their first row from the strata marked in
        #  that row, then the (time1, time2) interval of their second row,
        #  then the third.
        #  * Then stitch the results together.
        #  * The second component of the variance has to computed after we
        #    stitch the dt matrix together.
        #  * Because risk2 can vary within a subject, final scaling cannot
        #  *  wait until the end.
        #  * Call onecurve below once per unique id.
        # See survival:time dependent coefficients in the methods document.
        #
        onecurve <- function(survlist, x2, y2, strata2,  risk2, se.fit) {
            ntarget <- nrow(x2)  #number of different time intervals
            surv <- vector('list', ntarget)
            n.event <- n.risk <- n.censor <- varh1 <- dt <-  time <- surv
            hazard  <- vector('list', ntarget)
            stemp <- as.integer(strata2)
            for (i in 1:ntarget) {
                if (i==1) toffset <- 0
                else toffset <- toffset + y2[i-1,2]- y2[i,1]

                slist <- survlist[[stemp[i]]]
                indx <- which(slist$time > y2[i,1] & slist$time <= y2[i,2])
                if (length(indx)==0) {
                    # No deaths or censors in user interval.  Possible
                    # user error, but not uncommon at the tail of the curve.
                } else {
                    time[[i]] <- toffset + slist$time[indx]
                
                    hazard[[i]] <- slist$hazard[indx]*risk2[i]
                    if (survtype==1) surv[[i]] <- slist$surv[indx]^risk2[i]
                    
                    n.event[[i]] <- slist$n.event[indx]
                    n.risk[[i]]  <- slist$n.risk[indx]
                    n.censor[[i]]<- slist$n.censor[indx]
                    dt[[i]] <- (outer(slist$hazard[indx], x2[i,], '*') -
                        slist$xbar[indx,,drop=F]) * risk2[i]
                    varh1[[i]] <- slist$varhaz[indx] *risk2[i]^2 
                }
            }

            cumhaz <- cumsum(unlist(hazard))
            if (survtype==1) surv <- cumprod(unlist(surv))  #increments (K-M)
            else surv <- exp(-cumhaz)

            if (se.fit) {
                # paste together dt, integrate, and compute term 2
                dt <- do.call(rbind, dt)
                dt <- apply(dt, 2, cumsum)
                varh2 <- rowSums((dt %*% varmat)* dt)
 
                list(n=as.vector(table(strata)[stemp[1]]),
                       time= unlist(time),
                       n.risk = unlist(n.risk),
                       n.event= unlist(n.event),
                       n.censor= unlist(n.censor),
                       surv = surv,
                       cumhaz= cumhaz,
                       std.err = sqrt(cumsum(unlist(varh1)) + varh2))
            } else list(n=as.vector(table(strata)[stemp[1]]),
                       time= unlist(time),
                       n.risk = unlist(n.risk),
                       n.event= unlist(n.event),
                       n.censor= unlist(n.censor),
                       surv = surv,
                       cumhaz= cumhaz)
        }

        if (all(id2 ==id2[1])) {
            result <- list(onecurve(survlist, x2, y2, strata2, risk2, se.fit))
        }
        else {
            uid <- unique(id2)
            result <- vector('list', length=length(uid))
            for (i in 1:length(uid)) {
                indx <- which(id2==uid[i])
                result[[i]] <- onecurve(survlist, x2[indx,,drop=FALSE], 
                                         y2[indx,,drop=FALSE], 
                                         strata2[indx],  risk2[indx], se.fit)
            }
            names(result) <- uid
        }
    }

    if (unlist) {
        if (length(result)==1) { # the no strata case
            if (se.fit)
                result[[1]][c("n", "time", "n.risk", "n.event", "n.censor",
                          "surv", "cumhaz", "std.err")]
            else result[[1]][c("n", "time", "n.risk", "n.event", "n.censor",
                          "surv", "cumhaz")]
        }
        else {
            temp <-list(n   =    unlist(lapply(result, function(x) x$n),
                                        use.names=FALSE),
                        time=    unlist(lapply(result, function(x) x$time),
                                        use.names=FALSE),
                        n.risk=  unlist(lapply(result, function(x) x$n.risk),
                                        use.names=FALSE),
                        n.event= unlist(lapply(result, function(x) x$n.event),
                                        use.names=FALSE),
                        n.censor=unlist(lapply(result, function(x) x$n.censor),
                                        use.names=FALSE),
                        strata = sapply(result, function(x) length(x$time)))
            names(temp$strata) <- names(result)

            if ((missing(id2) || is.null(id2)) && nrow(x2)>1) {
                 temp$surv <- t(matrix(unlist(lapply(result, 
                                   function(x) t(x$surv)), use.names=FALSE),
                                       nrow= nrow(x2)))
                 dimnames(temp$surv) <- list(NULL, row.names(x2))
                 temp$cumhaz <- t(matrix(unlist(lapply(result, 
                                   function(x) t(x$cumhaz)), use.names=FALSE),
                                       nrow= nrow(x2)))
                 if (se.fit) 
                     temp$std.err <- t(matrix(unlist(lapply(result,
                                    function(x) t(x$std.err)), use.names=FALSE),
                                             nrow= nrow(x2)))
                 }
            else {             
                temp$surv <- unlist(lapply(result, function(x) x$surv),
                                    use.names=FALSE)
                temp$cumhaz <- unlist(lapply(result, function(x) x$cumhaz),
                                    use.names=FALSE)
                if (se.fit) 
                    temp$std.err <- unlist(lapply(result, 
                                   function(x) x$std.err), use.names=FALSE)
                }
            temp
        } 
    }
    else {
        names(result) <- ustrata
        result
    }
}    
