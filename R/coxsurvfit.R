# Automatically generated from the noweb directory
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
        
    for (i in 1:nstrata) {
        indx <- which(strata== ustrata[i])
        survlist[[i]] <- agsurv(y[indx,,drop=F], x[indx,,drop=F], 
                                wt[indx], risk[indx],
                                survtype, vartype)
        }
    expand <- function(fit, x2, varmat, se.fit) {
        if (survtype==1) 
            surv <- cumprod(fit$surv)
        else surv <- exp(-fit$cumhaz)

        if (is.matrix(x2) && nrow(x2) >1) {  #more than 1 row in newdata
            fit$surv <- outer(surv, risk2, '^')
            dimnames(fit$surv) <- list(NULL, row.names(x2))
            if (se.fit) {
                varh <- matrix(0., nrow=length(fit$varhaz), ncol=nrow(x2))
                for (i in 1:nrow(x2)) {
                    dt <- outer(fit$cumhaz, x2[i,], '*') - fit$xbar
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
                dt <-  outer(fit$cumhaz, c(x2)) - fit$xbar
                varh <- (cumsum(fit$varhaz) + rowSums((dt %*% varmat)* dt)) * 
                    risk2^2
                fit$std.err <- sqrt(varh)
                }
            fit$cumhaz <- fit$cumhaz * risk2
            }
        fit
        }
    if (missing(id2) || is.null(id2)) 
        result <- lapply(survlist, expand, x2, varmat, se.fit)
    else {
        onecurve <- function(slist, x2, y2, strata2,  risk2, se.fit) {
            ntarget <- nrow(x2)  #number of different time intervals
            surv <- vector('list', ntarget)
            n.event <- n.risk <- n.censor <- varh1 <- varh2 <-  time <- surv
            hazard  <- vector('list', ntarget)
            stemp <- as.integer(strata2)
            timeforward <- 0
            for (i in 1:ntarget) {
                slist <- survlist[[stemp[i]]]
                indx <- which(slist$time > y2[i,1] & slist$time <= y2[i,2])
                if (length(indx)==0) {
                    timeforward <- timeforward + y2[i,2] - y2[i,1]
                    # No deaths or censors in user interval.  Possible
                    # user error, but not uncommon at the tail of the curve.
                }
                else {
                    time[[i]] <- diff(c(y2[i,1], slist$time[indx])) #time increments
                    time[[i]][1] <- time[[i]][1] + timeforward
                    timeforward <- y2[i,2] - max(slist$time[indx])
                
                    hazard[[i]] <- slist$hazard[indx]*risk2[i]
                    if (survtype==1) surv[[i]] <- slist$surv[indx]^risk2[i]
                    
                    n.event[[i]] <- slist$n.event[indx]
                    n.risk[[i]]  <- slist$n.risk[indx]
                    n.censor[[i]]<- slist$n.censor[indx]
                    dt <-  outer(slist$cumhaz[indx], x2[i,]) - slist$xbar[indx,,drop=F]
                    varh1[[i]] <- slist$varhaz[indx] *risk2[i]^2
                    varh2[[i]] <- rowSums((dt %*% varmat)* dt) * risk2[i]^2
                }
            }

            cumhaz <- cumsum(unlist(hazard))
            if (survtype==1) surv <- cumprod(unlist(surv))  #increments (K-M)
            else surv <- exp(-cumhaz)

            if (se.fit) 
                list(n=as.vector(table(strata)[stemp[1]]),
                       time=cumsum(unlist(time)),
                       n.risk = unlist(n.risk),
                       n.event= unlist(n.event),
                       n.censor= unlist(n.censor),
                       surv = surv,
                       cumhaz= cumhaz,
                       std.err = sqrt(cumsum(unlist(varh1)) + unlist(varh2)))
            else list(n=as.vector(table(strata)[stemp[1]]),
                       time=cumsum(unlist(time)),
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
