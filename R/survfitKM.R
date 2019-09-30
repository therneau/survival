# Automatically generated from the noweb directory
survfitKM <- function(x, y, weights=rep(1.0,length(x)), 
                      stype=1, ctype=1,
                      se.fit=TRUE,
                      conf.int= .95,
                      conf.type=c('log',  'log-log',  'plain', 'none', 
                                  'logit', "arcsin"),
                      conf.lower=c('usual', 'peto', 'modified'),
                      start.time, id, cluster, influence=FALSE,
                      type) {
    
    if (!missing(type)) {
        if (!is.character(type)) stop("type argument must be character")
        # older style argument is allowed
        temp <- charmatch(type, c("kaplan-meier", "fleming-harrington", "fh2"))
        if (is.na(temp)) stop("invalid value for 'type'")
        type <- c(1,3,4)[temp]
    }
    else {
        if (!(ctype %in% 1:2)) stop("ctype must be 1 or 2")
        if (!(stype %in% 1:2)) stop("stype must be 1 or 2")
        type <- as.integer(2*stype + ctype  -2)
    }
 
    conf.type <- match.arg(conf.type)
    conf.lower<- match.arg(conf.lower)
    if (is.logical(conf.int)) {
        # A common error is for users to use "conf.int = FALSE"
        #  it's not correct, but allow it
        if (!conf.int) conf.type <- "none"
        conf.int <- .95
    }
      
    # The user can call with cluster, id, both, or neither
    # If only id, treat it as the cluster too
    has.cluster <-  !(missing(cluster) || length(cluster)==0) 
    has.id <-       !(missing(id) || length(id)==0)
    if (has.id) id <- as.factor(id)
    if (has.cluster) {
        if (is.factor(cluster)) {
            clname <- levels(cluster)
            cluster <- as.integer(cluster)
        } else {
            clname  <- sort(unique(cluster))
            cluster <- match(cluster, clname)
        }
        ncluster <- length(clname)
    } else {
        if (has.id) {
            # treat the id as both identifier and clustering
            clname <- levels(id)
            cluster <- as.integer(id)
            ncluster <- length(clname)
        }
        else {
            ncluster <- 0  # has neither
            clname <- NULL
        }
    }
 
    if (is.logical(influence)) {
        # TRUE/FALSE is treated as all or nothing
        if (!influence) influence <- 0L
        else influence <- 3L
    }
    else if (!is.numeric(influence))
        stop("influence argument must be numeric or logical")
    if (!(influence %in% 0:3)) stop("influence argument must be 0, 1, 2, or 3")
    else influence <- as.integer(influence)
 
    if (!se.fit) {
        # if the user asked for no standard error, skip any robust computation
        ncluster <- 0L
        influence <- 0L
    }

    if (!is.Surv(y)) stop("y must be a Surv object")
    if (attr(y, 'type') != 'right' && attr(y, 'type') != 'counting')
            stop("Can only handle right censored or counting data")
    ny <- ncol(y)       # Will be 2 for right censored, 3 for counting
    # The calling routine has used 'strata' on x, so it is a factor with
    #  no unused levels.  But just in case a user called this...
    if (!is.factor(x)) stop("x must be a factor")
    xlev <- levels(x)   # Will supply names for the curves
    x <- as.integer(x)  # keep the integer index

    if (missing(start.time)) time0 <- min(0, y[,ny-1])
    else time0 <- start.time

    if (ny==3 & has.id) position <- survflag(y, id)
    else position <- integer(0)

    if (length(xlev) ==1) {# only one group
        if (ny==2) {
            sort1 <- NULL
            sort2 <- order(y[,1]) 
        }
        else {
            sort2 <- order(y[,2])
            sort1 <- order(y[,1])
        }
        toss <- (y[sort2, ny-1] < time0)
        if (any(toss)) {
            # Some obs were removed by the start.time argument
            sort2 <- sort2[!toss]
            if (ny ==3) {
                index <- match(which(toss), sort1)
                sort1 <- sort1[-index]
            }
        }  
        n.used <- length(sort2)
        if (ncluster > 0)
            cfit <- .Call(Csurvfitkm, y, weights, sort1-1L, sort2-1L, type, 
                                   cluster-1L, ncluster, position, influence)
        else cfit <- .Call(Csurvfitkm, y, weights, sort1-1L, sort2-1L, type,
                                  0L, 0L, position, influence)
    } else {
        # multiple groups
        ngroup <- length(xlev)
        cfit <- vector("list", ngroup)
        n.used <- integer(ngroup)
        if (influence) clusterid <- cfit # empty list of group id values
        for (i in 1:ngroup) {
            keep <- which(x==i & y[,ny-1] >= time0)
            if (length(keep) ==0) next;  # rare case where all are < start.time
            ytemp <- y[keep,]
            n.used[i] <- nrow(ytemp)
            if (ny==2) {
                sort1 <- NULL
                sort2 <- order(ytemp[,1]) 
            }
            else {
                sort2 <- order(ytemp[,2])
                sort1 <- order(ytemp[,1])
            }
     
            # Cluster is a nuisance: every curve might have a different set
            #  We need to relabel them from 1 to "number of unique clusters in this
            #  curve for the C routine
            if (ncluster > 0) {
                c2 <- cluster[keep]
                c.unique <- sort(unique(c2))
                nc <- length(c.unique)
                c2 <- match(c2, c.unique)  # renumber them
                if (influence >0) {
                    clusterid[[i]] <-c.unique
                }
            }
            
            if (ncluster > 0) 
                cfit[[i]] <- .Call(Csurvfitkm, ytemp, weights[keep], sort1 -1L, 
                               sort2 -1L, type,
                               c2 -1L, length(c.unique), position, influence)
            else cfit[[i]] <- .Call(Csurvfitkm, ytemp, weights[keep], sort1 -1L, 
                               sort2 -1L, type,
                               0L, 0L, position, influence)
        }
    }
    # create the survfit object
    if (length(n.used) == 1) {
        rval <- list(n= length(x),
                     time= cfit$time,
                     n.risk = cfit$n[,4],
                     n.event= cfit$n[,5],
                     n.censor=cfit$n[,6],
                     surv = cfit$estimate[,1],
                     std.err = cfit$std[,1],
                     cumhaz  = cfit$estimate[,2],
                     std.chaz = cfit$std[,2])
     } else {
         strata <- sapply(cfit, function(x) nrow(x$n))
         names(strata) <- xlev
         # we need to collapse the curves
         rval <- list(n= as.vector(table(x)),
                      time =   unlist(lapply(cfit, function(x) x$time)),
                      n.risk=  unlist(lapply(cfit, function(x) x$n[,4])),
                      n.event= unlist(lapply(cfit, function(x) x$n[,5])),
                      n.censor=unlist(lapply(cfit, function(x) x$n[,6])),
                      surv =   unlist(lapply(cfit, function(x) x$estimate[,1])),
                      std.err =unlist(lapply(cfit, function(x) x$std[,1])),
                      cumhaz  =unlist(lapply(cfit, function(x) x$estimate[,2])),
                      std.chaz=unlist(lapply(cfit, function(x) x$std[,2])),
                      strata=strata)
          if (ny==3) rval$n.enter <- unlist(lapply(cfit, function(x) x$n[,8]))
    }
        
    if (ny ==3) {
            rval$n.enter <- cfit$n[,8]
            rval$type <- "counting"
    }
    else rval$type <- "right"

    if (se.fit) {
        rval$logse = (ncluster==0 || (type==2 || type==4))  # se(log S) or se(S)
        rval$conf.int   = conf.int
        rval$conf.type= conf.type
        if (conf.lower != "usual") rval$conf.lower = conf.lower

        if (conf.lower == "modified") {
            nstrat = length(n.used)
            events <- rval$n.event >0
            if (nstrat ==1) events[1] <- TRUE
            else           events[1 + cumsum(c(0, rval$strata[-nstrat]))] <- TRUE
            zz <- 1:length(events)
            n.lag <- rep(rval$n.risk[events], diff(c(zz[events], 1+max(zz))))
            #
            # n.lag = the # at risk the last time there was an event (or
            #   the first time of a strata)
            #
        }
        std.low <- switch(conf.lower,
                          'usual' = rval$std.err,
                          'peto' = sqrt((1-rval$surv)/ rval$n.risk),
                          'modified' = rval$std.err * sqrt(n.lag/rval$n.risk))
            
        if (conf.type != "none") {
            ci <- survfit_confint(rval$surv, rval$std.err, logse=rval$logse,
                                  conf.type, conf.int, std.low)
            rval <- c(rval, list(lower=ci$lower, upper=ci$upper))
         }
    } else {
        # for consistency don't return the se if std.err=FALSE
        rval$std.err <- NULL  
        rval$std.chaz <- NULL
    }

    # Add the influence, if requested by the user
    #  remember, if type= 3 or 4, the survival influence has to be constructed.
    if (influence > 0) {
        if (type==1 | type==2) {
            if (influence==1 || influence ==3) {
                if (length(xlev)==1) {
                    rval$influence.surv <- cfit$influence1
                    row.names(rval$influence.surv) <- clname
                } 
                else {
                    temp <- vector("list", ngroup)
                    for (i in 1:ngroup) {
                        temp[[i]] <- cfit[[i]]$influence1
                        row.names(temp[[i]]) <- clname[clusterid[[i]]]
                    }
                    rval$influence.surv <- temp
                }
            }
            if (influence==2 || influence==3) {
                if (length(xlev)==1) {
                    rval$influence.chaz <- cfit$influence2
                    row.names(rval$influence.chaz) <- clname
                }
                else {
                    temp <- vector("list", ngroup)
                    for (i in 1:ngroup) {
                        temp[[i]] <- cfit[[i]]$influence2
                        row.names(temp[[i]]) <- clname[clusterid[[i]]]
                    }
                    rval$influence.chaz <- temp
                }
            }
        }
        else {
            # everything is derived from the influence of the cumulative hazard
            if (length(xlev) ==1) {
                temp <- cfit$influence2
                row.names(temp) <- clname
            } else {
                temp <- vector("list", ngroup)
                for (i in 1:ngroup) {
                    temp[[i]] <- cfit[[i]]$influence2
                    row.names(temp[[i]]) <- clname[clusterid[[i]]]
                }
            }
            
            if (influence==2 || influence ==3)
                rval$influence.chaz <- temp
          
            if (influence==1 || influence==3) {
                # if an obs moves the cumulative hazard up, then it moves S down
                if (length(xlev) ==1) 
                    rval$influence.surv <- -temp * rep(rval$surv, each=nrow(temp))
                else {
                    for (i in 1:ngroup)
                        temp[[i]] <- -temp[[i]] * rep(cfit[[i]]$estimate[,1],
                                                     each=nrow(temp[[i]]))
                    rval$influence.surv <- temp
                }
            }
        }
    }

    if (!missing(start.time)) rval$start.time <- start.time
    rval  
}
