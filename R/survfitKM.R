# Working routine for survival data with a single endpoint
# 
survfitKM <- function(x, y, weights=rep(1.0,length(x)), 
                      stype=1, ctype=1,
                      se.fit=TRUE,
                      conf.int= .95,
                      conf.type=c('log',  'log-log',  'plain', 'none', 
                                  'logit', "arcsin"),
                      conf.lower=c('usual', 'peto', 'modified'),
                      start.time, id, cluster, robust, influence=FALSE,
                      type, entry=FALSE, time0=FALSE) {
    
    # ctype and stype are preferred, 'type' is for backward compatability
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
    if (!is.logical(entry)) stop("entry argument must be TRUE/FALSE")

    if (!is.Surv(y)) stop("y must be a Surv object")
    if (attr(y, 'type') != 'right' && attr(y, 'type') != 'counting')
            stop("Can only handle right censored or counting data")
    ny <- ncol(y)       # Will be 2 for right censored, 3 for counting
    # The calling routine has used 'strata' on x, so it is a factor with
    #  no unused levels.  But just in case a user called this routine...
    if (!is.factor(x)) stop("x must be a factor")
    xlev <- levels(x)   # Will supply names for the curves
    x <- as.integer(x)  # keep the integer index

    if (missing(start.time)) t0 <- min(0, y[,-ny])
    else {
        if (!is.numeric(start.time) || length(start.time) > 1)
            stop("start.time must be a single numeric value")
        t0 <- start.time
    }
    # delete obs if necessary that fall entirely before the staring time
    #  (do this before the robust/cluster/id logic
    #  below, else we could mess up the clname variable
    keep <- y[,ny-1] >= t0
    if (!all(keep)) {
        if (!any(keep)) stop("all observations removed by start.time")
        y <- y[keep,,drop=FALSE]
        if (length(id) >0) id <- id[keep]
        if (length(cluster) >0) cluster <- cluster[keep]
        x <- x[keep]
        weights <- weights[keep]
    }

    # The user can call with cluster, id, robust, or any combination thereof
    # Default for robust: if cluster or any id with > 1 event or 
    #  any weights that are not 0 or 1, then TRUE
    # If only id is present, treat it as the clustering variable
    has.cluster <- !(missing(cluster) || length(cluster)==0) 
    has.id <-      !(missing(id) || length(id)==0)
    has.rwt<-      (!missing(weights) && any(weights != floor(weights)))
    #has.rwt <- FALSE   # we are rethinking this
    has.robust <-  !missing(robust) && !is.null(robust)

    # Residuals will be in the order of integer(id), and we want them to
    #  be in the same order as the data.  So prevent the default sorted levels
    if (has.id) id <- factor(id, unique(id))

    if (missing(robust) || is.null(robust)) {
        if (influence) {
            robust <- TRUE
            if (!(has.cluster || has.id)) {
                cluster <- seq_along(x)
                clname <- cluster
                has.cluster <- TRUE
            }
        }
        else if (has.cluster || has.rwt ||
                 (has.id && anyDuplicated(id[y[,ncol(y)]==1])))
            robust <- TRUE 
        else robust <- FALSE
    }
    if (!is.logical(robust)) stop("robust must be TRUE/FALSE")

    if (has.cluster) {
        if (!robust) {
            warning("cluster specified with robust=FALSE, cluster ignored")
            ncluster <- 0
            clname <- NULL
        }
        else {
            cluster <- factor(cluster, unique(cluster)) # same arg as id above
            clname <- levels(cluster)
            cluster <- as.integer(cluster)
            ncluster <- length(clname)
        }
    } else if (robust) {
        if (has.id) {
            # treat the id as both identifier and clustering
            clname <- levels(id)  # a prior line has ensured it is a factor
            cluster <- as.integer(id)
            ncluster <- length(clname)
        }
        else if (ncol(y)==2 || !has.robust) {
            # create our own clustering
            n <- nrow(y)
            cluster <- 1:n
            ncluster <- n
            clname <- 1:n
        }   
        else stop("id or cluster option required")
    } else ncluster <- 0
 
    if (is.logical(influence)) {
        # TRUE/FALSE is treated as all or nothing
        if (!influence) influence <- 0L
        else influence <- 3L
    }
    else if (!is.numeric(influence))
        stop("influence argument must be numeric or logical")
    if (!(influence %in% 0:3)) stop("influence argument must be 0, 1, 2, or 3")
    else influence <- as.integer(influence)
    if (!robust && influence >0) {
        warning("robust=FALSE implies influence=FALSE")
        influence <- 0L
    }       
  
    if (!se.fit) {
        # if the user asked for no standard error, skip any robust computation
        ncluster <- 0L
        influence <- 0L
    }
    if (!has.id || ncol(y)==2) entry <- FALSE

    if (ny==3 && has.id) position <- survflag(y, id, x)
    else position <- rep.int(3L, nrow(y))  # every observation stands alone

    if (length(xlev) ==1) {# only one group
        n.used <- nrow(y)
        if (ny==2) {
            sort1 <- NULL
            sort2 <- order(y[,1])
        }
        else {
            sort2 <- order(y[,2])
            sort1 <- order(y[,1])
        }
        if (has.id) n.id <- length(unique(id))

        if (ncluster > 0) {
            # cluster is an integer vector, clname the levels
            cfit <- .Call(Csurvfitkm, y, weights, sort1-1L, sort2-1L, type,
                          cluster- 1L, ncluster, position, influence, 0L, entry)
        }
        else cfit <- .Call(Csurvfitkm, y, weights, sort1-1L, sort2-1L, 
                           type, 0L, 0L, position, influence, 0L, entry)
    } else {
        # multiple groups
        ngroup <- length(xlev)
        n.used <- integer(ngroup)
        cfit <- vector("list", ngroup)
        if (influence) clusterid <- cfit # empty list, fill later with group ids
        if (has.id) n.id <- integer(ngroup)

        # The C routine will be called once per curve (values of x)
        # The y, weights, & position values stay the same across calls,
        #  sort1 and sort2 select out the rows that we need
        if (ncluster >0) ctemp <- integer(length(cluster))
        for (i in 1:ngroup) {
            keep <- which(x==i)
            n.used[i] <- length(keep)
            if (length(keep) ==0) next;  # rare case where all are < start.time
            if (has.id) n.id[i] <- length(unique(id[keep]))
            if (ny==2) {
                sort1 <- NULL
                sort2 <- keep[order(y[keep,1])]
            }
            else {
                sort2 <- keep[order(y[keep,2])]
                sort1 <- keep[order(y[keep,1])]
            }
      
            # Cluster is a nuisance: each curve will often have a different set
            #  Say curve1 had id 1,3,5,...99 and curve2 2,4,...,100  
            # We don't want to add up over the 50 extra zeros for curve1, and 
            #  more importantly shouldn't return rows for all thos
            # Use ctemp to give clusters of 0, 1, 2, ...; it can be full length
            #  since the .Call only looks at the sort1/sort2 rows.
            if (ncluster > 0) {
                c.unique <- unique(cluster[keep])
                ctemp[keep] <- match(cluster[keep], c.unique) -1L 
                if (influence >0) 
                    clusterid[[i]] <- c.unique

                cfit[[i]] <- .Call(Csurvfitkm, y, weights, sort1 -1L, 
                               sort2 -1L, type, ctemp, length(c.unique),
                               position, influence, 0L, entry)
            }
            else cfit[[i]] <- .Call(Csurvfitkm, y, weights, sort1 -1L, 
                               sort2 -1L, type, 0L, 0L, 
                               position, influence, 0L, entry)
        }
    }
    # Create the survfit object by 'stacking' the curves one after the
    #  other. Each is likely to have a unique set of times.
    addcounts <- !isTRUE(all.equal(weights, rep(1.0, length(weights)))) 
    if (length(n.used) == 1) {
        rval <- list(n= n.used,
                     time= cfit$time,
                     n.risk = cfit$n[,4],
                     n.event= cfit$n[,5],
                     n.censor=cfit$n[,6],
                     surv = cfit$estimate[,1],
                     std.err = cfit$std.err[,1],
                     cumhaz  = cfit$estimate[,2],
                     std.chaz = cfit$std.err[,2])
        if (entry) rval$n.enter <- cfit$n[,8]
        if (addcounts) {
            if (entry) {
                rval$counts <- cfit$n[,c(1:3, 7), drop=FALSE]
                colnames(rval$counts) <- c("nrisk", "nevent", "ncensor","nenter")
            } else {
                rval$counts <- cfit$n[,1:3, drop=FALSE]
                colnames(rval$counts) <- c("nrisk", "nevent", "ncensor")
            }
        }
     } else {
         strata <- sapply(cfit, function(x) if (is.null(x$n)) 0L else nrow(x$n))
         names(strata) <- xlev
         # we need to collapse the curves
         rval <- list(n= n.used,
                      time =   unlist(lapply(cfit, function(x) x$time)),
                      n.risk=  unlist(lapply(cfit, function(x) x$n[,4])),
                      n.event= unlist(lapply(cfit, function(x) x$n[,5])),
                      n.censor=unlist(lapply(cfit, function(x) x$n[,6])),
                      surv =   unlist(lapply(cfit, function(x) x$estimate[,1])),
                      std.err =unlist(lapply(cfit, function(x) x$std.err[,1])),
                      cumhaz  =unlist(lapply(cfit, function(x) x$estimate[,2])),
                      std.chaz=unlist(lapply(cfit, function(x) x$std.err[,2])),
                      strata=strata)
         if (entry) rval$n.enter <- unlist(lapply(cfit, function(x) x$n[,8]))
         if (addcounts) {
            if (entry) {
                rval$counts <- do.call("rbind", lapply(cfit, function(x)
                                                              x$n[,c(1:3, 7)]))
                colnames(rval$counts) <- c("nrisk", "nevent", "ncensor","nenter")
            } else {
                rval$counts <- do.call("rbind", lapply(cfit,function(x) 
                                                               x$n[,1:3]))
                colnames(rval$counts) <- c("nrisk", "nevent", "ncensor")
            }
        }
     }
    if (has.id) rval$n.id <- n.id   
    if (ny ==3) rval$type <- "counting"
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

    rval$t0 <- t0
    rval  
}
