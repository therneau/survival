#  Do expected survival based on a Cox model
#  This version relies on the survfit routine to do most of
#   the work.
survexp.cfit <- function(group, ndata, y, method, coxfit, weights) {
    # If it is individual survival, call the predict method
    if (method=='individual') {
        temp <- predict(coxfit, newdata=ndata, type='expect', se=FALSE)
        return(list(surv= exp(-temp)))
    }
 
    # Get the set of survival curves on which I'll base my work
    # There is no id statement allowed yet, so no survexp for time-dependent
    #  covariates
    sfit <- survfit.coxph(coxfit, newdata=ndata, se.fit=FALSE, censor=FALSE)
    
    # rare case: someone called survexp with a single-obs newdata
    #  The average of n curves is just the curve, when n=1
    if (length(group)==1) return(sfit)
    
    # number of curves to create & number of subjects
    ncurve <- max(group)  #group was preset to contain integer group number
    n <- length(group)    # matches nrow(ndata)

    # If the Cox model had strata then the newdata object also had to contain
    #  the strata (needed to fully identify the new subjects), and the 
    #  n survival curves will be "strung out" as a single surv vector in
    #  sfit, along with a strata component saying how many points for each.
    # If the Cox model did not have strata, sfit$surv and sfit$cumhaz will be
    #  matrices with n columns.
    # The output should be a list with components time, n, and surv.
    #   time = vector of unique time points
    #   surv = matrix with 1 column per created curve (often just 1)
    #   n = same shape as surv, containing the number of obs from ndata
    #     that contribute to each row.
    
    newtime <- sort(unique(sfit$time)) # all of the unique times
    ntime <- length(newtime)
    newsurv <- list(time=newtime)

    # Each row of the input data is part of one and only one of the output
    #  curves.  Each column of gmat will contain the weights we need.
    #  Each col sums to 1, and has zeros for those who belong to another curve
    gmat <- matrix(0., nrow=n, ncol=ncurve)
    for (i in 1:ncurve) {
        temp <- weights[group==i]
        gmat[group==i, i] <- temp/sum(temp)
    }

    # If the result is a set of curves with strata rather than a matrix, we
    #  need to index into it, using a code trick taken from summary.survfit
    # Note that is is possible (though odd) for someone to specify a population
    #  of subjects in survexp whose individual members come from different
    #  strata in sfit.  The result curves could have any of the times
    #  that appear in any stratum.  So we create a regular matrix of survivals.
    if (is.null(sfit$strata)) ssurv <- sfit$surv
    else {
        ssurv <- matrix(0., nrow=ntime, ncol=n)
        indx <- rep(1:length(sfit$strata), sfit$strata)
        for (i in 1:n) {
            itemp <- which(indx==i)
            ssurv[,i] <- approx(sfit$time[itemp], sfit$surv[itemp], newtime,
                                 yleft=0, method="constant",
                                 f=0, rule=2)$y
        }
    }

    if (method=="ederer") {
         # This is the most common call. We can work directly
        #  with the returned survival curves, taking weighed averages.
        newsurv$n <- matrix(rep(table(group), each=ntime), nrow=ntime)
        newsurv$surv <- ssurv %*% gmat
    }
    
    else {
        # These are rarely used, so are implemented in S code rather than
        #   C, even though it involves a loop over time points.
        # We need the hazard at each of the new time points, from which
        #  a weighted average at each time point is computed
        #  the Hakulinen also the survival at each time point.
        hazard <- apply(rbind(0, sfit$cumhaz), 2, diff)
        cmat <- matrix(0, ntime, ncurve) # Holds the result

        if (method== "conditional") {
            for (i in 1:ntime) {
                tmat <- ifelse(y >= newtime[i],1,0) * gmat #zero if not at risk
                cmat[i,] <- hazard[i,] %*% tmat / colSums(tmat)
            }
        }
        else { #Hakulinen method
            # Weights in this case are S(newtime) * I(newtime >=y) * gmat
            lsurv <- rbind(1.0, ssurv)  #right continuous time
            for (i in 1:ntime) {
                tmat <- (ifelse(y>=newtime[i],1,0) * lsurv[i,]) * gmat
                cmat[i,] <- hazard[i,] %*% tmat / colSums(tmat)
            }
        }
        newsurv$surv <- exp(-apply(cmat, 2, cumsum))
    }
    newsurv
}

