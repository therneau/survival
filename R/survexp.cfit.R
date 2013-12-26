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
    #  the strata (needed to fully identify the new subjects), and in this
    #  the n survival curves will be "strung out" as a single surv vector in
    #  the result, along with a strata component saying how many points for each.
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

    # Matrix of 0/1 weights, identify the output curve for each subject.
    gmat <- matrix(0., nrow=n, ncol=ncurve)
    for (i in 1:ncurve) {
        temp <- weights[group==i]
        gmat[group==i, i] <- temp/sum(temp)
    }

    # If the result is a set of curves with strata rather than a matrix, we
    #  need to index into it, using a code trick taken from summary.survfit
    if (!is.null(sfit$strata)) {
        sindex <- matrix(0L, nrow=ntime, ncol=n)
        indx <- rep(1:length(sfit$strata), sfit$strata)
        for (i in 1:n) {
            itemp <- which(indx==i)
            sindex[,i] <- approx(sfit$time[itemp], itemp, newtime,
                                 yleft=0, method="constant",
                                 f=0, rule=2)$y
        }
    }

    if (method=="ederer") {
         # This is the most common call. We can work directly
        #  with the returned survival curves, taking weighed averages.
        newsurv$n <- matrix(rep(table(group), each=ntime), nrow=ntime)
  
        # If there are no strata things are particularly easy
        if(is.null(sfit$strata)) 
            newsurv$surv <- sfit$surv %*% gmat
        else {
            newsurv$surv <- matrix(c(1, sfit$surv)[sindex+1], ntime) %*% gmat
        }
    }
    
    else {
        cmat <- matrix(0., nrow=ntime, ncol=ncurve)
        if (is.null(sfit$strata))
            hazard <- apply(rbind(0, sfit$cumhaz), 2, diff)
        else {
            temp <-  matrix(c(0, sfit$cumhaz)[sindex+1], ntime)
            hazard <-  apply(rbind(0, temp), 2, diff)
        }

        tgt <- outer(newtime, y, function(x,y) x <= y)
        # Weights are now 3 dimensional: subject, output curve, time point
        # The second likely has the smallest dimension
        if (method== "conditional") {
            # modify the weights matrix so that someone drops out once their
            #  "y" value exceeds the time
            for (i in 1:ncurve) {
                wt <- tgt * rep(gmat[,i], each=ntime)
                cmat[,i] <- rowSums(hazard * wt)/rowSums(wt)
            }
        }
        else { #Hakulinen method
            # Weights in this case are S(newtime) * I(newtime >=y) * gmat
            for (i in 1:n) {
                temp <- approx(sfit$time, newtime, sfit$surv[,i],
                               method='constant', rule=2, f=0, yleft=1)$y
                gmat[,i] <- ifelse(y[i]>newtime, 0, gmat[,i]*temp) 
                gmat[,i] <- gmat[,i]/sum(gmat[,i])  #rescale to add to 1
            }
        }
 
       if (is.null(sfit$strata))
            cmat <- hazard %*% gmat
        else
            cmat <- matrix(c(0, hazard)[sindex+1], ntime) %*% gmat
        newsurv$surv <- exp(-apply(cmat, 2, cumsum))
    }
    newsurv
}

