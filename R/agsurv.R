# Automatically generated from the noweb directory
agsurv <- function(y, x, wt, risk, survtype, vartype) {
    nvar <- ncol(as.matrix(x))
    status <- y[,ncol(y)]
    dtime <- y[,ncol(y) -1]
    death <- (status==1)

    time <- sort(unique(dtime))
    nevent <- as.vector(rowsum(wt*death, dtime))  
    ncens  <- as.vector(rowsum(wt*(!death), dtime))
    wrisk <- wt*risk
    rcumsum <- function(x) rev(cumsum(rev(x))) # sum from last to first
    nrisk <- rcumsum(rowsum(wrisk, dtime))
    irisk <- rcumsum(rowsum(wt, dtime))
    if (ncol(y) ==2) {
        temp2  <- rowsum(wrisk*x, dtime)
        xsum   <- apply(temp2, 2, rcumsum)
        }
    else {
        delta <- min(diff(time))/2
        etime <- c(sort(unique(y[,1])), max(y[,1])+delta)  #unique entry times
        indx  <- approx(etime, 1:length(etime), time, method='constant',
                        rule=2, f=1)$y   
        esum <- rcumsum(rowsum(wrisk, y[,1]))  #not yet entered
        nrisk <- nrisk - c(esum,0)[indx]
        irisk <- irisk - c(rcumsum(rowsum(wt, y[,1])),0)[indx]
        xout   <- apply(rowsum(wrisk*x, y[,1]), 2, rcumsum) #not yet entered
        xin  <- apply(rowsum(wrisk*x, dtime), 2, rcumsum) # dtime or alive
        xsum  <- xin - (rbind(xout,0))[indx,,drop=F]
        }
        
    ndeath <- rowsum(status, dtime)  #unweighted death count
    ntime  <- length(time)        
    if (survtype ==1) {  #Kalbfleisch-Prentice
        indx <- (which(status==1))[order(dtime[status==1])] #deaths
        km <- .C(Cagsurv4,
             as.integer(ndeath),
             as.double(risk[indx]),
             as.double(wt[indx]),
             as.integer(ntime),
             as.double(nrisk),
             inc = double(ntime))
    }

    if (survtype==3 || vartype==3) {  # Efron approx
        xsum2 <- rowsum((wrisk*death) *x, dtime)
        erisk <- rowsum(wrisk*death, dtime)  #risk score sums at each death
        tsum  <- .C(Cagsurv5, 
                    as.integer(length(nevent)),
                    as.integer(nvar),
                    as.integer(ndeath),
                    as.double(nrisk),
                    as.double(erisk),
                    as.double(xsum),
                    as.double(xsum2),
                    sum1 = double(length(nevent)),
                    sum2 = double(length(nevent)),
                    xbar = matrix(0., length(nevent), nvar))
    }
    haz <- switch(survtype,
                     nevent/nrisk,
                     nevent/nrisk,
                     nevent* tsum$sum1)
    varhaz <- switch(vartype,
                     nevent/(nrisk * 
                               ifelse(nevent>=nrisk, nrisk, nrisk-nevent)),
                     nevent/nrisk^2,
                     nevent* tsum$sum2)
    xbar <- switch(vartype,
                   (xsum/nrisk)*haz,
                   (xsum/nrisk)*haz,
                   nevent * tsum$xbar)

    result <- list(n= nrow(y), time=time, n.event=nevent, n.risk=irisk, 
                   n.censor=ncens, hazard=haz, 
                   cumhaz=cumsum(haz), varhaz=varhaz, ndeath=ndeath, 
                   xbar=apply(matrix(xbar, ncol=nvar),2, cumsum))
    if (survtype==1) result$surv <- km$inc
    result
}
