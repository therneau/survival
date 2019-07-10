# Automatically generated from the noweb directory
# Replaced by concordance.R.  This code is now frozen and will be deprecated
survConcordance.fit <- function(y, x, strata, weight) { 
    # The coxph program may occassionally fail, and this will kill the C
    #  routine below
    if (any(is.na(x)) || any(is.na(y))) return(NULL)   
    btree <- function(n) {
        ranks <- rep(0L, n)  #will be overwritten
        yet.to.do <- 1:n
        depth <- floor(logb(n,2))
        start <- as.integer(2^depth)
        lastrow.length <- 1+n-start
        indx <- seq(1L, by=2L, length= lastrow.length)
        ranks[yet.to.do[indx]] <- start + 0:(length(indx)-1L)
        yet.to.do <- yet.to.do[-indx]

        while (start >1) {
            start <- as.integer(start/2)
            indx <- seq(1L, by=2L, length=start)
            ranks[yet.to.do[indx]] <- start + 0:(start-1L)
            yet.to.do <- yet.to.do[-indx]
        }
        ranks
    }
        
    docount <- function(stime, risk, wts) {
        if (attr(stime, 'type') == 'right') {
            ord <- order(stime[,1], -stime[,2])
            ux <- sort(unique(risk))
            n2 <- length(ux)
            index <- btree(n2)[match(risk[ord], ux)] - 1L
             .Call(Cconcordance1, stime[ord,], 
                   as.double(wts[ord]), 
                   as.integer(index), 
                   as.integer(length(ux)))
        }
        else if (attr(stime, 'type') == "counting") {
            sort.stop <- order(-stime[,2], stime[,3])
            sort.start <- order(-stime[,1])
            ux <- sort(unique(risk))
            n2 <- length(ux)
            index <- btree(n2)[match(risk, ux)] - 1L

            .Call(Cconcordance2, stime, 
                  as.double(wts), 
                  as.integer(index), 
                  as.integer(length(ux)),
                  as.integer(sort.stop-1L), 
                  as.integer(sort.start-1L))
        }
        else stop("Invalid survival type for concordance")
    }
        
    if (missing(weight) || length(weight)==0)
        weight <- rep(1.0, length(x))
    storage.mode(y) <- "double"
    
    if (missing(strata) || length(strata)==0) {
        count <- docount(y, x, weight)
        if (count[1]==0 && count[2]==0) count[5]<-0
        else count[5] <- 2*sqrt(count[5])
        names(count) <- c("concordant", "discordant", "tied.risk", "tied.time",
                          "std(c-d)")
    }
    else {
        strata <- as.factor(strata)
        ustrat <- levels(strata)[table(strata) >0]  #some strata may have 0 obs
        count <- matrix(0., nrow=length(ustrat), ncol=5)
        for (i in 1:length(ustrat)) {
            keep <- which(strata == ustrat[i])
            count[i,] <- docount(y[keep,,drop=F], x[keep], weight[keep])
        }
        
        count[,5] <- 2*sqrt(ifelse(count[,1]+count[,2]==0, 0, count[,5]))
        dimnames(count) <- list(ustrat,  c("concordant", "discordant",
                                           "tied.risk", "tied.time",
                                           "std(c-d)"))
    }
    count
}
