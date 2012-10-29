#
# Create a data set with a particular time-dependent covariate,
#  the default corresponds to a test proposed by Peter O'Brien
#
survobrien <- function(formula, data, subset,
                       na.action, transform) {

    Call <- match.call()

    if (missing(transform)) transform <-function(x) {
        r <- rank(x, na.last="keep")
        temp <- (r -.5)/ length(r[!is.na(r)]) #percentiles
        log(temp/(1-temp)) #logits
    }
    else if (length(transform(1:10)) != 10)
        stop("Transform function must be 1 to 1")

    # create a call to model.frame() that contains the formula (required)
    #  and any other of the relevant optional arguments
    # then evaluate it in the proper frame
    indx <- match(c("formula", "data", "subset", "na.action"),
                  names(Call), nomatch=0) 
    if (indx[1] ==0) stop("A formula argument is required")
    temp <- Call[c(1,indx)]  # only keep the arguments we wanted
    temp[[1]] <- as.name('model.frame')  # change the function called

    special <- c("strata", "cluster", "tt")
    temp$formula <- if(missing(data)) terms(formula, special)
                    else              terms(formula, special, data=data)
    m <- eval(temp, parent.frame())

    if (nrow(m) ==0) stop("No (non-missing) observations")
    n <- nrow(m)
    Terms <- attr(m, 'terms')

    y <- model.extract(m, 'response')
    if (!inherits(y, "Surv")) stop ("Response must be a survival object")
    if (!attr(y, "type") %in% c("right", "counting"))
        stop("Response must be right censored or (start, stop] data")

    cluster <- untangle.specials(Terms, "cluster")
    if (length(cluster$terms) > 0) {
        if (length(cluster$terms) >1) stop ("Can have only 1 cluster term")
        idvar <- m[[cluster$vars]]
        Terms2 <- Terms[-cluster$terms]
    }
    else {
        idvar <- 1:n
        Terms2 <- Terms
    }

    if (length(attr(Terms, "specials")$strata)) {
        stemp <- untangle.specials(Terms2, 'strata', 1)
        if (length(stemp$terms) >0) #beware strata by covariate interactions
            Terms2 <- Terms2[-stemp$terms] #not needed for model.matrix later
        if (length(stemp$vars)==1) strata.keep <- m[[stemp$vars]]
        else strata.keep <- strata(m[,stemp$vars], shortlabel=TRUE)
        }
    else strata.keep <- NULL

    if (any(attr(Terms2, "order") > 1))
        stop("This function cannot deal with iteraction terms")

    # Figure out which are the continuous predictor variables
    myvars <- attr(Terms2, "term.labels")
    factors <- sapply(m[myvars], is.factor)
    protected <- sapply(m[myvars], function(x) inherits(x, "AsIs"))
    keepers <- factors | protected  #variables to be left alone

    if (all(keepers)) stop ("No continuous variables to modify")
    
    if (ncol(y) ==3) {
        # counting process data
        if (is.null(strata.keep)) {
            etime <- sort(unique(y[y[,3]==1, 2])) #unique event times
            indx <- lapply(etime, function(x) which(y[,1]<x & y[,2] >= x))
        }
        else {
            temp <- unique(data.frame(y[,2], strata.keep)[y[,3]==1,])
            etime <- temp[,1]
            indx <- lapply(1:nrow(temp), function(x)
                          which(y[,1] < temp[x,1] & y[,2]>= temp[x,1] &
                                !strata.keep == temp[x,2]))
        }
        
    }
    else {
        # Simple survival data
        if (is.null(strata.keep)) {
            etime <- sort(unique(y[y[,2]==1,1])) #unique event times
            indx <- lapply(etime, function(x) which(y[,1] >=x))
        }
        else {
            temp <- unique(data.frame(y[,1], strata.keep)[y[,2]==1,])
            etime <- temp[,1]
            indx <- lapply(1:nrow(temp), function(x)
                          which(y[,2] >= temp[x,1] & strata.keep == temp[x,2]))
        }
    }

    # The indx list now has an entry for each event time containing the
    # row numbers of those at risk
    indx2 <- unlist(indx)

    #  Create the new survival variables
    nrisk <- unlist(sapply(indx, length))  #number of obs at risk
    if (ncol(y)==3) {
        newdata <- list(y[indx2, 1], y[indx2,2])
        newdata <- c(newdata, 
                     list(1L*(newdata[[2]]==rep(etime, nrisk) &y[indx2,3]==1)))
    }
    else {
        newdata <- list(y[indx2,1])
        newdata <- c(newdata, 
                     list(1L*(newdata[[1]]==rep(etime, nrisk) &y[indx2,2]==1)))
    }
    names(newdata) <- dimnames(y)[[2]]

    # Add any untransformed variables
    if (any(keepers)) {
        temp <- lapply(myvars[keepers], function(x) all.vars(parse(text=x)))
        knames <- unlist(temp)
    }
    else knames <- NULL
    if (length(strata.keep)) {
        knames <- c(knames, unlist(lapply(names(m)[stemp$vars],
                                         function(x) all.vars(parse(text=x)))))
    }
    if (length(knames))         
        newdata <- c(newdata, lapply(data[knames], function(x) x[indx2]))

    # Add the identifier variable
    if (length(cluster$vars) >0) {
        clname <- all.vars(parse(text=names(m)[cluster$vars]))
        newdata <- c(newdata, lapply(data[clname], function(x) x[indx2]))
    }
    else newdata <- c(newdata, list(".id."=idvar[indx2]))

    # Add transformed variables
    tvars <- myvars[!keepers]  
    newx <- lapply(m[tvars],
                   function(z) unlist(lapply(indx, 
                                      function(x) transform(z[x]))))
    data.frame(c(newdata, newx, list(".strata."=rep(1:length(indx),
                                                    sapply(indx, length)))))
}
