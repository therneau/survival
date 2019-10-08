# Automatically generated from the noweb directory
parsecovar1 <- function(flist, statedata) {
    if (any(sapply(flist, function(x) !inherits(x, "formula"))))
        stop("an element of the formula list is not a formula")
    if (any(sapply(flist, length) != 3))
        stop("all formulas must have a left and right side")
    
    # split the formulas into a right hand and left hand side
    lhs <- lapply(flist, function(x) x[-3])
    rhs <- lapply(flist, function(x) x[-2])
    
    # take apart the right hand side of each formula
    # the function below is applied to each formula in turn.
    rh2 <- lapply(rhs, function(form) {
        parts <- strsplit(deparse(form, width.cutoff=300, control=NULL), 
                          '/', fixed=TRUE)[[1]]
        if (length(parts) ==1)  { # nothing after a /
            ival <- NULL; common <- FALSE; fixed <- FALSE; clear <- FALSE;
        }
        else{
            # treat the right hand side as though it were a formula
            optterms <- terms(formula(paste("~", parts[2])))
            ff <- rownames(attr(optterms, "factors"))
            index <- match(ff, c("common", "fixed", "init", "clear"))
            if (any(is.na(index)))
                stop("option not recognized in a covariates formula: ",
                     paste(ff[is.na(index)], collapse=", "))
            common <- any(index==1)
            fixed  <- any(index==2)
            clear  <- any(index==3)
            if (any(index==3)) {
                optatt <- attributes(optterms)
                j <- optatt$variables[1 + which(index==3)]
                j[[1]] <- as.name("list")
                ival <- unlist(eval(j, parent.frame()))
            } else ival <- NULL
        }

        # and now the terms before the slash, which is the actual formula
        #  a formula of -1 +1 is recorded as intercept=TRUE, pasting a -1 on
        #  allows us to tell if a 1 was included.  (But don't add a second ~).
        form <- formula(paste("~ -1 +", substring(parts[1], 2, nchar(parts[1]))))
        list(common=common, fixed=fixed, clear=clear, ival=ival, 
           formula = form) 
    })
    # deal with the left hand side of the formula
    # the next routine cuts at '+' signs
    pcut <- function(form) {
        if (length(form)==3) {
            if (form[[1]] == '+') 
                c(pcut(form[[2]]), pcut(form[[3]]))
            else if (form[[1]] == '~') pcut(form[[2]])
            else list(form)
        }
        else list(form)
    }
    lcut <- lapply(lhs, function(x) pcut(x[[2]]))
    env1 <- new.env(parent= parent.frame(2))
    env2 <- new.env(parent= env1)
    if (missing(statedata)) {
        assign("state", function(...) list(stateid= "state", 
                                           values=c(...)), env1)
        assign("state", list(stateid="state"))
    }
    else {
        for (i in statedata) {
            assign(i, eval(list(stateid=i)), env2)
            tfun <- eval(parse(text=paste0("function(...) list(stateid='"
                                           , i, "', values=c(...))")))
            assign(i, tfun, env1)
        }
    }
    lterm <- lapply(lcut, function(x) {
        lapply(x, function(z) {
            if (length(z)==1) {
                temp <- eval(z, envir= env2)
                if (is.list(temp) && names(temp)[[1]] =="stateid") temp
                else temp
            }
            else if (length(z) ==3 && z[[1]]==':')
                list(left=eval(z[[2]], envir=env2), right=eval(z[[3]], envir=env2))
            else stop("invalid term: ", deparse(z))
        })
    })
    list(rhs = rh2, lhs= lterm)
}
parsecovar2 <- function(covar1, statedata, dformula, Terms, transitions,states) {
    if (is.null(statedata))
        statedata <- data.frame(state = states, stringsAsFactors=FALSE)
    else {
        if (is.null(statedata$state)) 
            stop("the statedata data set must contain a variable 'state'")
        indx1 <- match(states, statedata$state, nomatch=0)
        if (any(indx1==0))
            stop("statedata does not contain all the possible states: ", 
                 states[indx1==0])
        statedata <- statedata[indx1,]   # put it in order
    }
    
    # Statedata might have rows for states that are not in the data set,
    #  for instance if the coxph call had used a subset argument.  Any of
    #  those were eliminated above.
    # Likewise, the formula list might have rules for transitions that are
    #  not present.  Don't worry about it at this stage.
    allterm <- attr(Terms, 'term.labels')
    nterm <- length(allterm)

    # create the map and fill it in with the default formula
    nstate <- length(states)
    tmap <- array(0L, dim=c(nterm+1, nstate, nstate))
    dterms <- match(attr(terms.formula(dformula), "term.labels"), allterm)
    dterms <- c(1L, 1L + dterms)  # add the intercept
    k <- seq(along=dterms)
    for (i in 1:nstate) {
        for (j in 1:nstate) {
            tmap[dterms,j,i] <- k    # fill in in column major order
            k <- k + length(k)
        }
    }
    ncoef <- max(tmap)  # number of coefs used so far
    inits <- NULL
    
    # if there is no formula extension, the middle part of the work is skipped
    if (!is.null(covar1)) {
        for (i in 1:length(covar1$rhs)) {  
            rhs <- covar1$rhs[[i]]
            lhs <- covar1$lhs[[i]]  # the two are the same length
            rterm <- terms.formula(rhs$formula)
            rindex <- 1L + match(attr(rterm, "term.labels"), allterm, nomatch=0)
            if (any(rindex== 1L)) stop("dterm mismatch bug 2")
            if (attr(rterm, "intercept")==1) rindex <- c(1L, rindex)
            
            state1 <- state2 <- NULL
            for (x in lhs) {
                # x is one term
                if (is.null(x$left)) stop("term found without a :", x)
                # left of the colon
                if (!is.list(x$left) && length(x$left) ==1 && x$left==0) 
                    temp1 <- 1:nrow(statedata)
                else if (is.numeric(x$left)) {
                    temp1 <- as.integer(x$left)
                    if (any(temp1 != x$left)) stop("non-integer state number")
                    if (any(temp1 <1 | temp1> nstate))
                        stop("numeric state is out of range")
                }
                else if (is.list(x$left) && names(x$left)[1] == "stateid"){
                    if (is.null(x$left$value)) 
                        stop("state variable with no list of values: ",x$left$stateid)
                    else {
                        if (any(k= is.na(match(x$left$stateid, names(statedata)))))
                            stop(x$left$stateid[k], ": state variable not found")
                        zz <- statedata[[x$left$stateid]]
                        if (any(k= is.na(match(x$left$value, zz))))
                            stop(x$left$value[k], ": state value not found")
                        temp1 <- which(zz %in% x$left$value)
                    }
                }
                else {
                    k <- match(x$left, statedata$state)
                    if (any(is.na(k))) stop(x$left[k], ": state not found")
                    temp1 <- which(statedata$state %in% x$left)
                }
                
                # right of colon
                if (!is.list(x$right) && length(x$right) ==1 && x$right ==0) 
                    temp2 <- 1:nrow(statedata)
                else if (is.numeric(x$right)) {
                    temp2 <- as.integer(x$right)
                    if (any(temp2 != x$right)) stop("non-integer state number")
                    if (any(temp2 <1 | temp2> nstate))
                        stop("numeric state is out of range")
                }
                else if (is.list(x$right) && names(x$right)[1] == "stateid") {
                    if (is.null(x$right$value))
                        stop("state variable with no list of values: ",x$right$stateid)
                    else {
                        if (any(k= is.na(match(x$right$stateid, names(statedata)))))
                            stop(x$right$stateid[k], ": state variable not found")
                        zz <- statedata[[x$right$stateid]]
                        if (any(k= is.na(match(x$right$value, zz))))
                            stop(x$right$value[k], ": state value not found")
                        temp2 <- which(zz %in% x$right$value)
                    }
                }
                else {
                    k <- match(x$right, statedata$state)
                    if (any(is.na(k))) stop(x$right[k], ": state not found")
                    temp2 <- which(statedata$state %in% x$right)
                }


                state1 <- c(state1, rep(temp1, length(temp2)))
                state2 <- c(state2, rep(temp2, each=length(temp1)))
            }           
            npair <- length(state1)
            if (rhs$clear) {
                for(k in 1:npair) tmap[-1, state1[k], state2[k]] <- 0
            }
            if (length(rhs$ival)) 
                inits <- c(inits, list(term=rindex, state1=state1, 
                                       state2= state2, init= rhs$ival))
            j <- ncoef + seq_len(length(rindex))
            if (rhs$common) {
                for(k in 1:npair) tmap[rindex, state1[k], state2[k]] <-j
            }
            else {
                for(k in 1:npair){
                    tmap[rindex, state1[k], state2[k]] <-j
                    j <- j + length(rindex)
                }  
            }       
            ncoef <- max(j)
        }    
    }
    i <- match("(censored)", colnames(transitions), nomatch=0)
    if (i==0) t2 <- transitions
    else t2 <- transitions[,-i, drop=FALSE]   # transitions to 'censor' don't count
    indx1 <- match(rownames(t2), states)
    indx2 <- match(colnames(t2), states)
    tmap2 <- matrix(0L, nrow= 1+nterm, ncol= sum(t2>0))

    trow <- row(t2)[t2>0]
    tcol <- col(t2)[t2>0]
    for (i in 1:nrow(tmap2)) {
        for (j in 1:ncol(tmap2))
            tmap2[i,j] <- tmap[i, indx1[trow[j]], indx2[tcol[j]]]
    }

    # relabel as 1-k
    tmap2[1,] <- match(tmap2[1,], unique(c(0L, tmap2[1,]))) -1L
    if (nrow(tmap2) > 1)
        tmap2[-1,] <- match(tmap2[-1,], unique(c(0L, tmap2[-1,]))) -1L
      
    dimnames(tmap2) <- list(c("(Baseline)", allterm),
                                paste(indx1[trow], indx2[tcol], sep=':')) 
    list(tmap = tmap2, inits=inits, mapid= cbind(indx1[trow], indx2[tcol]))
}
parsecovar3 <- function(tmap, Xcol, Xassign) {
    # sometime X will have an intercept, sometimes not, tmap and cmap
    #  always do
    hasintercept <- (Xassign[1] ==0)

    cmap <- matrix(0L, length(Xcol) + !hasintercept, ncol(tmap))
    cmap[1,] <- tmap[1,]

    xcount <- table(factor(Xassign, levels=1:max(Xassign)))
    mult <- 1+ max(xcount)  #used to keep the coefs in the same oder

    ii <- 1
    for (i in 2:nrow(tmap)) {
        k <- seq_len(xcount[i-1])
        for (j in 1:ncol(tmap))
            cmap[ii+k, j] <- if(tmap[i,j]==0) 0 else tmap[i,j]*mult +k

        ii <- ii + max(k)
    }

    # renumber coefs as 1, 2, 3, ...
    cmap[-1,] <- match(cmap[-1,], sort(unique(c(0L, cmap[-1,])))) -1L
    
    colnames(cmap) <- colnames(tmap)
    if (hasintercept) rownames(cmap) <- Xcol
    else rownames(cmap) <- c("(Baseline)", Xcol)

    cmap
}
