# Automatically generated from the noweb directory
parsecovar1 <- function(flist, statedata) {
    if (any(sapply(flist, function(x) !inherits(x, "formula"))))
        stop("an element of the formula list is not a formula")
    if (any(sapply(flist, length) != 3))
        stop("all formulas must have a left and right side")
    
    # split the formulas into a right hand and left hand side
    lhs <- lapply(flist, function(x) x[-3])   # keep the ~
    rhs <- lapply(flist, function(x) x[[3]])  # don't keep the ~
    
    rhs <- parse_rightside(rhs)
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
    list(rhs = rhs, lhs= lterm)
}
rightslash <- function(x) {
    if (class(x) != 'call') return(x)
    else {
        if (x[[1]] == as.name('/')) return(list(x[[2]], x[[3]]))
        else if (x[[1]]==as.name('+') || (x[[1]]==as.name('-') && length(x)==3)||
                 x[[1]]==as.name('*') || x[[1]]==as.name(':')  ||
                 x[[1]]==as.name('%in%')) {
                     temp <- rightslash(x[[3]])
                     if (is.list(temp)) {
                         x[[3]] <- temp[[1]]
                         return(list(x, temp[[2]]))
                     } else {
                         temp <- rightslash(x[[2]])
                         if (is.list(temp)) {
                             x[[2]] <- temp[[2]]
                             return(list(temp[[1]], x))
                         } else return(x)
                     }
                 }
        else return(x)
    }
}
parse_rightside <- function(rhs) {
    parts <- lapply(rhs, rightslash)
    new <- lapply(parts, function(opt) {
        tform <- ~ x    # a skeleton, "x" will be replaced
        if (!is.list(opt)) { # no options for this line
            tform[[2]] <- opt
            list(formula = tform, ival = NULL, common = FALSE,
                 shared = FALSE, prop = FALSE)
        }
        else{
            # treat the option list as though it were a formula
            temp <- ~ x
            temp[[2]] <- opt[[2]]
            optterms <- terms(temp)
            ff <- rownames(attr(optterms, "factors"))
            index <- match(ff, c("common", "shared", "prop", "init"))
            if (any(is.na(index)))
                stop("option not recognized in a covariates formula: ",
                     paste(ff[is.na(index)], collapse=", "))
            common <- any(index==1)
            shared  <- any(index==2)
            prop    <- any(index==3)
            if (shared & prop) shared <- FALSE
            if (any(index==4)) {
                optatt <- attributes(optterms)
                j <- optatt$variables[1 + which(index==3)]
                j[[1]] <- as.name("list")
                ival <- unlist(eval(j, parent.frame()))
            } 
            else ival <- NULL
            tform[[2]] <- opt[[1]] 
            list(formula= tform, ival= ival, common= common, shared=shared,
                 prop=prop)
        }
    })
    new
}
termmatch <- function(f1, f2) {
    # look for f1 in f2, each the factors attribute of a terms object
    if (length(f1)==0) return(NULL)   # a formula with only ~1
    irow <- match(rownames(f1), rownames(f2))
    if (any(is.na(irow))) stop ("termmatch failure 1") 
    hashfun <- function(j) sum(ifelse(j==0, 0, 2^(seq(along=j))))
    hash1 <- apply(f1, 2, hashfun)
    hash2 <- apply(f2[irow,,drop=FALSE], 2, hashfun)
    index <- match(hash1, hash2)
    if (any(is.na(index))) stop("termmatch failure 2")
    index
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
    allterm <- attr(Terms, 'factors')
    nterm <- ncol(allterm)

    # create a map for every transition, even ones that are not used.
    # at the end we will thin it out
    # It has an extra first row for intercept (baseline)
    # Fill it in with the default formula
    nstate <- length(states)
    tmap <- array(0, dim=c(nterm+1, nstate, nstate))
    dmap <- array(seq_len(length(tmap)), dim=c(nterm+1, nstate, nstate)) #unique values
    dterm <- termmatch(attr(terms(dformula), "factors"), allterm)
    dterm <- c(1L, 1L+ dterm)  # add intercept
    tmap[dterm,,] <- dmap[dterm,,]
    inits <- NULL

    if (!is.null(covar1)) {
        for (i in 1:length(covar1$rhs)) {  
            rhs <- covar1$rhs[[i]]
            lhs <- covar1$lhs[[i]]  # one rhs and one lhs per formula
          
            state1 <- state2 <- NULL
            for (x in lhs) {
                # x is one term
                if (!is.list(x) || is.null(x$left)) stop("term found without a ':' ", x)
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
                    if (any(is.na(k))) stop(x$left[is.na(k)], ": state not found")
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
            npair <- length(state1)  # number of state:state pairs for this line

            # update tmap for this set of transitions
            # first, what variables are mentioned, and check for errors
            rterm <- terms(rhs$formula)
            rindex <- 1L + termmatch(attr(rterm, "factors"), allterm)

            # the update.formula function is good at identifying changes
            # formulas that start with  "- x" have to be pasted on carefully
            temp <- substring(deparse(rhs$formula, width.cutoff=500), 2)
            if (substring(temp, 1,1) == '-') dummy <- formula(paste("~ .", temp))
            else dummy <- formula(paste("~. +", temp))

            rindex1 <- termmatch(attr(terms(dformula), "factors"), allterm)
            rindex2 <- termmatch(attr(terms(update(dformula, dummy)), "factors"),
                             allterm)
            dropped <- 1L + rindex1[is.na(match(rindex1, rindex2))] # remember the intercept
            if (length(dropped) >0) {
                for (k in 1:npair) tmap[dropped, state1[k], state2[k]] <- 0
            }
           
            # grab initial values
            if (length(rhs$ival)) 
                inits <- c(inits, list(term=rindex, state1=state1, 
                                       state2= state2, init= rhs$ival))
            
            # adding -1 to the front is a trick, to check if there is a "+1" term
            dummy <- ~ -1 + x
            dummy[[2]][[3]] <- rhs$formula
            if (attr(terms(dummy), "intercept") ==1) rindex <- c(1L, rindex)
         
            # an update of "- sex" won't generate anything to add
            if (length(rindex) > 0) {
                if (rhs$common) {
                    j <- dmap[rindex, state1[1], state2[1]] 
                    for(k in 1:npair) tmap[rindex, state1[k], state2[k]] <- j
                }
                else {
                    for (k in 1:npair)
                        tmap[rindex, state1[k], state2[k]] <- dmap[rindex, state1[k], state2[k]]
                }
            }

            # Deal with the prop argument
            if (rhs$prop && npair>1) {
                j <- dmap[1, state1[1], state2[1]]
                for (k in 2:npair) tmap[1, state1[k], state2[k]] <- -j
            }
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

    # Add a map for any hazards that have PH
    phbaseline <- ifelse(tmap2[1,]<0, -tmap2[1,], 0)
    tmap2[1,] <- abs(tmap2[1,])
    utemp <-  unique(c(0L, tmap2[1,]))
    tmap2[1,] <- match(tmap2[1,], utemp) -1L
    phbaseline <- match(phbaseline, utemp) -1L
                       
    if (nrow(tmap2) > 1)
        tmap2[-1,] <- match(tmap2[-1,], unique(c(0L, tmap2[-1,]))) -1L
      
    dimnames(tmap2) <- list(c("(Baseline)", colnames(allterm)),
                                paste(indx1[trow], indx2[tcol], sep=':')) 
    # mapid gives the from,to for each realized state
    list(tmap = tmap2, inits=inits, mapid= cbind(from=indx1[trow], to=indx2[tcol]),
         phbaseline = phbaseline)
}
parsecovar3 <- function(tmap, Xcol, Xassign, phbaseline=NULL) {
    # sometime X will have an intercept, sometimes not; cmap never does
    hasintercept <- (Xassign[1] ==0)

    nph <- sum(phbaseline > 0)
    cmap <- matrix(0L, length(Xcol) + nph - hasintercept, ncol(tmap))
    uterm <- unique(Xassign[Xassign != 0])   # terms that will have coefficients
    
    xcount <- table(factor(Xassign, levels=1:max(Xassign)))
    mult <- 1+ max(xcount)  # temporary scaling

    ii <- 0
    for (i in uterm) {
        k <- seq_len(xcount[i])
        for (j in 1:ncol(tmap)) 
            cmap[ii+k, j] <- if(tmap[i+1,j]==0) 0 else tmap[i+1,j]*mult +k
        ii <- ii + max(k)
    }

    if (nph > 0) {
        k <- seq_len(nph)
        i <- length(Xcol) + k - hasintercept # extra rows in cmap
        j <- which(phbaseline >0)            # coefficients to add
        cmap[cbind(i, j)] <- k + max(cmap)
        
        # I have changed my mind, twice, about a good name
        #newname <- paste0("(", colnames(tmap)[j], ', ',
        #      colnames(tmap)[phbaseline[j]], ")")
        #newname <- paste0("baseline(", j, " vs ", phbaseline[j], ")")
        newname <- paste0("ph(", colnames(tmap)[j], ",", 
                                 colnames(tmap)[phbaseline[j]], ")")
    } else newname <- NULL

    # renumber coefs as 1, 2, 3, ...
    cmap[,] <- match(cmap, sort(unique(c(0L, cmap)))) -1L
    
    colnames(cmap) <- colnames(tmap)
    if (hasintercept) rownames(cmap) <- c(Xcol[-1], newname)
    else rownames(cmap) <- c(Xcol, newname)

    cmap
}
