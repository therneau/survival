# Automatically generated from the noweb directory
cmatrix <- function(fit, term, 
                    test =c("global", "trend", "pairwise", "mean"),
                    levels, assign) {
    # Make sure that "fit" is present and isn't missing any parts.
    if (missing(fit)) stop("a fit argument is required")
    Terms <- try(terms(fit), silent=TRUE)

    if (inherits(Terms, "try-error"))
        stop("the fit does not have a terms structure")
    else Terms <- delete.response(Terms)   # y is not needed
    Tatt <- attributes(Terms)
    # a flaw in delete.response: it doesn't subset dataClasses
    Tatt$dataClasses <- Tatt$dataClasses[row.names(Tatt$factors)]
    test <- match.arg(test)

    if (missing(term)) stop("a term argument is required")
    if (is.character(term)) term <- formula(paste("~", term))
    else if (is.numeric(term)) {
        if (all(term == floor(term) & term >0 & term < length(Tatt$term.labels)))
            term <- formula(paste("~", 
                                  paste(Tatt$term.labels[term], collapse='+')))
        else stop("a numeric term must be an integer between 1 and max terms in the fit")
        }
    else if (!inherits(term, "formula"))
        stop("the term must be a formula or integer")
    fterm <- delete.response(terms(term))
    fatt <- attributes(fterm)
    user.name <- fatt$term.labels  # what the user called it
    termname <- all.vars(fatt$variables)
    indx <- match(termname, all.vars(Tatt$variables))
    if (any(is.na(indx))) 
        stop(gettextf("variable %s not found in the formula", termname[is.na(indx)]))
    
    # What kind of term is being tested?  It can be categorical, continuous,
    #  an interaction of only categorical terms, interaction of only continuous
    #  terms, or a mixed interaction.
    # Key is a trick to get "zed" from ns(zed, df= dfvar)
    key <- sapply(Tatt$variables[-1], function(x) all.vars(x)[1])
    parts <- names(Tatt$dataClasses)[match(termname, key)]
    types <- Tatt$dataClasses[parts]
    iscat <- as.integer(types=="factor" | types=="character")
    if (length(iscat)==1) termtype <- iscat
    else  termtype <- 2 + any(iscat) + all(iscat)

    # Were levels specified?  If so we either simply accept them (continuous),
    #  or double check them (categorical)
    if (missing(levels)) {
        temp <- fit$xlevels[match(parts, names(fit$xlevels), nomatch=0)]
        if (length(temp) < length(parts))
            stop("continuous variables require the levels argument")
        levels <- do.call(expand.grid, c(temp, stringsAsFactors=FALSE))
    }
    else {  #user supplied
        if (is.list(levels)) {
            if (is.null(names(levels))) {
                if (length(termname)==1) names(levels)== termname
                else stop("levels list requires named elements")
            }
        }
        if (is.data.frame(levels) || is.list(levels)) {
            index1 <- match(termname, names(levels), nomatch=0)
            # Grab the cols from levels that are needed (we allow it to have
            #  extra, unused columns)
            levels <- as.list(levels[index1])
            # now, levels = the set of ones that the user supplied (which might
            #   be none, if names were wrong)
            if (length(levels) < length(termname)) {
                # add on the ones we don't have, using fit$xlevels as defaults
                temp <- fit$xlevels[parts[index1==0]]
                if (length(temp) > 0) {
                    names(temp) <- termname[index1 ==0]
                    levels <- c(levels, temp)
                }
            } 
            index2 <- match(termname, names(levels), nomatch=0)
            if (any(index2==0)) 
                stop("levels information not found for: ", termname[index2==0])
            levels <- expand.grid(levels[index2], stringsAsFactors=FALSE)
            if (any(duplicated(levels))) stop("levels data frame has duplicates")
        }
        else if (is.matrix(levels)) {
            if (ncol(levels) != length(parts))
                stop("levels matrix has the wrong number of columns")
            if (!is.null(dimnames(levels)[[2]])) {
                index <- match(termname, dimnames(levels)[[2]], nomatch=0)
                if (index==0)
                    stop("matrix column names do no match the variable list")
                else levels <- levels[,index, drop=FALSE]
            } else if (ncol(levels) > 1) 
                stop("multicolumn levels matrix requires column names")
            if (any(duplicated(levels)))
                stop("levels matrix has duplicated rows")
            levels <- data.frame(levels, stringsAsFactors=FALSE)
            names(levels) <- termname
         }
        else if (length(parts) > 1)
            stop("levels should be a data frame or matrix")
        else {
            levels <- data.frame(x=unique(levels), stringsAsFactors=FALSE)
            names(levels) <- termname
        }       
    }

    # check that any categorical levels are legal
    for (i in which(iscat==1)) {
        xlev <- fit$xlevels[[parts[i]]]
        if (is.null(xlev))
            stop(gettextf("xlevels attribute not found for %s", termname[i]))
        temp <- match(levels[[i]], xlev)
        if (any(is.na(temp)))
            stop(gettextf("invalid level for term %s", termname[i]))
    }
    
    rval <- list(levels=levels, termname=termname)
    # Now add the contrast matrix between the levels, if needed
    if (test=="global") {
        if (TRUE) {
        #if (length(parts) ==1) {
            cmat <- diag(nrow(levels))
            cmat[, nrow(cmat)] <- -1   # all equal to the last
            cmat <- cmat[-nrow(cmat),, drop=FALSE]
        }
        else if (termtype== 4) { # anova type
            stop("not yet done 1")
        }
        else stop("not yet done 2")
    }
    else if (test=="pairwise") {
        nlev <- nrow(levels)  # this is the number of groups being compared
        if (nlev < 2) stop("pairwise tests need at least 2 groups")
        npair <- nlev*(nlev-1)/2
        if (npair==1) cmat <- matrix(c(1, -1), nrow=1)
        else {
            cmat <- vector("list", npair)
            k <- 1
            cname <- rep("", npair)
            for (i in 1:(nlev-1)) {
                temp <- double(nlev)
                temp[i] <- 1
                for (j in (i+1):nlev) {
                    temp[j] <- -1
                    cmat[[k]] <- matrix(temp, nrow=1)
                    temp[j] <- 0
                    cname[k] <- paste(i, "vs", j)
                    k <- k+1
                }
            }
            names(cmat) <- cname
        }
    }
    else if (test=="mean") {
        ntest <- nrow(levels)
        cmat <- vector("list", ntest)
        for (k in 1:ntest) {
            temp <- rep(-1/ntest, ntest)
            temp[k] <- (ntest-1)/ntest
            cmat[[k]] <- matrix(temp, nrow=1)
        }
        names(cmat) <- paste(1:ntest, "vs mean")
    }
    else {
        cmat <- vector("list", 2)
        cmat[[1]] <- matrix(1:ntest, 1, ntest)
        cmat[[2]] <- diag(ntest)
        attr(cmat, "nested") <- TRUE
        if (is.null(levels[[1]])) {
            # a continuous variable, and the user didn't give levels for the test
            #  look up the call and use the knots
            tcall <- Tatt$predvars[[indx + 1]]  # skip the 'call' 
            if (tcall[[1]] == as.name("pspline")) {
                bb <- tcall[["Boundary.knots"]]
                levels[[1]] <- seq(bb[1], bb[2], length=ntest)
            }
            else if (tcall[[1]] %in% c("ns", "bs")) {
                bb <- c(tcall[["Boundary.knots"]], tcall[["knots"]])
                levels[[1]] <- sort(bb)
            }
            else stop("don't know how to do a linear contrast for this term")
        }
    }
    # the user can say "age" when the model has "ns(age)", but we need
    #   the more formal label going forward
    rval <- list(levels=levels, termname=parts, cmat=cmat, iscat=iscat)
    class(rval) <- "cmatrix"
    rval
}
gsolve <- function(mat, y, eps=sqrt(.Machine$double.eps)) {
    # solve using a generalized inverse
    # this is very similar to the ginv function of MASS
    temp <- svd(mat, nv=0)
    dpos <- (temp$d > max(temp$d[1]*eps, 0))
    dd <- ifelse(dpos, 1/temp$d, 0)
    # all the parentheses save a tiny bit of time if y is a vector
    if (all(dpos)) x <- drop(temp$u %*% (dd*(t(temp$u) %*% y)))
    else if (!any(dpos)) x <- drop(temp$y %*% (0*y)) # extremely rare
    else x <-drop(temp$u[,dpos] %*%(dd[dpos] * (t(temp$u[,dpos, drop=FALSE]) %*% y)))
    attr(x, "df") <- sum(dpos)
    x
}

qform <- function(var, beta) { # quadratic form b' (V-inverse) b
    temp <- gsolve(var, beta)
    list(test= sum(beta * temp), df=attr(temp, "df"))
}
estfun <- function(cmat, beta, varmat) {
    nabeta <- is.na(beta)
    if (any(nabeta)) {
        k <- which(!nabeta)  #columns to keep
        estimate <- drop(cmat[,k] %*% beta[k])  # vector of predictions
        evar <- cmat[,k] %*% varmat %*% t(cmat[,k, drop=FALSE])
        list(estimate = estimate, var=evar)
    }
    else {
        list(estimate = drop(cmat %*% beta),
             var = cmat %*% varmat %*% t(cmat))
    }
}
             
testfun <- function(cmat, beta, varmat, sigma2) {
    nabeta <- is.na(beta)
    if (any(nabeta)) {
        k <- which(!nabeta)  #columns to keep
        estimate <- drop(cmat[,k] %*% beta[k])  # vector of predictions
        temp <- qform(cmat[,k] %*% varmat %*% t(cmat[,k,drop=FALSE]), estimate)
        rval <- c(chisq=temp$test, df=temp$df)
    }
    else {
       estimate <- drop(cmat %*% beta)
       temp <- qform(cmat %*% varmat %*% t(cmat), estimate)
       rval <- c(chisq=temp$test, df=temp$df)
       }
    if (!is.null(sigma2)) rval <- c(rval, ss= unname(rval[1]) * sigma2)
    rval
}

nafun <- function(cmat, est) {
    used <- apply(cmat, 2, function(x) any(x != 0))
    any(used & is.na(est))
    }
yates <- function(fit, term, population=c("data", "factorial", "sas"),
                  levels, test =c("global", "trend", "pairwise"),
                  predict="linear", options, nsim=200,
                  method=c("direct", "sgtt")) {
    Call <- match.call()
    if (missing(fit)) stop("a fit argument is required")
    Terms <- try(terms(fit), silent=TRUE)
    if (inherits(Terms, "try-error"))
        stop("the fit does not have a terms structure")
    else Terms <- delete.response(Terms)   # y is not needed
    Tatt <- attributes(Terms)
    # a flaw in delete.response: it doesn't subset dataClasses
    Tatt$dataClasses <- Tatt$dataClasses[row.names(Tatt$factors)]
    
    if (inherits(fit, "coxphms")) stop("multi-state coxph not yet supported")
    if (is.list(predict) || is.function(predict)) { 
        # someone supplied their own
        stop("user written prediction functions are not yet supported")
    }
    else {  # call the method
        indx <- match(c("fit", "predict", "options"), names(Call), nomatch=0)
        temp <- Call[c(1, indx)]
        temp[[1]] <- quote(yates_setup)
        mfun <- eval(temp, parent.frame())
    }
    if (is.null(mfun)) predict <- "linear"

   # we will need the original model frame and X matrix
    mframe <- fit$model
    if (is.null(mframe)) mframe <- model.frame(fit)
    Xold <- model.matrix(fit)
    if (is.null(fit$assign)) { # glm models don't save assign
        xassign <- attr(Xold, "assign")
    }
    else xassign <- fit$assign 
    

    nvar <- length(xassign)
    nterm <- length(Tatt$term.names)
    termname <- rownames(Tatt$factors)
    iscat <- sapply(Tatt$dataClasses, 
                    function(x) x %in% c("character", "factor"))
    
    method <- match.arg(casefold(method), c("direct", "sgtt")) #allow SGTT
    if (method=="sgtt" && missing(population)) population <- "sas"

    if (inherits(population, "data.frame")) popframe <- TRUE
    else if (is.character(population)) {
        popframe <- FALSE
        population <- match.arg(tolower(population[1]),
                                c("data", "factorial", "sas",
                                  "empirical", "yates"))
        if (population=="empirical") population <- "data"
        if (population=="yates") population <- "factorial"
    }
    else stop("the population argument must be a data frame or character")
    test <- match.arg(test)
    
    if (popframe || population != "data") weight <- NULL
    else {
        weight <- model.extract(mframe, "weights")
        if (is.null(weight)) {
            id <- model.extract(mframe, "id")
            if (!is.null(id)) { # each id gets the same weight
                count <- c(table(id))
                weight <- 1/count[match(id, names(count))]
            }
        }
    }       

    if (method=="sgtt" && (population !="sas" || predict != "linear"))
        stop("sgtt method only applies if population = sas and predict = linear")

    beta <-  coef(fit, complete=TRUE)
    nabeta <- is.na(beta)  # undetermined coefficients
    vmat <-  vcov(fit, complete=FALSE)
    if (nrow(vmat) > sum(!nabeta)) {
        # a vcov method that does not obey the complete argument
        vmat <- vmat[!nabeta, !nabeta]
    }
    
    # grab the dispersion, needed for the writing an SS in linear models
    if (class(fit)[1] =="lm") sigma <- summary(fit)$sigma
    else sigma <- NULL   # don't compute an SS column
    
    # process the term argument and check its legality
    if (missing(levels)) 
        contr <- cmatrix(fit, term, test, assign= xassign)
    else contr <- cmatrix(fit, term, test, assign= xassign, levels = levels)
    x1data <- as.data.frame(contr$levels)  # labels for the PMM values
    
    # Make the list of X matrices that drive everything: xmatlist
    #  (Over 1/2 the work of the whole routine)
    xmatlist <- yates_xmat(Terms, Tatt, contr, population, mframe, fit,
                                iscat)
 
    # check rows of xmat for estimability
    if (any(is.na(beta)) && (popframe || population != "none")) {
        Xu <- unique(Xold)  # we only need unique rows, saves time to do so
        if (inherits(fit, "coxph")) X.qr <- qr(t(cbind(1.0,Xu)))
        else  X.qr <- qr(t(Xu))   # QR decomposition of the row space
        estimcheck <- function(x, eps= sqrt(.Machine$double.eps)) {
            temp <- abs(qr.resid(X.qr, t(x)))
            # apply(abs(temp), 1, function(x) all(x < eps)) # each row estimable
            all(temp < eps)
        }
        estimable <- sapply(xmatlist, estimcheck)
    } else estimable <- rep(TRUE, length(xmatlist))
    
    # Drop missing coefficients, and use xmatlist to compute the results
    beta <- beta[!nabeta]
    if (predict == "linear" || is.null(mfun)) {
        # population averages of the simple linear predictor
        #temp <- match(contr$termname, colnames(Tatt$factors)) 
        #if (any(is.na(temp)))
        #    stop(gettextf("term '%s' not found in the model", contr$termname[is.na(temp)]))

        meanfun <- if (is.null(weight)) colMeans else function(x) {
            colSums(x*weight)/ sum(weight)}
        Cmat <- t(sapply(xmatlist, meanfun))[,!nabeta]
                  
        # coxph model: the X matrix is built as though an intercept were there (the
        #  baseline hazard plays that role), but then drop it from the coefficients
        #  before computing estimates and tests.  If there was a strata * covariate
        #  interaction there will be many more colums to drop.
        if (inherits(fit, "coxph")) {
            nkeep <- length(fit$means)  # number of non-intercept columns
            col.to.keep <- seq(to=ncol(Cmat), length= nkeep)
            Cmat <- Cmat[,col.to.keep, drop=FALSE]
            offset <- -sum(fit$means[!nabeta] * beta)  # recenter the predictions too
            }
        else offset <- 0
            
        # Get the PMM estimates, but only for estimable ones
        estimate <- cbind(x1data, pmm=NA, std=NA)
        if (any(estimable)) {
            etemp <- estfun(Cmat[estimable,,drop=FALSE], beta, vmat)
            estimate$pmm[estimable] <- etemp$estimate + offset
            estimate$std[estimable] <- sqrt(diag(etemp$var))
        }
            
        # Now do tests on the PMM estimates, one by one
        if (method=="sgtt") {
                # It would be simplest to have the contrasts.arg to be a list of function names.
                # However, model.matrix plays games with the calling sequence, and any function
                #  defined at this level will not be seen.  Instead create a list of contrast
                #  matrices.
                temp <- sapply(fit$contrasts, function(x) (is.character(x) &&
                                           x %in% c("contr.SAS", "contr.treatment")))
                if (!all(temp)) 
                        stop("yates sgtt method can only handle contr.SAS or contr.treatment")
                temp <- vector("list", length(fit$xlevels))
                names(temp) <- names(fit$xlevels)
                for (i in 1:length(fit$xlevels)) {
                    cmat <- diag(length(fit$xlevels[[i]]))
                    dimnames(cmat) <- list(fit$xlevels[[i]], fit$xlevels[[i]])
                    if (i>1 || Tatt$intercept==1) {
                        if (fit$contrasts[[i]] == "contr.treatment")
                            cmat <- cmat[, c(2:ncol(cmat), 1)]
                    }
                    temp[[i]] <- cmat
                }
                sasX <- model.matrix(formula(fit),  data=mframe, xlev=fit$xlevels,
                                      contrasts.arg=temp)
                sas.assign <- attr(sasX, "assign")
                    
                # create the dependency matrix D.  The lm routine is unhappy if it thinks
                #  the right hand and left hand sides are the same, fool it with I().
                # We do this using the entire X matrix even though only categoricals will
                #  eventually be used; if a continuous variable made it NA we need to know.
                D <- coef(lm(sasX ~ I(sasX) -1))
                dimnames(D)[[1]] <- dimnames(D)[[2]] #get rid if the I() names
                zero <- is.na(D[,1])  # zero rows, we'll get rid of these later
                D <- ifelse(is.na(D), 0, D) 
                    
                # make each row orthagonal to rows for other terms that contain it
                #  Containing blocks, if any, will always be below
                # this is easiest to do with the transposed matrix
                # Only do this if both row i and j are for a categorical variable
                if (!all(iscat)) {
                    # iscat marks variables in the model frame as categorical
                    # tcat marks terms as categorical.  For x1 + x2 + x1:x2 iscat has
                    # 2 entries and tcat has 3.
                    tcat <- (colSums(Tatt$factors[!iscat,,drop=FALSE]) == 0)
                }
                else tcat <- rep(TRUE, max(sas.assign)) # all vars are categorical
                   
                B <- t(D)
                dimnames(B)[[2]] <- paste0("L", 1:ncol(B))  # for the user
                if (ncol(Tatt$factors) > 1) {
                    share <- t(Tatt$factors) %*% Tatt$factors
                    nc <- ncol(share)
                    for (i in which(tcat[-nc])) {
                        j <- which(share[i,] > 0 & tcat)
                        k <- j[j>i]  # terms that I need to regress out
                        if (length(k)) {
                            indx1 <- which(sas.assign ==i)
                            indx2 <- which(sas.assign %in% k)
                            B[,indx1] <- resid(lm(B[,indx1] ~ B[,indx2]))
                        }
                    }
                }

                # Cut B back down to the non-missing coefs of the original fit
                Smat <- t(B)[!zero, !zero]
                Sassign <- xassign[!nabeta]
                keep <- match(contr$termname, colnames(Tatt$factors))
                if (length(keep) > 1) { # more than 1 term in the model
                    test <- t(sapply(keep, function(i)
                                   testfun(Smat[Sassign==i,,drop=FALSE], beta, vmat, sigma^2)))
                    rownames(test) <- contr$termname
                }  else {
                    test <- testfun(Smat[Sassign==keep,, drop=FALSE], beta, vmat, sigma^2)
                    test <- matrix(test, nrow=1, 
                                   dimnames=list(contr$termname, names(test)))
                }
        }
        else {
            if (is.list(contr$cmat)) {
                test <- t(sapply(contr$cmat, function(x)
                                 testfun(x %*% Cmat, beta, vmat, sigma^2)))
                natest <- sapply(contr$cmat, nafun, estimate$pmm)
            }
            else {
                test <- testfun(contr$cmat %*% Cmat, beta, vmat, sigma^2)
                test <- matrix(test, nrow=1, 
                               dimnames=list("global", names(test)))
                natest <- nafun(contr$cmat, estimate$pmm)
            }
            if (any(natest)) test[natest,] <- NA
        }
        if (any(estimable)){
        #    Cmat[!estimable,] <- NA
            result <- list(estimate=estimate, test=test, mvar=etemp$var, cmat=Cmat)
            }
        else  result <- list(estimate=estimate, test=test, mvar=NA)
        if (method=="sgtt") result$SAS <- Smat
    }
    else {
        xall <- do.call(rbind, xmatlist)[,!nabeta, drop=FALSE]
        if (inherits(fit, "coxph")) {
            xall <- xall[,-1, drop=FALSE]  # remove the intercept
            eta <- xall %*% beta -sum(fit$means[!nabeta]* beta)
        }
        else eta <- xall %*% beta
        n1 <- nrow(xmatlist[[1]])  # all of them are the same size
        index <- rep(1:length(xmatlist), each = n1)
        if (is.function(mfun)) predfun <- mfun
        else {  # double check the object
            if (!is.list(mfun) || 
                any(is.na(match(c("predict", "summary"), names(mfun)))) ||
                !is.function(mfun$predic) || !is.function(mfun$summary))
                stop("the prediction should be a function, or a list with two functions")
            predfun <- mfun$predict
            sumfun  <- mfun$summary
        }
        pmm <- predfun(eta, xall)
        n2 <- length(eta)
        if (!(is.numeric(pmm)) || !(length(pmm)==n2 || nrow(pmm)==n2))
            stop("prediction function should return a vector or matrix")
        pmm <- rowsum(pmm, index, reorder=FALSE)/n1
        pmm[!estimable,] <- NA

        # get a sample of coefficients, in order to create a variance
        # this is lifted from the mvtnorm code (can't include a non-recommended
        # package in the dependencies)
        tol <- sqrt(.Machine$double.eps)
        if (!isSymmetric(vmat, tol=tol, check.attributes=FALSE))
            stop("variance matrix of the coefficients is not symmetric")
        ev <- eigen(vmat, symmetric=TRUE)
        if (!all(ev$values >= -tol* abs(ev$values[1])))
            warning("variance matrix is numerically not positive definite")
        Rmat <- t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
        bmat <- matrix(rnorm(nsim*ncol(vmat)), nrow=nsim) %*% Rmat
        bmat <- bmat + rep(beta, each=nsim)  # add the mean

        # Now use this matrix of noisy coefficients to get a set of predictions
        # and use those to create a variance matrix
        # Since if Cox we need to recenter each run
        sims <- array(0., dim=c(nsim, nrow(pmm), ncol(pmm)))
        if (inherits(fit, 'coxph')) offset <- bmat %*% fit$means[!nabeta]
        else offset <- rep(0., nsim)
           
        for (i in 1:nsim)
            sims[i,,] <- rowsum(predfun(xall %*% bmat[i,] - offset[i]), index, 
                                reorder=FALSE)/n1
        mvar <- var(sims[,,1])  # this will be used for the tests
        estimate <- cbind(x1data, pmm=unname(pmm[,1]), std= sqrt(diag(mvar)))

        # Now do the tests, on the first column of pmm only
        if (is.list(contr$cmat)) {
            test <- t(sapply(contr$cmat, function(x)
                testfun(x, pmm[,1], mvar[estimable, estimable], NULL)))
            natest <- sapply(contr$cmat, nafun, pmm[,1])
        }
        else {
            test <- testfun(contr$cmat, pmm[,1], mvar[estimable, estimable], NULL)
            test <- matrix(test, nrow=1, 
                           dimnames=list(contr$termname, names(test)))
            natest <- nafun(contr$cmat, pmm[,1])
        }
        if (any(natest)) test[natest,] <- NA
        if (any(estimable))
            result <- list(estimate=estimate,test=test, mvar=mvar)
        else  result <- list(estimate=estimate, test=test, mvar=NA)

        # If there were multiple columns from predfun, compute the matrix of
        #  results and variances 
        if (ncol(pmm) > 1 && any(estimable)){
            pmm <-  apply(sims, 2:3, mean)
            mvar2 <- apply(sims, 2:3, var)
            # Call the summary function, if present
            if (is.list(mfun)) result$summary <- sumfun(pmm, mvar2)
            else {
                result$pmm <- pmm
                result$mvar2 <- mvar2
            }
        }
    }
    result$call <- Call
    class(result) <- "yates"
    result
}
yates_xmat <- function(Terms, Tatt, contr, population, mframe, fit, 
                       iscat, weight) {
    # which variables(s) are in x1 (variables of interest)
    # First a special case of strata(grp):x, which causes strata(grp) not to
    #  appear as a column
    if (any(is.na(match(contr$termname, colnames(Tatt$factors))))) {
        #tis rare
        if (length(contr$termname) > 1) stop("incomplete code 1")
        x1indx <- (contr$termname== rownames(Tatt$factors))
        names(x1indx) <- rownames(Tatt$factors)
        if (!any(x1indx)) stop(gettextf("variable %s not found", contr$termname))
    } else x1indx <- apply(Tatt$factors[,contr$termname,drop=FALSE] >0, 1, any)  
    x2indx <- !x1indx  # adjusters
    if (inherits(population, "data.frame")) pdata <- population  #user data
    else if (population=="data") pdata <- mframe  #easy case
    else if (population=="factorial") 
        pdata <- yates_factorial_pop(mframe, Terms, x2indx, fit$xlevels)
    else if (population=="sas") {
        if (all(iscat[x2indx])) 
            pdata <- yates_factorial_pop(mframe, Terms, x2indx, fit$xlevels)
        else if (!any(iscat[x2indx])) pdata <- mframe # no categoricals
        else { # mixed population
            pdata <- yates_factorial_pop(mframe, Terms, x2indx & iscat, 
                                         fit$xlevels)
            n2 <- nrow(pdata)
            pdata <- pdata[rep(1:nrow(pdata), each=nrow(mframe)), ]
            row.names(pdata) <- 1:nrow(pdata)
            # fill in the continuous
            k <- rep(1:nrow(mframe), n2)
            for (i in which(x2indx & !iscat)) {
                j <- names(x1indx)[i]
                if (is.matrix(mframe[[j]])) 
                    pdata[[j]] <- mframe[[j]][k,, drop=FALSE]
                else pdata[[j]] <- (mframe[[j]])[k]
                attributes(pdata[[j]]) <- attributes(mframe[[j]])
            }
        }
    }
    else stop("unknown population")  # this should have been caught earlier

    # Now create the x1 data set, the unique rows we want to test
    if (is.null(contr$levels)) stop("levels are missing for this contrast")
    x1data <- as.data.frame(contr$levels)  # in case it is a list
    x1name <- names(x1indx)[x1indx]
    for (i in 1:ncol(x1data)) {
        if (is.character(x1data[[i]])) {
            if (is.null(fit$xlevels[[x1name[i]]])) 
                x1data[[i]] <- factor(x1data[[i]])
            else x1data[[i]] <- factor(x1data[[i]], fit$xlevels[[x1name[i]]])
        }
    }

    xmatlist <- vector("list", nrow(x1data))
    if (is.null(attr(pdata, "terms"))) {
        np <- nrow(pdata)
        k <- match(x1name, names(pdata), nomatch=0)
        if (any(k>0)) pdata <- pdata[, -k, drop=FALSE]  # toss out yates var
        for (i in 1:nrow(x1data)) {
            j <- rep(i, np)
            tdata <- cbind(pdata, x1data[j,,drop=FALSE]) # new data set
            xmatlist[[i]] <- model.matrix(Terms, tdata, xlev=fit$xlevels,
                                          contrast.arg= fit$contrasts)
        }
    } else {
        # pdata is a model frame, convert x1data
        # if the name and the class agree we go forward simply
        index <- match(names(x1data), names(pdata), nomatch=0)
            
        if (all(index >0) && 
            identical(lapply(x1data, class), lapply(pdata, class)[index]) &
            identical(sapply(x1data, ncol) , sapply(pdata, ncol)[index]))
                { # everything agrees
            for (i in 1:nrow(x1data)) {
                j <- rep(i, nrow(pdata))
                tdata <- pdata
                tdata[,names(x1data)] <- x1data[j,]
                xmatlist[[i]] <- model.matrix(Terms, tdata,
                                               contrasts.arg= fit$contrasts)
            }
        }
        else {
            # create a subset of the terms structure, for x1 only
            #  for instance the user had age=c(75, 75, 85) and the term was ns(age)
            # then call model.frame to fix it up
            x1term <- Terms[which(x1indx)]
            x1name <- names(x1indx)[x1indx]
            attr(x1term, "dataClasses") <- Tatt$dataClasses[x1name] # R bug
            x1frame <- model.frame(x1term, x1data, xlev=fit$xlevels[x1name])
            for (i in 1:nrow(x1data)) {
                j <- rep(i, nrow(pdata))
                tdata <- pdata
                tdata[,names(x1frame)] <- x1frame[j,]
                xmatlist[[i]] <- model.matrix(Terms, tdata, xlev=fit$xlevels,
                                          contrast.arg= fit$contrasts)
            }
        }
    }      
    
    xmatlist
}
yates_factorial_pop <- function(mframe, terms, x2indx, xlevels) {
    x2name <- names(x2indx)[x2indx]
    dclass <- attr(terms, "dataClasses")[x2name]
    if (!all(dclass %in% c("character", "factor")))
        stop("population=factorial only applies if all the adjusting terms are categorical")
   
    nvar <- length(x2name)
    n2 <- sapply(xlevels[x2name], length)  # number of levels for each
    n <- prod(n2)                          # total number of rows needed
    pdata <- mframe[rep(1, n), -1]  # toss the response
    row.names(pdata) <- NULL        # throw away funny names
    n1 <- 1
    for (i in 1:nvar) {
        j <- rep(rep(1:n2[i], each=n1), length=n)
        xx <- xlevels[[x2name[i]]]
        if (dclass[i] == "factor") 
            pdata[[x2name[i]]] <- factor(j, 1:n2[i], labels= xx)
        else pdata[[x2name[i]]] <- xx[j]
        n1 <- n1 * n2[i]
    }
    attr(pdata, "terms") <- terms
    pdata
}
print.yates <- function(x, digits = max(3, getOption("digits") -2),
                        dig.tst = max(1, min(5, digits-1)),
                        eps=1e-8, ...) {
    temp1 <- x$estimate
    temp1$pmm <- format(temp1$pmm, digits=digits)
    temp1$std <- format(temp1$std, digits=digits)

    # the spaces help separate the two parts of the printout
    temp2 <- cbind(test= paste("    ", rownames(x$test)), 
                   data.frame(x$test), stringsAsFactors=FALSE)
    row.names(temp2) <- NULL

    temp2$Pr <- format.pval(pchisq(temp2$chisq, temp2$df, lower.tail=FALSE),
                            eps=eps, digits=dig.tst)
    temp2$chisq <- format(temp2$chisq, digits= dig.tst)
    temp2$df <- format(temp2$df)
    if (!is.null(temp2$ss)) temp2$ss <- format(temp2$ss, digits=digits)
    
    if (nrow(temp1) > nrow(temp2)) {
        dummy <- temp2[1,]
        dummy[1,] <- ""
        temp2 <- rbind(temp2, dummy[rep(1, nrow(temp1)-nrow(temp2)),])
        }
    if (nrow(temp2) > nrow(temp1)) {
        # get rid of any factors before padding
        for (i in which(sapply(temp1, is.factor))) 
            temp1[[i]] <- as.character(temp1[[i]])
        
        dummy <- temp1[1,]
        dummy[1,] <- ""
        temp1 <- rbind(temp1, dummy[rep(1, nrow(temp2)- nrow(temp1)),])
        }
    print(cbind(temp1, temp2), row.names=FALSE)
    invisible(x)
}
yates_setup <- function(fit, ...)
    UseMethod("yates_setup", fit)

yates_setup.default <- function(fit, type, ...) {
    if (!missing(type) && !(type %in% c("linear", "link")))
        warning(gettextf("no yates_setup method exists for a model of class %s and estimate type %s, linear predictor estimate used by default",
                dQuote(class(fit)[1]), type))
    NULL
}

yates_setup.glm <- function(fit, predict = c("link", "response", "terms", 
                                          "linear"), ...) {
    type <- match.arg(predict)
    if (type == "link" || type== "linear") NULL # same as linear
    else if (type == "response") {
        finv <- family(fit)$linkinv
        function(eta, X) finv(eta)
    }
    else if (type == "terms")
        stop("type terms not yet supported")
}
yates_setup.coxph <- function(fit, predict = c("lp", "risk", "expected",
                                     "terms", "survival", "linear"), 
                              options, ...) {
    type <- match.arg(predict)
    if (type=="lp" || type == "linear") NULL  
    else if (type=="risk") function(eta, X) exp(eta)
    else if (type == "survival") {
        # If there are strata we need to do extra work
        # if there is an interaction we want to suppress a spurious warning
        suppressWarnings(baseline <- survfit(fit, censor=FALSE))
        if (missing(options) || is.null(options$rmean)) 
            rmean <- max(baseline$time)  # max death time
        else rmean <- options$rmean

        if (!is.null(baseline$strata)) 
            stop("stratified models not yet supported")
        cumhaz <- c(0, baseline$cumhaz)
        tt <- c(diff(c(0, pmin(rmean, baseline$time))), 0)
         
        predict <- function(eta, ...) {
            c2 <- outer(exp(drop(eta)), cumhaz)  # matrix of values
            surv <- exp(-c2)
            meansurv <- apply(rep(tt, each=nrow(c2)) * surv, 1, sum)
            cbind(meansurv, surv)
        }
        summary <- function(surv, var) {
            bsurv <- t(surv[,-1])
            std <- t(sqrt(var[,-1]))
            chaz <- -log(bsurv)
            zstat <- -qnorm((1-baseline$conf.int)/2)
            baseline$lower <- exp(-(chaz + zstat*std))
            baseline$upper <- exp(-(chaz - zstat*std))
            baseline$surv <- bsurv
            baseline$std.err  <- std/bsurv
            baselinecumhaz <- chaz
            baseline
        }
        list(predict=predict, summary=summary)
     }
    else stop("type expected is not supported")
}
    
