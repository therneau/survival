summary.coxph.penal <-  function(object, conf.int = 0.95, scale=1,
                                 terms=FALSE, ...) {
    
    beta <- object$coefficients
    if (length(beta)==0 && length(object$frail)==0)
            stop("Penalized summary function can't be used for a null model")

    if (length(beta) > 0) { #has non-penalized coefs
	nacoef <- !(is.na(beta))          #non-missing coefs
	beta2 <- beta[nacoef]
	if(is.null(beta2) | is.null(object$var))
		stop("Input is not valid")
	se <- sqrt(diag(object$var))
	}
    #
    # Map terms to special print functions, and the list of iteration histories
    #
    pterms <- object$pterms
    nterms <- length(pterms)
    npenal <- sum(pterms>0)
    print.map <- rep(0,nterms)
    if (!is.null(object$printfun)) {
	temp <- unlist(lapply(object$printfun, is.null))  #which ones are missing
	print.map[pterms>0] <- (1:npenal) * (!temp)
	}

    # Tedious, but build up the coef matrix a term at a time
    print1 <- NULL
    pname1 <- NULL
    if (is.null(object$assign2)) alist <- object$assign[-1]
    else alist <- object$assign2

    print2 <- NULL
    for (i in 1:nterms) {
	kk <- alist[[i]]
	if (print.map[i] >0) {
	    j <- print.map[i]	
	    if (pterms[i]==2) 
		 temp <- (object$printfun[[j]])(object$frail, object$fvar, ,
				  object$df[i], object$history[[j]])
	    else temp <- (object$printfun[[j]])(beta[kk], object$var[kk,kk], 
					   object$var2[kk,kk], 
					   object$df[i], object$history[[j]])
	    print1 <- rbind(print1, temp$coef)
	    if (is.matrix(temp$coef)) {
		xx <- dimnames(temp$coef)[[1]]
		if (is.null(xx))
			xx <- rep(names(pterms)[i], nrow(temp$coef))
		else    xx <- paste(names(pterms)[i], xx, sep=', ')
		pname1 <- c(pname1, xx)
		}
	    else  pname1 <- c(pname1, names(pterms)[i])
	    print2 <- c(print2, temp$history)
	    }

	else if (terms && length(kk)>1) {
	    pname1 <- c(pname1, names(pterms)[i])
	    temp <- coxph.wtest(object$var[kk,kk], beta[kk])$test
	    print1 <- rbind(print1, c(NA, NA, NA,
		       temp, object$df[i], pchisq(temp, 1, lower.tail=FALSE)))
	    }
	else {
	    pname1 <- c(pname1, names(beta)[kk])
	    tempe<- (diag(object$var))[kk]
	    temp <- beta[kk]^2/ tempe
	    print1 <- rbind(print1, cbind(beta[kk], sqrt(tempe),
				      sqrt((diag(object$var2))[kk]), 
				   temp, 1, pchisq(temp, 1, lower.tail=FALSE)))
	    }
	}


    dimnames(print1) <- list(pname1, 
                             c("coef","se(coef)", "se2", "Chisq","DF","p"))

    rval <- object[match(c("call", "fail", "na.action", "n", "nevent", "loglik",
                     "iter", "df"), names(object), nomatch=0)]
    rval$coefficients <- print1
    rval$print2 <- print2
                   
    if(conf.int & length(beta) >0 ) {
        z <- qnorm((1 + conf.int)/2, 0, 1)
        beta <- beta * scale
        se <- se * scale
        tmp <- cbind(exp(beta), exp(-beta), exp(beta - z * se),
            exp(beta + z * se))
        dimnames(tmp) <- list(names(beta), 
			      c("exp(coef)", "exp(-coef)",
            paste("lower .", round(100 * conf.int, 2), sep = ""),
            paste("upper .", round(100 * conf.int, 2), sep = "")))
        rval$conf.int <- tmp
        }

    df <- sum(object$df)
    logtest <- -2 * (object$loglik[1] - object$loglik[2])
    rval$logtest <- c(test = logtest, df=df, 
                      pvalue= pchisq(logtest,df, lower.tail=FALSE))
    
    if (!is.null(object$waldtest)) 
        rval$waldtest <- c(test= object$wald.test, df=df,
                       pvalue = pchisq(object$wald.test, df, lower.tail=FALSE))

    if (!is.null(object$concordance)) {
       ctemp <- object$concordance
       # C= (concordant + tied/2)/(concordant + discordant + tied)
       rval$concordance <- ctemp[c("concordance", "std")]
       names(rval$concordance) <- c("C", "se(C)")
    }
 
    class(rval) <- "summary.coxph.penal"
    rval
    }
