summary.coxph <- function(object,  conf.int = 0.95, scale = 1, ...) {
    cox<-object
    beta <- cox$coefficients * scale
    if (is.null(cox$coefficients)) {   # Null model
	return(object)  #The summary method is the same as print in this case
	}
    nabeta <- !(is.na(beta))          #non-missing coefs
    beta2 <- beta[nabeta]
    if(is.null(beta) | is.null(cox$var))
         stop("Input is not valid")
    se <- sqrt(diag(cox$var)) * scale
    if (!is.null(cox$naive.var)) nse <- sqrt(diag(cox$naive.var))

    rval<-list(call=cox$call,fail=cox$fail, na.action=cox$na.action,
                n=cox$n, loglik=cox$loglik)
    if (!is.null(cox$nevent)) rval$nevent <- cox$nevent

    if (is.null(cox$naive.var)) {
        tmp <- cbind(beta, exp(beta), se, beta/se,
                     pchisq((beta/ se)^2, 1, lower.tail=FALSE))
        dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
                                             "se(coef)", "z", "Pr(>|z|)"))
        }
    else {
        tmp <- cbind(beta, exp(beta), nse, se, beta/se,
                     pchisq((beta/ se)^2, 1, lower.tail=FALSE))
        dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
                               "se(coef)", "robust se", "z", "Pr(>|z|)"))
        }
    rval$coefficients <- tmp

    if (conf.int) {
        z <- qnorm((1 + conf.int)/2, 0, 1)
        tmp <- cbind(exp(beta), exp(-beta), exp(beta - z * se),
            exp(beta + z * se))
        dimnames(tmp) <- list(names(beta), c("exp(coef)", "exp(-coef)",
            paste("lower .", round(100 * conf.int, 2), sep = ""),
            paste("upper .", round(100 * conf.int, 2), sep = "")))
        rval$conf.int <- tmp
        }

    df <- length(beta2)
    logtest <- -2 * (cox$loglik[1] - cox$loglik[2])
    rval$logtest <- c(test=logtest,
                      df=df,
                      pvalue= pchisq(logtest, df, lower.tail=FALSE))
    rval$sctest <- c(test=cox$score,
                     df=df,
                     pvalue= pchisq(cox$score, df, lower.tail=FALSE))
    rval$rsq<-c(rsq=1-exp(-logtest/cox$n),
                 maxrsq=1-exp(2*cox$loglik[1]/cox$n))
    if (!is.null(cox$wald.test)) {
        rval$waldtest<-c(test=as.vector(round(cox$wald.test, 2)),
                         df=df,
                         pvalue=pchisq(as.vector(cox$wald.test), df,
                                       lower.tail=FALSE))
    }
    if (!is.null(cox$rscore))
         rval$robscore<-c(test=cox$rscore,
                          df=df,
                          pvalue= pchisq(cox$rscore, df, lower.tail=FALSE))
    rval$used.robust<-!is.null(cox$naive.var)

    if (!is.null(cox$concordance)) {
        # throw away the extra info, in the name of backwards compatability
        rval$concordance <- cox$concordance[6:7]
        names(rval$concordance) <- c("C", "se(C)")
    }
    if (inherits(cox, "coxphms")) {
        rval$cmap <- cox$cmap
        rval$states <- cox$states
    }

    class(rval)    <-"summary.coxph"
    rval
    }
