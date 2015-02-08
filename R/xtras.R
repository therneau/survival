vcov.coxph<-function (object, ...) {
        rval<-object$var
        dimnames(rval)<-list(names(coef(object)),names(coef(object)))
        rval
    }

vcov.survreg<-function (object, ...) {
        object$var
    }

# The extractAIC methods for coxph and survreg objects are defined
#  in the stats package.  Don't reprise them here.
extractAIC.coxph.penal<- function(fit,scale,k=2,...){
        edf<-sum(fit$df)
        loglik <- fit$loglik[length(fit$loglik)]
        c(edf, -2 * loglik + k * edf)
    }


labels.survreg <- function(object, ...) attr(object,"term.labels")

rep.Surv <- function(x, ...) {
    indx <- rep(1:nrow(x), ...)
    x[indx,]
}

# This function is just like all.vars -- except that it does not recur
#  on the $ sign, it follows both arguments of +, * and - in order to
#  track formulas, all arguments of Surv, and only the first of things 
#  like ns().
# This is used to generate a warning in coxph if the same variable is used
#  on both sides, so perfection is not required.
terms.inner <- function(x) {
    if (class(x) == "formula") c(terms.inner(x[[2]]), terms.inner(x[[3]]))
    else if (class(x)== "call" && 
             (x[[1]] != as.name("$") && x[[1]] != as.name("["))) {
        if (x[[1]] == '+' || x[[1]]== '*' || x[[1]] == '-') {
            # terms in a model equation
            c(terms.inner(x[[2]]), terms.inner(x[[3]]))
        }
        else if (x[[1]] == as.name("Surv"))
                 unlist(lapply(x[-1], terms.inner))
        else terms.inner(x[[2]])
    }
    else(deparse(x))
}

