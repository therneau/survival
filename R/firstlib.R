.onLoad <- function(lib, pkg) {
    ## survfit.print.n=="start" is compatible with previous R
    ##     and with MASS
    if (is.null(getOption("survfit.print.n")))
        options(survfit.print.n="start")
    
    ## survfit.print.mean==TRUE is compatible with previous R/SPLUS
    if (is.null(getOption("survfit.print.mean")))
        options(survfit.print.mean=FALSE)
}

.onUnload <- function(libpath)
    library.dynam.unload("survival", libpath)


is.category <- function(x) inherits(x,"factor") || is.factor(x)



labels.survreg <- function(object, ...) attr(object,"term.labels")

