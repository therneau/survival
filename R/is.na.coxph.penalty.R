#  $Id: is.na.coxph.penalty.S 11447 2010-11-12 15:10:18Z therneau $
# The subscript function for coxph.penalty objects
#  without it the "subset" arg of a model statement tosses
#  away all of the attributes
#
"[.coxph.penalty" <- function(x, ..., drop=FALSE) {
    attlist <- attributes(x)
    attributes(x) <- attlist[match(c('dim', 'dimnames', 'levels', 'class'), 
                                   names(attlist), 0)] 
    x <- NextMethod('[')  #let the default method do actual subscripting

    # Tack back on all of the old attributes, except dim and dimnames
    #   which will have been properly modified by the standard [ method,
    #   "levels" which may have dropped some, and "class" which is special
    attributes(x) <- c(attributes(x),
                       attlist[is.na(match(names(attlist),
                                 c("dim", "dimnames", "levels", "class")))])
    # The class will have lost it's first level
    oldClass(x) <- attlist$class
    return(x)
}
			  

is.na.coxph.penalty <- function(x) {
    if (is.matrix(x))
        is.na(c(unclass(x) %*% rep(1,ncol(x))))
    else
        is.na(unclass(x))
}
