# When X is a model matrix, Splus and R have a different format
#   for the "assign" attribute
# For instance 
#           survreg(Surv(time, status) ~ age + sex + factor(ph.ecog), lung)
# R gives the compact form, a vector (0, 1, 2, 3, 3, 3); which can be
#   read as "the first column of the X matrix (intercept) goes with none of
#   the terms', 'the second column goes with term 1', etc.
# Splus gives a list
#      $(Intercept)     1
#      $age             2
#      $sex             3
#      $factor(ph.ecog) 4 5 6  
#
# This function creates the Splus style of output from the R style.  Several
#  of the routines in the package use this, as it is somewhat easier (more
#  transparent) to work with.  
#   
attrassign <- function(object, tt){
        if (!inherits(tt,"terms"))
                stop("need terms object")
        aa<-attr(object,"assign")
        if (is.null(aa))
                stop("argument is not really a model matrix")
        ll<-attr(tt,"term.labels")
        temp <- c("(Intercept)", ll)[aa+1]  #vector of term names
        # Don't put them in alphabetical order, retain the order we inherited
        split(seq(along.with=temp), factor(temp, levels=unique(temp)))
}


# For a multistate model, create the assign attribute for the
#  expanded X matrix (returned from stacker), using cmap and the original assign
# Say for instance that cmap looks like this, x2 a factor with 3 levels of
# a, b, c
#
#     1:2   1:3   2:3
# x1   1     1     7
# x2b  2     4     4
# x2c  3     5     5
# x3   0     6     8
#
# The new assign has an element for each unique variable (8 elements)
#  while the old had 4
# If there is a term argument, return the Splus style, which in this case
# would be a list with elements of
#   x1_1:2 = 1
#   x2_1:2 = 2,3 
#   x2_1:3 = 4,5
#   x3_1:3 = 6
#   x1_2:3 = 7
#   x3_2:3 = 8
# The latter style is used by the cox.zph routine
#
expandassign <- function(assign, cmap, terms) {
    if (is.list(assign)) {
        # the input is the Splus style, convert to R style
        assign <- rep(1:length(assign), sapply(assign, length))
    }
    # multistate coxph never have an intercept, so assign has no zeros
    
    if (!missing(terms)) { 
        # if terms is present, assume that the user wants the expanded (Splus)
        # style result.  
        if (!inherits(terms, "terms"))
            stop("terms argument is not a terms object")
        tname <- attr(terms, "term.labels")
        allname <- outer(tname[assign], colnames(cmap), paste, sep="_")
        # above provides a name for every element of cmap,
        # we only need/want the unique elementsQ

        j <- (cmap>0 & !duplicated(c(cmap)))
        split(1:max(cmap), factor(allname[j], unique(allname[j])))
    } else {
        # idup: what col of cmap had the first appearance of each coefficient?
        idup <- t(apply(cmap, 1, function(x) match(x, c(0,x)))) - 1L
        # remember, idup will have been transposed from cmap, hence the t()
        # for the above example, idup=  0,0,0,0,  0,1,1,1, 2,1,1,2

        first <- match(unique(cmap[cmap>0]), cmap) # the unique coefs
        new <- assign[row(idup)] + (pmax(0, idup-1)* nrow(cmap))[first]

        # the values of new are in the right order, but might not be sequential
        # The intercept gets mapped to 0
        if (any(new ==0)) match(new, unique(new)) -1L  # model had an intercept
        else  match(new, unique(new))
    }
}
                    
    
