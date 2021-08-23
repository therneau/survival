# This routine is loosely based in drop.terms.
#   In a terms structure, the factors attribute is a matrix with row and column
# names.  The predvars and dataClasses attributites, if present, index to the
# row names; as do values of the specials attribute.  The term.labels attribute
# aligns with the column names.
#  For most model formula the row and column names nicely align, but not always.
# [.terms, unfortunately, implicitly assumes that they do align.
#
#  Unlike drop.terms, do not remove offset terms in the process
drop.special <- function(termobj, i, addparen= TRUE) {
    # First step is to rebuild the formula using term.labels and reformulate
    # Map row name to the right column name
    ff <- attr(termobj, "factors")
    index <- match(rownames(ff)[i], colnames(ff))
    if (any(is.null(index))) stop("failure in drop.specials function")
    
    newterms <- attr(termobj, "term.labels")[-index]
    # the above ignores offsets, add them back in
    if (length(attr(termobj, "offset")) > 0)
        newterms <- c(newterms, rownames(ff)[attr(termobj, "offset")])

    rvar <- if (attr(termobj, "response") ==1) termobj[[2L]]
    # Adding () around each term is for a formula containing  + (sex=='male')
    #   A more sophisticated version might look for strings with
    #  no parenthesis and one of ==, < or >
    if (addparen)
        newformula <- reformulate(paste0("(", newterms, ")"), response= rvar,
                              intercept = attr(termobj, "intercept"),
                              env = environment(termobj))
    else  newformula <- reformulate(newterms, response= rvar,
                              intercept = attr(termobj, "intercept"),
                              env = environment(termobj))
    if (length(newformula) == 0L) newformula <- "1"
    
    # addition of an extra specials label causes no harm
    result <- terms(newformula, specials = names(attr(termobj, "specials")))
    
    # now add back the predvars and dataClasses attributes; which do contain
    # the response and offset.
    index2 <- seq.int(nrow(ff))[-i]
    if (!is.null(attr(termobj, "predvars")))
        attr(result, "predvars") <- attr(termobj, "predvars")[c(1, index2 +1)]
    if (!is.null(attr(termobj, "dataClasses")))
        attr(result, "dataClasses") <- attr(termobj, "dataClasses")[index2]

    result
}
