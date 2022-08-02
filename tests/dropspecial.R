# Tests for the drop.special function, which is currently not exported
library(survival)
library(splines)
dfun <- survival:::drop.special

# In a terms structure, the factors attribute is a matrix with row and column
# names.  The predvars and dataClasses attributites, if present, index to the
# row names; as do values of the specials attribute.  The term.labels attribute
# aligns with the column names.

#  For most model formula the row and column names nicely align, but not always.
# [.terms, unfortunately, implicitly assumes that they do align.  This means
# that the result is not always proper.  

# The coxph routine needs to drop terms such as strata, fraily, and cluster
# from the terms structure BEFORE calling model.matrix, otherwise a huge
# matrix could result, from which most of the columns would need to be
# dropped.  (A matched case-control fit will have n/2 strata).

# Things are okay for form1, fail for the others
form0 <- terms(Surv(time, status) ~ age + ns(wt.loss) + strata(sex), 
               specials="strata")
form1 <- terms(Surv(time, status) ~ age + ns(wt.loss) + strata(sex) -1, 
               specials="strata")
form2 <- terms(Surv(time, status) ~ age + offset(ph.ecog) + ns(wt.loss) + 
                   strata(sex), specials= "strata")
form3 <- terms(Surv(time, status) ~ age + strata(sex)/ph.ecog + ns(wt.loss), 
               specials= "strata")
form4 <- terms(~ age + meal.cal*ph.ecog - ph.ecog + strata(sex) + ns(wt.loss), 
               specials = "strata")
form5 <- terms(~ age + ph.ecog - ph.ecog + strata(sex) + ns(wt.loss), 
               specials = "strata")

test0 <- terms(model.frame(form0, data= lung))
test1 <- terms(model.frame(form1, data= lung))
test2 <- terms(model.frame(form2, data= lung))
test3 <- terms(model.frame(form3, data= lung))
test4 <- terms(model.frame(form4, data= lung))
test5 <- terms(model.frame(form5, data= lung))

ccheck <- function(term) {
    aa <- attributes(term)
    cname <- colnames(aa$factors)
    rname <- rownames(aa$factors)
    sindx <- aa$specials$strata
    offset <- aa$offset
    vname <- sapply(aa$variables, deparse)[-1]
    pname <- sapply(aa$predvars,  function(x) deparse(x, nlines=1))[-1]

    # predvars contains an expanded form of ns(wt.loss), so only use the first
    #  8 chars.  The goal is to see if it is the right terms, in the right order
    test <- c(vname= identical(rname, vname),
              labels = identical(cname, aa$term.labels),
              data = identical(rname, names(aa$dataClasses)),
              predvars= identical(substr(pname, 1,8), substr(rname, 1,8)),
              strata = any(grepl("strata", aa$term.labels)))
    
    test
}

# the untangle.specials was the first attempt to get subscripts right
u0 <- untangle.specials(test0, 'strata')
u1 <- untangle.specials(test1, 'strata')
u2 <- untangle.specials(test2, 'strata')
u3 <- untangle.specials(test3, 'strata')
u4 <- untangle.specials(test4, 'strata')
u5 <- untangle.specials(test5, 'strata')

# All is well
ccheck(test0[-u0$terms])  # last element should be FALSE
ccheck(dfun(test0, attr(test0, "specials")$strata))

ccheck(test1[-u1$terms])
ccheck(dfun(test1, attr(test1, "specials")$strata))

# In most of these, the wrong variable gets dropped (serious), the dataClasses
#  attribute has too many elements (not serious).  Predvars may be out of
#  order (serious).
ccheck(test2[-u2$terms])  # correct is T,T,T, T, F
ccheck(dfun(test2, attr(test2, "specials")$strata))

ccheck(test3[-u3$terms]) # correct is all TRUE
ccheck(dfun(test3, attr(test3, "specials")$strata))

ccheck(test4[-u4$terms]) # correct is T, T, T, T, F
ccheck(dfun(test4, attr(test4, "specials")$strata))

ccheck(test5[-u5$terms]) # correct is T, T, T, T, F
ccheck(dfun(test5, attr(test5, "specials")$strata))

