# Tests for the drop.special function, which is not exported
# 
library(survival)
library(splines)
dfun <- survival:::drop.special

# In a terms structure, the factors attribute is a matrix with row and column
# names.  The predvars and dataClasses attributites, if present, index to the
# row names; as do values of the specials attribute.  The term.labels attribute
# aligns with the column names.

#  For most model formula the row and column names nicely align, but not always.
# [.terms, unfortunately, implicitly assumes that they do align.  This means
# that the result of subscripting is not always proper.  

# The coxph routine needs to drop terms such as strata, fraily, and cluster
# from the terms structure BEFORE calling model.matrix, otherwise a huge
# matrix could result, from which most of the columns would need to be
# dropped.  (A matched case-control fit will have n/2 strata).

# The assumptions are correct for form1, fail for the others
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

# The survival package, in various places, had used a [.terms to remove
#   the strata, via a newTerms <- Terms[-u0$terms] assignment where
#   Terms = terms(model-frame)
# The newer drop.specials function, assigned as 'dfun' above, addresses
#   the issue.
# All is well for test0 and test1: [.terms works as expected
# In the remainder, some combination of a/b/c occurs with [.terms
#   a. the wrong variable gets dropped (serious), 
#   b. the dataClasses attribute has too many elements (not serious).,
#   c.  Predvars may be out of order (serious).
# 
#
ans1 <- c(TRUE, TRUE, TRUE, TRUE, FALSE) # correct answer for all but test3
ans2 <- rep(TRUE, 5) # correct answer for test3
names(ans1) <- names(ans2) <- c("vname", "labels", "data", "predvars", "strata")

all.equal(ccheck(test0[-u0$terms]), ans1) 
all.equal(ccheck(test1[-u1$terms]), ans1)
all.equal(ccheck(test2[-u2$terms]), ans1)
all.equal(ccheck(test3[-u3$terms]), ans2)
all.equal(ccheck(test4[-u4$terms]), ans1)
all.equal(ccheck(test5[-u5$terms]), ans1)

# The dropterms function works in all 6 cases
all.equal(ccheck(dfun(test0, attr(test0, "specials")$strata)), ans1)
all.equal(ccheck(dfun(test1, attr(test1, "specials")$strata)), ans1)
all.equal(ccheck(dfun(test2, attr(test2, "specials")$strata)), ans1)
all.equal(ccheck(dfun(test3, attr(test3, "specials")$strata)), ans2)
all.equal(ccheck(dfun(test4, attr(test4, "specials")$strata)), ans1)
all.equal(ccheck(dfun(test5, attr(test5, "specials")$strata)), ans1)



