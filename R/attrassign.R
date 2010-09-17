# $Id$
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

attrassign<-function (object, ...) UseMethod("attrassign")

attrassign.lm<-function(object, ...){
	attrassign(model.matrix(object), terms(object))}

attrassign.default<-function(object, tt, ...){
        if (!inherits(tt,"terms"))
                stop("need terms object")
        aa<-attr(object,"assign")
        if (is.null(aa))
                stop("argument is not really a model matrix")
        ll<-attr(tt,"term.labels")
        if (attr(tt,"intercept")>0)
                ll<-c("(Intercept)",ll)
        aaa<-factor(aa,labels=ll)
        split(order(aa),aaa)
}

