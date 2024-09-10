# $Id $
residuals.coxph.null <-
  function(object, type=c("martingale", "deviance", "score", "schoenfeld"),
           collapse=FALSE, weighted=FALSE,  ...)
    {
    type <- match.arg(type)
    if (type=='martingale' || type=='deviance') NextMethod()
    else stop(gettextf("'%s' residuals are not defined for a null model", type))
    }
