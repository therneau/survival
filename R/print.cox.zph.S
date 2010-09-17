# $Id: print.cox.zph.S 11166 2008-11-24 22:10:34Z therneau $
print.cox.zph <- function(x, digits = max(options()$digits - 4, 3),...)
    invisible(print(x$table, digits=digits))
