summary.pyears <- function(object, ...) {
    if (!is.null(cl<- object$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
        }

    if (is.null(object$data)) {
        cat (   "Total number of person-years tabulated:", 
             format(sum(object$pyears)),
             "\nTotal number pf person-years off table:",
             format(object$offtable), "\n")
        }
    else {
        cat (   "Total number of person-years tabulated:", 
             format(sum(object$data$pyears)),
             "\nTotal number of person-years off table:",
             format(object$offtable), "")
        }
    if (!is.null(object$summary)) {
        cat("Matches to the chosen rate table:\n  ", 
            object$summary)
        }
    cat("\n")
    }

                               
