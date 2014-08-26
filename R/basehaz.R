#
# This function is simply an alias for "survfit".  In the Cox model
#  case users often look for the words "baseline hazard"
#
basehaz <- function(...) survfit(...)

