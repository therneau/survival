# $Id: format.Surv.S 11273 2009-03-20 15:02:17Z tlumley $
#
format.Surv <- function(x, ...) format(as.character.Surv(x), ...)

# The function to "make something suitable for inclusion in a data frame"
#   was "as.data.frame.x" in versions <5, now it is "data.frameAux.x",
#   so here we have a version specific definition.
# This is needed for both S-plus and R
if (!is.R()) {
    if (version$major >= 5) {
	data.frameAux.Surv <- function(x, ...) data.frameAux.AsIs(x, ...)
	} 
    else as.data.frame.Surv <- as.data.frame.model.matrix
  } else {
    as.data.frame.Surv <- as.data.frame.model.matrix
  }
