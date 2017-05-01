
format.Surv <- function(x, ...) format(as.character.Surv(x), ...)

# The function to place a Surv object into a data frame
as.data.frame.Surv <- as.data.frame.model.matrix

