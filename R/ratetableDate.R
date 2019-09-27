#
# Support for older style ratetables: if the type attribute for the dimension
#  is 3 or 4 (a date) and the associated cutpoint is a vector of integers,
#  then the date has a baseline of 1/1/1960.  (Ratetables predate the
#  Date class).
# The newer and simpler form uses a Date vector for the cutpoints.
#
ratetableDate <- function(x) {
    UseMethod("ratetableDate", x)
    }

# This function places a fake "rtabledate" class on the object, for recognition
rtfun <- function(x) {
    y <- as.vector(x)  # as.integer makes sense, but Dates are double
    class(y) <- "rtabledate"
    y
}
    
# Normally used in R
ratetableDate.Date <- function(x) 
    rtfun(x)

# POSIXt includes both POSIXlt and POSIXct
ratetableDate.POSIXt <- function(x)
    rtfun(as.Date(x))

# Normally Splus
#ratetableDate.timeDate <- function(x)
#    rtfun(x - timeDate('1/1/1970'))

# Therneau's old "date" class (will someday wither away?)
ratetableDate.date <- function(x)  rtfun(x - 3653)

# David James's old "chron" class (will someday wither away)
# Support it without using the chron library, which may not be loaded.
ratetableDate.chron <- function(x) {
    origin <- attr(x, "origin")
    x<- as.numeric(x) + as.Date(paste(origin["year"], origin["month"],
                                          origin["day"], sep='/'))
    rtfun(x)
}
ratetableDate.dates <- ratetableDate.chron

# Old ratetables had an integer based on 1/1/1960
ratetableDate.integer <- function(x) 
    rtfun(x - 3653)    # number of days from 1/1/1960 to 1/1/1970

# leave other data types alone
ratetableDate.default <- function(x) x 
