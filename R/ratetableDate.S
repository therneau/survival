#
# survexp/pyears ratetables keep all dates as number of days since 1/1/1960
#   convert other types of objects to this form
#
ratetableDate <- function(x) {
    UseMethod("ratetableDate", x)
    }

# Normally used in R
ratetableDate.Date <- function(x) 
    as.numeric(x - as.Date("1960/01/01"))

ratetableDate.POSIXt <- function(x)
    as.numeric(as.Date(x) - as.Date("1960/01/01"))

# Normally Splus
#ratetableDate.timeDate <- function(x)
#    as.numeric(x - timeDate('1/1/1960'))

# Therneau's old "date" class (will someday wither away)
ratetableDate.date <- function(x)  as.numeric(x)

# David James's old "chron" class (will someday wither away)
# Support it without using the chron library, which may not be loaded.
ratetableDate.chron <- function(x) {
    origin <- attr(x, "origin")
    x<- as.numeric(x) + as.Date(paste(origin["year"], origin["month"],
                                          origin["day"], sep='/'))
    ratetableDate(x)
}
ratetableDate.dates <- ratetableDate.chron

# the routines that call this are responsible for a useful error message
ratetableDate.default <- function(x) NULL 

