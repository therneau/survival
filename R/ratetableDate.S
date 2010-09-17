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
ratetableDate.timeDate <- function(x)
    as.numeric(x - timeDate('1/1/1960'))

# Therneau's old "date" class (will someday wither away)
ratetableDate.date <- function(x)  as.numeric(x)

# David James's old "chron" class (will someday wither away)
ratetableDate.chron <- function(x)
    as.numeric(x - chron("01/01/1960"))

# the routines that call this are responsible for a useful error message
ratetableDate.default <- function(x) NULL 

