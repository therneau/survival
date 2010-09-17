#
# Create the Minnesota tables
#
makemn <- function() {
    decennial <- read.table('minndecennial.dat', 
                            col.names=c("year", "age", "sex", "race", "q"))
    data1 <- decennial[decennial$race=='total',]
    dyear <- unique(data1$year)
    ndecen <- length(dyear)

    singleyear <- c(2000, 2004)
    nsingle <- length(singleyear)

    data2 <- array(0., dim=c(100, 2, nsingle))
    for (i in 1:nsingle) {
        data2[,,i] <- scan(paste('minn', singleyear[i], '.dat', sep=''), 
                           what=0, skip=4, quiet=TRUE)
        }

    # Create the table of raw rates
    rates <- array(0., dim=c(110,2, ndecen+nsingle))
    rates[,,1:ndecen] <- -log(1- data1$q)/365.25

    # For annual Minnesota rates from the state (2000 and 2004 data)
    #  are very bumpy.  Smooth them out using
    #  smoothing splines on the hazard ratios.  The biggest advantage of
    #  using ratios wrt 1990 (currently the last state decennial) is that
    #  a direct smooth has to "jump up" for infant mortality of age 0, which
    #  a spline doesn't like to do.  It also makes our smooth somewhat
    #  similar to what the decennial census does.
    #
    for (sex in 1:2) {
        baseline <- rates[,sex,ndecen]  #the oldest decennial
        for (j in 1:nsingle) {
            y <- -log(1- data2[,sex,j]/100000) / 365.25
            fit <- smooth.spline(0:99, y/baseline[1:100], df=20)
            rates[, sex, j+ndecen] <- predict(fit, 0:109)$y * baseline
            }
        }
    
    # Now create the actual rate table
    years <- seq(min(dyear), max(singleyear))
    if (exists('as.Date')) { # R
        datecut <- as.Date(paste(years, '/01/01', sep=''))-
                   as.Date('1960/01/01')
        datecut <- as.integer(datecut)
        }
    else if (exists('month.day.year')) { #Splus
        datecut <- julian(1,1, years, origin=c(1,1,1960))
        }

    rtable <- array(0., dim=c(110, 2, length(years)))
    tempx <- c(dyear, singleyear)
    for (i in 1:110) {
        for (j in 1:2) {
            rtable[i,j,] <- approx(tempx, rates[i,j,], years)$y
            }
        }

    attributes(rtable) <- list(
        dim= c(110, 2, length(years)),
        dimnames = list(0:109, c('male', 'female'), years),
        dimid    =c("age", "sex", "year"),
        type  = c(2,1,4),
        cutpoints= list(0:109 * 365.25, NULL,  datecut),
        summary = function(R) {
            x <- c(format(round(min(R[,1]) /365.25, 1)),
                   format(round(max(R[,1]) /365.25, 1)),
                   sum(R[,2]==1), sum(R[,2]==2))
            if (is.R()) x2<- as.Date(c(min(R[,3]), max(R[,3])), 
                                    origin='1960/01/01')
            else x2 <- timeDate(julian=c(min(R[,3]), max(R[,3])))
            paste("  age ranges from", x[1], "to", x[2], "years\n",
                  " male:", x[3], " female:", x[4], "\n",
                  " date of entry from", x2[1], "to", x2[2], "\n")
            })

    if (is.R()) class(rtable) <- 'ratetable'
    else        oldclass(rtable) <- 'ratetable'
    rtable
    }

survexp.mn <- makemn()

# The Minnesota white table
makemn2 <- function() {
    decennial <- read.table('minndecennial.dat', 
                            col.names=c("year", "age", "sex", "race", "q"))
    data1 <- decennial[decennial$race=='white',]
    dyear <- unique(data1$year)
    ndecen <- length(dyear)

    # Extrapolated 2000 data
    extrap <- read.table('minnpred.dat', header=T, sep=',', skip=5)

    # Create the table of raw rates
    rates <- array(0., dim=c(110,2, ndecen+1))
    rates[,,1:ndecen] <- -log(1- data1$q)/365.25
    rates[,,ndecen+1] <- as.matrix(extrap)

    # Make the rate table
    years <- seq(min(dyear), 2000)
    if (exists('as.Date')) { # R
        datecut <- as.Date(paste(years, '/01/01', sep='')) -
                   as.Date('1960/01/01')
        datecut <- as.integer(datecut)
        }
    else if (exists('month.day.year')) { #Splus
        datecut <- julian(1,1, years, origin=c(1,1,1960))
        }

    rtable <- array(0., dim=c(110, 2, length(years)))
    tempx <- c(dyear, 2000)
    for (i in 1:110) {
        for (j in 1:2) {
            rtable[i,j,] <- approx(tempx, rates[i,j,], years)$y
            }
        }
                    
    attributes(rtable) <- list(
        dim= c(110, 2, length(years)),
        dimnames = list(0:109, c('male', 'female'), years),
        dimid    =c("age", "sex", "year"),
        type  = c(2,1,4),
        cutpoints= list(0:109 * 365.25, NULL,  datecut),
        summary = function(R) {
            x <- c(format(round(min(R[,1]) /365.25, 1)),
                   format(round(max(R[,1]) /365.25, 1)),
                   sum(R[,2]==1), sum(R[,2]==2))
            if (is.R()) x2<- as.Date(c(min(R[,3]), max(R[,3])), 
                                    origin='1960/01/01')
            else x2 <- timeDate(julian=c(min(R[,3]), max(R[,3])))
            paste("  age ranges from", x[1], "to", x[2], "years\n",
                  " male:", x[3], " female:", x[4], "\n",
                  " date of entry from", x2[1], "to", x2[2], "\n")
            })

    if (is.R()) class(rtable) <- 'ratetable'
    else        oldclass(rtable) <- 'ratetable'
    rtable
    }
  
survexp.mnwhite <- makemn2()

rm(makemn, makemn2)

                              
            
