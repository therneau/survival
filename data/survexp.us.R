# $Id$
#
# Create the US and US-by-race hazards tables
#
#  We make it a function so that all the temp variables disappear
usfun1 <- function() {
    usd <- read.table('usdecennial.dat', header=T)
    usd <- usd[usd$year <=1990,]  # Annual tables are used from 1996 forward
    temp <- array(usd$q[usd$race=='total'], dim=c(113,2,6))
    #
    # At this point temp is an array of age, sex, year.  
    #  The first 3 ages are the q's for days 0-1, 1-7, and 7-28, the fourth
    #  is the q for the entire first year.
    # Change the array to one of daily hazard rates
    #  For the 4th row, make it so the sum of the first year's hazard is 
    #  correct,  i.e., 1*row1 + 6*row2 + 21*row3 + 337.25* row4 = -log(1-q)
    temp2      <- -log(1- temp)
    temp2[4,,] <- temp2[4,,] - (temp2[1,,] + temp2[2,,] + temp2[3,,])
    temp2[2,,] <- temp2[2,,] /6     #days 1-7
    temp2[3,,] <- temp2[3,,] /21    #days 7-28
    temp2[4,,] <- temp2[4,,] /337.25
    temp2[5:113,,] <- temp2[5:113,,]/365.25
    # Note the change from prior releases of the code.  There are 36424 days
    # per century, so I used 365.24.  However the year 2000 was the
    # exception to the exception "Every 4 years is a leap year, unless
    # divisible by 100; unless divisible by 1000".  So over the lifetime
    # that these tables will be used 365.25 is the right number.  And using
    # .24 confused everyone.

    # Decennial tables go through 2000
    #  In 1996 the single year tables were extended to higher ages and became 
    #  useful to us.
    singleyear <- 1996:2004    #change this as more data is added
    nsingle <- length(singleyear)
    temp3 <- array(0, dim=c(113, 2, 3, nsingle)) #age, sex, race, year
    for (i in 1:nsingle) {
        temp3[4:113,,,i] <- scan(paste("us", singleyear[i], ".dat", sep=''),
                                 quiet=TRUE,
                                 what=0, skip=5)/100000
        }
    temp3[4:113,,,] <- -log(1-temp3[4:113,,,]) #  hazard = -log(1-q)

    # US annual data does not divide the first year of life,
    #   we do so using another data source
    infant <- read.csv('usinfant.dat')
    iyears <- unique(infant$year)

    # extract the deaths for total (all races)
    #  then scale each year/sex group out as proportions
    deaths <- array(as.matrix(infant[,3:4]), dim=c(4, length(iyears), 2),
                dimnames=list(c('0', '1-6','7-27', '28-365'), iyears,
                              c("Male", "Female")))
    for (i in 1:length(iyears)) {
        deaths[,i,1] <- deaths[,i,1]/sum(deaths[,i,1])
        deaths[,i,2] <- deaths[,i,2]/sum(deaths[,i,2])
        }

    # Now rescale to daily hazards
    temp4 <- array(0., dim=c(113, 2, nsingle))
    indx <- match(singleyear, iyears)
    indx[singleyear < min(iyears)] <- 1
    for (i in 1:nsingle) {
        temp4[1:4, ,i] <- temp3[4, ,1,i] * deaths[,indx[i],]/c(1,6,21, 337.25)
        temp4[5:113,,i] <- temp3[5:113,,1,i] /365.25
        }
    
    #
    # Put together the rate table
    #  In prior releases the C code in pyears[123].c did interpolation on the
    #  fly.  But updating that for the current case where the index years are
    #  not evenly spaced at 10's was deemed too much work.  I now interpolate
    #  the calendar years here, and store a larger rate table.  It's about .1
    #  MByte in double precision, which no longer seems large.
    #
    years <- seq(1940, max(singleyear))
    rtable <- array(0., dim=c(113, 2, length(years)))
    xtemp <- c(194:199*10, singleyear)
    ytemp <- xtemp*0.0
    for (i in 1:113) { 
        for (j in 1:2) {  #so what if loops are slow, we only do this once
            ytemp <- c(temp2[i,j,], temp4[i,j,])
            rtable[i,j,] <- approx(xtemp, ytemp, xout=years)$y
            }
        }
                       
    if (exists('as.Date')) { # R
        datecut <- as.Date(paste(years, '/01/01', sep=''))-
                   as.Date('1960/01/01')
        datecut <- as.integer(datecut)
        }
    else if (exists('month.day.year')) { #Splus
        datecut <- julian(1,1, years, origin=c(1,1,1960))
        }
    else stop("Cannot find appropriate routine for dates")
        
    attributes(rtable) <- list(
        dim= c(113,2, length(years)),
        dimnames = list(c('0-1d','1-7d', '7-28d', '28-365d', 
              as.character(1:109)), c("male", "female"), years),
        dimid    =c("age", "sex", "year"),
        type  = c(2,1,4),
        cutpoints= list(c(0,1,7,28,1:109 * 365.25), NULL,  datecut),
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

survexp.us <- usfun1()
rm (usfun1)

# Now the table by race    
usfun2 <- function() {
    usd <- read.table('usdecennial.dat', header=T)
    usd <- usd[usd$year < 2000,]   #Annual tables are used from 1996 forward
    temp <-array(0., dim=c(113,2,2,6)) #age, sex, race, year
    # US tables were not consistent
    #   1940: total, white, nonwhite, black
    #   50-60:total, white, nonwhite
    #   70-90:total, white, nonwhite, black
    #   2000 :total, white, black
    #   1997-2004 annual tables: total, white, black
    # We have used the 1950 and 1960 nonwhite as the rate for blacks
    temp[,,1,] <- usd$q[usd$race=='white']
    temp[,,2,c(1,4,5,6)] <- usd$q[usd$race=='black']
    temp[,,2, 2:3] <- usd$q[usd$race=='nonwhite' & 
                            (usd$year==1950 |usd$year==1970)]
    #
    # At this point temp is an array of age, sex, race, year.  
    #  The first 3 ages are the q's for days 0-1, 1-7, and 7-28, the fourth
    #  is the q for the entire first year.
    # Change the array to one of daily hazard rates
    #  For the 4th row, make it so the sum of the first year's hazard is 
    #  correct,  i.e., 1*row1 + 6*row2 + 21*row3 + 337.25* row4 = -log(1-q)
    temp2      <- -log(1- temp)
    temp2[4,,,] <- temp2[4,,,] - (temp2[1,,,] + temp2[2,,,] + temp2[3,,,])
    temp2[2,,,] <- temp2[2,,,]/6      #days 2-7
    temp2[3,,,] <- temp2[3,,,]/21     #days 8-28
    temp2[4,,,] <- temp2[4,,,]/337.25    #days 29-365.25
    temp2[5:113,,,] <- temp2[5:113,,,] / 365.25

    # 
    #  In 1996 the annual tables were extended to higher ages and became 
    #  useful to us.
    singleyear <- 1996:2004    #change this as more data is added
    nsingle <- length(singleyear)
    temp3 <- array(0, dim=c(113, 2, 3, nsingle)) #age, sex, race, year
    for (i in 1:nsingle) {
        temp3[4:113,,,i] <- scan(paste("us", singleyear[i], ".dat", sep=''),
                                 quiet=TRUE,
                                 what=0, skip=5)/100000
        }
    temp3[4:113,,,] <- -log(1-temp3[4:113,,,]) #  hazard = -log(1-q)

    # US annual data does not divide the first year of life
    infant <- read.csv('usinfant.dat')
    iyears <- unique(infant$year)

    # extract the deaths for whites and blacks
    #  then scale each year/sex group out as proportions
    deaths <- array(as.matrix(infant[,5:8]), dim=c(4, length(iyears), 4),
                dimnames=list(c('0', '1-6','7-27', '28-365'), iyears,
                              c("W Male", "W Female", "B Male", "B Female")))
    for (i in 1:length(iyears)) {
        deaths[,i,1] <- deaths[,i,1]/sum(deaths[,i,1])
        deaths[,i,2] <- deaths[,i,2]/sum(deaths[,i,2])
        deaths[,i,3] <- deaths[,i,3]/sum(deaths[,i,3])
        deaths[,i,4] <- deaths[,i,4]/sum(deaths[,i,4])
        }

    # Now rescale to daily hazards
    temp4 <- array(0., dim=c(113, 2, 2, nsingle))
    indx <- match(singleyear, iyears)
    indx[singleyear < min(iyears)] <- 1
    for (i in 1:nsingle) {
        temp4[1:4,,1 ,i] <- temp3[4, ,2,i] * deaths[,indx[i],1:2]/
                           c(1,6,21, 337.25)
        temp4[1:4,,2 ,i] <- temp3[4, ,3,i] * deaths[,indx[i],3:4]/
                           c(1,6,21, 337.25)
        temp4[5:113,,,i] <- temp3[5:113,,2:3,i] /365.25
        }
    
    #
    # Put together the rate table
    #
    years <- 1940:max(singleyear)
    rtable <- array(0., dim=c(113, 2, 2, length(years)))
    xtemp <- c(194:199*10, singleyear)
    ytemp <- xtemp*0.0
    for (i in 1:113) { 
        for (j in 1:2) {  #so what if loops are slow, we only do this once
            for (k in 1:2) { #race
                ytemp <- c(temp2[i,j,k,], temp4[i,j,k,])
                rtable[i,j,k,] <- approx(xtemp, ytemp, xout=years)$y
                }
            }
        }
    if (exists('as.Date')) { # R
        datecut <- as.Date(paste(years, '/01/01', sep=''))-
                   as.Date('1960/01/01')
        datecut <- as.integer(datecut)
        }
    else if (exists('month.day.year')) { #Splus
        datecut <- julian(1,1, years, origin=c(1,1,1960))
        }
        
    attributes(rtable) <- list(
        dim= c(113,2,2, length(years)),
        dimnames = list(c('0-1d','1-7d', '7-28d', '28-365d', 
              as.character(1:109)), c("male", "female"), c("white", "black"),
                        years),
        dimid    =c("age", "sex", "race", "year"),
        type  = c(2,1,1,4),
        cutpoints= list(c(0,1,7,28,1:109 * 365.25), NULL,  NULL, datecut),
        summary = function(R) {
            x <- c(format(round(min(R[,1]) /365.25, 1)),
                   format(round(max(R[,1]) /365.25, 1)),
                   sum(R[,2]==1), sum(R[,2]==2),
                   sum(R[,3]==1), sum(R[,3]==2))
            if (is.R()) x2<- as.Date(c(min(R[,4]), max(R[,4])), 
                                    origin='1960/01/01')
            else x2 <- timeDate(julian=c(min(R[,4]), max(R[,4])))
            paste("  age ranges from", x[1], "to", x[2], "years\n",
                  " male:", x[3], " female:", x[4], "\n",
                  " date of entry from", x2[1], "to", x2[2], "\n",
                  " white:",x[7], " black:", x[8], "\n")
            })

    if (is.R()) class(rtable) <- 'ratetable'
    else        oldclass(rtable) <- 'ratetable'
    rtable
    }
survexp.usr <- usfun2()
rm(usfun2)
