#
# Which is faster, a tapply or a pair of rowsums?
#

test1 <- data.frame(time=sample(1:1000000, 900000, replace=T),
                    status= rbinom(900000, 1, .3),
                    weight=runif(900000))
test1 <- test1[order(test1$time),]

fun1 <- function(x) {
    temp <- tapply(x$weight, list(x$time, x$status), sum)
    ifelse(is.na(temp), 0, temp)
    }

fun2 <- function(x) {
    temp1 <- rowsum(test1$weight[test1$status==0],
                    test1$time[test1$status==0], reorder=F)
    temp2 <- rowsum(test1$weight[test1$status==1],
                    test1$time[test1$status==1], reorder=F)
    alltime <- unique(test1$time)
    cbind(c(0, temp1)[1+match(alltime, unique(test1$time[test1$status==0]), 
                              nomatch=0)],
          c(0, temp2)[1+match(alltime, unique(test1$time[test1$status==1]), 
                               nomatch=0)])
    }

# This one doesn't work 
fun3 <- function(x) {
    ftime <- factor(x$time)
    temp1 <- rowsum(test1$weight[test1$status==0],
                    ftime[test1$status==0], reorder=F)
    temp2 <- rowsum(test1$weight[test1$status==1],
                    ftime[test1$status==1], reorder=F)
    cbind(temp1, temp2)
    }

system.time(fun1(test1))
system.time(fun2(test1))

# fun2 wins by a 15:1 ratio
