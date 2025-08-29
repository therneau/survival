# Better control of repeats
library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

# When there are multiple censored rows, don't add more than we need
baseline <- data.frame(idd=1:5,  sex=c('f', 'm')[c(1,2,1,2,2)], 
                       futime=c(20, 30, 40, 30, 20),
                       status= c(0, 1, 0, 1, 0))
d2 <- data.frame(idd= rep(2:5, each=3), 
                 day=c(0,15, 21, 0,25, 31, 0,20, 25, 0,10, 12),
                 cmc=c(2,3,3, 3,3,4, 2,3,3, 3,4,2))

test1 <- tmerge(baseline[,1:2], baseline, id=idd, death=event(futime, status))
test2 <- tmerge(test1, d2, id=idd, cmc= tdc(day, cmc))

# idd 2 has an obs at (0, 30), add new one at 21
# idd 3
aeq(test2$idd, c(1,2,2,3,3,4,4,5,5,5))
aeq(test2$cmc, c(NA, 2,3,3,4,2,3,3,4,2)) 
