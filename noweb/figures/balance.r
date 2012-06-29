# The balance figure for the concondance document
btree <- function(n) {
   tfun <- function(n, id, power) {
       if (n==1) id
       else if (n==2) c(2*id, id)
       else if (n==3) c(2*id, id, 2*id+1)
       else {
           nleft <- if (n== power*2) power  else min(power-1, n-power/2)
           c(tfun(nleft, 2*id, power/2), id,
             tfun(n-(nleft+1), 2*id +1, power/2))
           }
       }
   tfun(n, 1, 2^(floor(logb(n-1,2))))
   }

set.seed(12345)
temp <- sort(unique(floor(runif(40,0,30)))[1:13])/10
indx <- btree(13)

xpos <- 1:15
xpos[4:7] <- tapply(xpos[8:15], rep(1:4, each=2), mean)
xpos[2:3] <- tapply(xpos[4:7], rep(1:2, each=2),mean)
xpos[1] <- mean(xpos[2:3])
ypos <- rep(4:1, c(1,2,4,8))

pdf('balance.pdf', height=5, width=7)
par(mar=c(1,0,0,1) + .1)

plot(xpos, ypos, type='n', xaxt='n', yaxt='n', bty='n',
     xlab="", ylab="")
temp2 <-  c(13,7,5,3,3,3,1,1,1,1,1,1,1)
#text(xpos[indx], ypos[indx], paste(temp, " (", temp2[indx], ")",  sep=''))
text(xpos[indx], ypos[indx], temp)

delta=.1
for (i in 1:6) {
    segments(xpos[i]-delta, ypos[i]-delta,
             xpos[2*i]+delta, ypos[2*i]+delta)
    segments(xpos[i]+delta, ypos[i]-delta,
             xpos[2*i+1]-delta, ypos[2*i+1] +delta)
}
dev.off()
