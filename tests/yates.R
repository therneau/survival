library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

fit1 <- lm(skips ~ Opening + Solder + Mask + PadType + Panel,
              data=solder)
y1 <- yates(fit1, "Opening")

temp <- levels(solder$Opening)
tpred <- matrix(0., nrow(solder), 3)
for (i in 1:3) {
    tdata <- solder
    tdata$Opening <- temp[i]
    tpred[,i] <- predict(fit1, newdata=tdata)
 }
all.equal(y1$estimate$mpv, colMeans(tpred))

# This fit is deficient: there are no L*A6 obs
fit2 <- lm(skips ~ Opening*Mask + Solder,
              data=solder)
y2a <- yates(fit2, "Mask", population="factorial")
y2b <- yates(fit2, "Opening",  population="factorial")
y2c <- yates(fit2, "Solder",  population="factorial")

tdata <- do.call(expand.grid, fit2$xlevels[1:3])
temp <- levels(solder$Mask)
tpreda <- matrix(0., nrow(tdata), length(temp),
                 dimnames=list(NULL, temp))
for (i in seq(along=temp)) {
    tdata$Mask <- temp[i]
    suppressWarnings(tpreda[,i] <- predict(fit2, newdata=tdata))
 }
tpreda[,"A6"] <- NA  # the A6 estimate is deficient
aeq(y2a$estimate$mpv, colMeans(tpreda))

tdata <- do.call(expand.grid, fit2$xlevels[1:3])
temp <- levels(solder$Opening)
tpredb <- matrix(0., nrow(tdata), length(temp),
                 dimnames=list(NULL, temp))
for (i in seq(along=temp)) {
    tdata$Opening <- temp[i]
    suppressWarnings(tpredb[,i] <- predict(fit2, newdata=tdata))
 }
tpredb[,"L"] <- NA  
aeq(y2b$estimate$mpv, colMeans(tpredb))

tdata <- do.call(expand.grid, fit2$xlevels[1:3])
temp <- levels(solder$Solder)
tpredc <- matrix(0., nrow(tdata), length(temp),
                 dimnames=list(NULL, temp))
for (i in seq(along=temp)) {
    tdata$Solder <- temp[i]
    suppressWarnings(tpredc[,i] <- predict(fit2, newdata=tdata))
 }
aeq(y2c$estimate$mpv, colMeans(tpredc))


index <- 1 + (1:nrow(solder)) - match(solder$Mask, solder$Mask)
solder.balance <- droplevels(subset(solder, Mask != "A6" & index <= 180))

fit2 <- lm(skips ~ Opening + Solder + Mask + PadType + Panel,
              data=solder, subset=(Mask != "A6"))
y2 <- yates(fit2, "Opening", population="factorial")
t2 <- matrix(0., nrow(solder.balance), 3)
for (i in 1:3) {
    tdata <- solder.balance
    tdata$Opening <- temp[i]
    t2[,i] <- predict(fit2, newdata=tdata)
}
all.equal(y2$estimate$estimate, colMeans(t2))
