library(survival)

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
all.equal(y1$estimate$estimate, colMeans(tpred))

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
