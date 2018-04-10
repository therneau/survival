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
all.equal(y1$estimate[,"pmm"], colMeans(tpred))

# This fit is deficient: there are no Opening=L and Mask=A6 obs
# The MPV for Mask=A6 and Opening L will therefore be NA, as well
#   as for all levels of Solder, but we can compute the others.
# Solder will be NA for all levels
fit2 <- lm(skips ~ Opening*Mask + Solder,
              data=solder)
y2a <- yates(fit2, "Mask", population="factorial")
y2b <- yates(fit2, "Opening",  population="factorial")
y2c <- yates(fit2, "Solder",  population="factorial")

# The predict.lm function gives correct predictions for estimable
#  functions (all but L,A6) and nonsense for others.  It knows that
#  some are not estimable due to the NA coefficients, but not which ones,
#  so always prints a warning.  Hence the suppressWarnings call.
tdata <- do.call(expand.grid, fit2$xlevels[1:3])
temp <- levels(solder$Mask)
tpreda <- matrix(0., nrow(tdata), length(temp),
                 dimnames=list(NULL, temp))
for (i in seq(along=temp)) {
    tdata$Mask <- temp[i]
    suppressWarnings(tpreda[,i] <- predict(fit2, newdata=tdata))
 }
tpreda[,"A6"] <- NA  # the A6 estimate is deficient
aeq(y2a$estimate[,"pmm"], colMeans(tpreda))

tdata <- do.call(expand.grid, fit2$xlevels[1:3])
temp <- levels(solder$Opening)
tpredb <- matrix(0., nrow(tdata), length(temp),
                 dimnames=list(NULL, temp))
for (i in seq(along=temp)) {
    tdata$Opening <- temp[i]
    suppressWarnings(tpredb[,i] <- predict(fit2, newdata=tdata))
 }
tpredb[,"L"] <- NA  
aeq(y2b$estimate[,"pmm"], colMeans(tpredb))

# Solder should be all NA
all(is.na(y2c$estimate[,"pmm"]))

# Tests for Solder are defined for a non-factorial population, however.
# the [] below retains the factor structure of the variable, where the
#  runs above did not.  R gets prediction correct  both ways.
y2d <- yates(fit2, ~Solder)
temp <- levels(solder$Solder)
tdata <- solder
tpredd <- matrix(0, nrow(tdata), length(temp),
                 dimnames=list(NULL, temp))
for (i in seq(along=temp)) {
    tdata$Solder[] <- temp[i]
    suppressWarnings(tpredd[,i] <- predict(fit2, newdata=tdata))
}
aeq(y2d$estimate$pmm, colMeans(tpredd))

#
# Verify that the result is unchanged by how dummies are coded
#   The coefs move all over the map, but predictions are unchanged
fit3 <- lm(skips ~ C(Opening, contr.helmert)*Mask + C(Solder, contr.SAS),
           data=solder)
y3a <- yates(fit3, ~Mask, population='yates')
equal <- c("estimate", "test", "mvar")
all.equal(y3a[equal], y2a[equal])

tdata <- do.call(expand.grid, fit2$xlevels[1:3]) # use orignal variable names
temp <- levels(solder$Mask)
cpred <- matrix(0., nrow(tdata), length(temp),
                 dimnames=list(NULL, temp))
for (i in seq(along=temp)) {
    tdata$Mask <- temp[i]
    suppressWarnings(cpred[,i] <- predict(fit3, newdata=tdata))
 }
aeq(cpred[, temp!="A6"], tpreda[, temp!= "A6"])  # same predictions
all.equal(y3a$estimate, y2a$estimate)

y3b <- yates(fit3, ~Opening, population='yates')
# column names will differ
all.equal(y3b$estimate, y2b$estimate, check.attributes=FALSE)

y3d <- yates(fit3, ~Solder)
for (i in 1:3) {
    print(all.equal(y3d[[i]], y2d[[i]], check.attributes=FALSE))
}

# Reprise this with a character variable in the model
sdata <- solder
sdata$Mask <- as.character(sdata$Mask)
fit4 <-  lm(skips ~ Opening*Mask + Solder, data=sdata)
y4a <- yates(fit4, ~ Mask, population= "yates")
y4b <- yates(fit4, ~ Opening, population= "yates")
y4d <- yates(fit4, ~ Solder)
equal <- c("estimate", "tests", "mvar", "cmat")
all.equal(y2a[equal], y4a[equal])  # the "call" component differs
all.equal(y2b[equal], y4b[equal])
all.equal(y2d[equal], y4d[equal])
