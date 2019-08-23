#
# This code verified some oddities about model frames
#
options(contrasts=c("contr.treatment", "contr.poly")) # clean slate

tdata <- data.frame(y=1:10, x1=letters[c(1,2,3,1,2,3,1,2,3,1)],
                            x2=LETTERS[c(1,2,3,4,1,2,3,4,1,2)],
                    stringsAsFactors=TRUE)
tdata$x3 <- as.character(tdata$x1)

fit1 <- lm(y ~ x1 + x2, tdata, x=TRUE)
m1 <- fit1$model
t1 <- terms(fit1)

# Lesson 1: xlev is ignored when the variable already has levels
#  as an attribute

temp <- list(x1=c("a", "c", "b"), x2=LETTERS[4:1])
x2 <- model.matrix(t1, m1, xlev=temp)
x3 <- model.matrix(t1, tdata, xlev=temp)

all.equal(x2, fit1$x)
!is.logical(all.equal(x2, x3))  # x2 and x3 do not agree
attributes(m1$x1)


# Lesson 2: character variables do not have their levels
# remembered as attributes in the model frame, but these are
#   found in fit$xlevels and fit$contrasts.
# However, the xlev argument is still ignored for a model frame!
fit2 <- lm(y ~ x3 + x2, tdata, x=TRUE)
m2 <- fit2$model
x3 <- model.matrix(terms(fit2), m2, xlev=fit2$xlevels)
x4 <- model.matrix(terms(fit2), m2, xlev=list(x3=letters[3:1]))
x5 <- model.matrix(terms(fit2), tdata, xlev=list(x3=letters[3:1]))
all.equal(fit2$x, x3)
all.equal(x3, x4)
all.equal(x3, x5)  # FALSE


# Lesson 3: contrasts.arg is relevant, even when the model frame
#  has a saved contrast
ctemp <- list(x1="contr.SAS", x2= contr.helmert(LETTERS[1:4]))
x4 <- model.matrix(t1, m1, contrasts.arg=ctemp)  # no saved contrast

fit3 <- lm(y ~ x1 + C(x2, contr.SAS), tdata)   
m3 <- fit3$model
attr(m3[[3]], 'contr')   # the contrast is saved
c2 <- ctemp
names(c2) <- names(fit3$contrasts)
x5 <- model.matrix(terms(fit3), m3, contrasts.arg=c2)
all.equal(x4, x5, check.attributes=FALSE)

# Lesson 4: and this holds for a character variable as well
fit4 <- lm(y ~ x3 + C(x2, contr.SAS), tdata, x=TRUE)   
c4 <- ctemp
names(c4) <- names(fit4$contrasts)
x6 <- model.matrix(terms(fit4), fit4$model, xlev=fit4$xlevels, contrasts.arg=c4)
x7 <- model.matrix(terms(fit4), fit4$model, contrasts.arg=c4)
all.equal(x7, x6)
all.equal(x7, x5, check.attributes=FALSE)
