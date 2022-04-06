library(survival)
#
# A Cox model with a factor, followed by survexp.  
#
pfit2 <- coxph(Surv(time, status > 0) ~ trt + log(bili) +
          log(protime) + age + platelet + sex, data = pbc)
esurv <- survexp(~ trt, ratetable = pfit2, data = pbc)

temp <- pbc
temp$sex2 <- factor(as.numeric(pbc$sex), levels=2:0,
                    labels=c("f", "m", "unknown"))
esurv2 <- survexp(~ trt, ratetable = pfit2, data = temp, 
                  rmap=list(sex=sex2))

# The call components won't match, which happen to be first
all.equal(unclass(esurv)[-1], unclass(esurv2)[-1])


# Check that the ratetableDate function is okay
#
Datedate <- function(x) {
    # Dates have an origin of 1/1/1970, dates of 1/1/1960
    offset <- as.numeric(as.Date("1970-01-01") - as.Date("1960-01-01")) 
    y <- as.numeric(x) + offset
    class(y) <- "date"
    y
}
as.data.frame.date <- as.data.frame.vector # needed to make the functions work

n <- nrow(lung)
tdata <- data.frame(age=lung$age + (1:n)/365.25,
                    sex = c('male', 'female')[lung$sex],
                    ph.ecog = lung$ph.ecog,
                    time = lung$time*3, 
                    status = lung$status,
                    entry = as.Date("1940/01/01") + (n:1)*50)
tdata$entry2 <- as.POSIXct(tdata$entry)
tdata$entry3 <- Datedate(tdata$entry)

p1 <- pyears(Surv(time, status) ~ ph.ecog, data=tdata, ratetable=survexp.us,
             rmap= list(age=age*365.25, sex=sex, year=entry))
p2 <- pyears(Surv(time, status) ~ ph.ecog, data=tdata, ratetable=survexp.us,
             rmap= list(age=age*365.25, sex=sex, year=entry2))
p3 <- pyears(Surv(time, status) ~ ph.ecog, data=tdata, ratetable=survexp.us,
             rmap= list(age=age*365.25, sex=sex, year=entry3))

all.equal(p1$expected, p2$expected)
all.equal(p1$expected, p3$expected)


# Now a ratetable with ordinary dates rather than US census style year
trate <- survexp.us
attr(trate, 'type') <- c(2,1,3)
p4 <- pyears(Surv(time, status) ~ ph.ecog, data=tdata, ratetable=trate,
             rmap= list(age=age*365.25, sex=sex, year=entry))
p5 <- pyears(Surv(time, status) ~ ph.ecog, data=tdata, ratetable=trate,
             rmap= list(age=age*365.25, sex=sex, year=entry2))
p6 <- pyears(Surv(time, status) ~ ph.ecog, data=tdata, ratetable=trate,
             rmap= list(age=age*365.25, sex=sex, year=entry3))

#all.equal(p1$expected, p4$expected) # this won't be true, US special is special
all.equal(p4$expected, p5$expected)
all.equal(p5$expected, p6$expected)
