#
# This is a test of many plots.  They get saved as xxx.pdf, and then compared
#  to xxx.pdf.save using tools:Rdiff.  Pdf files have a header that contains
#  a date, but the Rdiff sniffs that out and ignores it.
# Per the example in the test directory of R sources, use compress=FALSE
#  on the pdf to keep it as ascii.
library(survival)
library(tools)

# survfit curves
pdf("plot1.pdf", compress=FALSE)
fit1 <- survfit(Surv(time, status) ~ ph.ecog, data=lung)
plot(fit1, mark=1:4)
plot(fit1, conf.int=T, fun="event", col=1:4)
plot(fit1, mark.time=FALSE, fun="cloglog")
lines(fit1[2], mark.time=FALSE, conf.int='only', col=3, lty=2,
      fun='cloglog')
#points(fit1[3], col=6)
dev.off()
if (Rdiff("plot1.pdf", "plot1.pdf.save")) cat("plot1 differs")
