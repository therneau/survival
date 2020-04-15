pdf("fig1.pdf", width=7, height=7)
ptplot(1:3 ~ strata(sex)/(age + trt) + ns(ht/wt, df=4) / common + fixed)
dev.off()
