#
# Some tests for the Yates function
#  Unbalanced, has a zero cell
tsex <- c('m','f')[lung$sex]
tfit <- lm(time ~ factor(ph.ecog)*tsex + age, lung)


# All two ways, with zero cells
tfit2 <- lm(time ~ factor(ph.ecog)*tsex + factor(ph.ecog)*factor(inst)  +
            tsex*factor(inst) + age + ns(wt.loss,3),
            lung)

