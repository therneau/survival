library(survival)
# Tests for glm

gfit1 <- glm(skips ~ Mask + Opening, data=solder, poisson)
