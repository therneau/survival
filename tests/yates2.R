library(survival)
# Tests for glm

gfit1 <- glm(skips ~ Mask* Opening + Solder, data=solder, poisson)
yg1 <- yates(gfit1, ~Mask)
yg2 <- yates(gfit1, ~ Mask, predict='response')
