# Verify the values found in the Royston paper
library(survival)

pbc2 <- na.omit(pbc[,-1])  # no id variable, no missings

pfit1 <- coxph(Surv(time, status==2) ~ . + log(bili) - bili, pbc2,
               ties="breslow")
# backwards elimination was used to eliminate all but 8
pfit2 <- coxph(Surv(time, status==2) ~ age + log(bili) + edema + albumin +
                   stage + copper, data=pbc2, ties="breslow")

temp <- rbind(royston(pfit1), royston(pfit1, adjust=TRUE),
              royston(pfit2), royston(pfit2, adjust=TRUE))
all.equal(round(temp[,1], 2), c(2.86, 2.56, 2.69, 2.59))
