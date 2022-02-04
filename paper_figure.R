library(interactionRCS)
library(survival)
library(rms)
sd_training7 <- read_excel("Z:/Minao/AHA Grant/AHA Abstract/data/sd_training7.xls")

myformula <- Surv(days2miistr, miistrfu) ~ curr_smoke*rcs(age, 3)
model_rcs <- cph(myformula , data = sd_training7 , x = TRUE , y=TRUE)

HR_rcs_delta <- rcsHR( var2values = c(50:80)
                       , model = model_rcs , data = sd_training7 , var1 ="curr_smoke", var2="age" ,ci.method = "delta")

plotHR(HR_rcs_delta , xlab = "Age")

myformula <- Surv(days2miistr, miistrfu) ~ curr_smoke*age
model_loglin <- cph(myformula , data = sd_training7 , x = TRUE , y=TRUE)

HR_loglin_delta <- loglinHR( var2values = c(50:80)
                             , model = model_loglin , data = sd_training7 , var1 ="curr_smoke", var2="age", ci.method = "delta")

plotHR(HR_loglin_delta , xlab = "Age")
