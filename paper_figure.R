library(interactionRCS)
library(survival)
library(rms)
library(readxl)
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



# with coxph, good
model_loglin <- coxph(myformula , data = sd_training7 )

HR_loglin_delta2 <- loglinHR( var2values = c(50:80)
                             , model = model_loglin , data = sd_training7 , var1 ="curr_smoke", var2="age", ci.method = "delta")

plotHR(HR_loglin_delta2 , xlab = "Age")


###### logistic


# rcs needs rms::lrm
modellog <- rms::lrm( hxcad~ curr_smoke*rcs(age, 3),data = sd_training7)

OR_rcs_delta <- rcsOR( var2values = c(50:80)
                          , model = modellog , data = sd_training7 , var1 ="curr_smoke", var2="age" ,ci.method = "delta")

plotINT(OR_rcs_delta , xlab = "Age")




# log linear good with both
modellog <- rms::lrm( hxcad~ curr_smoke*age,data = sd_training7)

OR_loglin_delta <- loglinOR( var2values = c(50:80)
                        , model = modellog , data = sd_training7 , var1 ="curr_smoke", var2="age" ,ci.method = "delta")

plotINT(OR_loglin_delta , xlab = "Age")



modellog <- glm( hxcad~ curr_smoke*age,data = sd_training7,family="binomial")



OR_loglin_delta2 <- loglinOR( var2values = c(50:80)
                          , model = modellog , data = sd_training7 , var1 ="curr_smoke", var2="age" ,ci.method = "delta")

plotINT(OR_loglin_delta2 , xlab = "Age")




#### linear model

# rcs needs rms::Glm

modellin <- rms::Glm( BLGFR~ curr_smoke*rcs(age, 3),data = sd_training7)



lin_rcs_delta <- rcsLIN( var2values = c(50:80)
                        , model = modellin , data = sd_training7 , var1 ="curr_smoke", var2="age" ,ci.method = "delta")

plotINT(lin_rcs_delta , xlab = "Age")

# linear good with the second one
modellin <- rms::Glm( BLGFR~ curr_smoke*age,data = sd_training7)

lin_delta <- linLIN( var2values = c(50:80)
                        , model = modellin , data = sd_training7 , var1 ="curr_smoke", var2="age" ,ci.method = "delta")

plotINT(lin_delta , xlab = "Age")

#
modellin <- glm( BLGFR~ curr_smoke*age,data = sd_training7)



lin_delta2 <- linLIN( var2values = c(50:80)
                        , model = modellin , data = sd_training7 , var1 ="curr_smoke", var2="age" ,ci.method = "delta")

plotINT(lin_delta2 , xlab = "Age")


