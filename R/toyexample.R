##########################
# Testing multiple knots #
##########################

# library(eventglm)
# library(survival)
# library(rms)
# data(cancer)
# cancer$sex <- ifelse(cancer$sex==2 , 0 , cancer$sex)
# myformula <- Surv(time, status) ~ age*sex
# cuminc <- cumincglm(myformula, data=cancer, time = 365, link = "identity", model.censoring="coxph")
# cumincINT <- intEST( var2values = 40:80
#                      , model = cuminc , data = cancer , var1 ="sex", var2="age"
#                      , ci=TRUE , conf = 0.95 , ci.method = "delta" , R=100)
# plotINT(cumincINT)
# library(rms)
# library(mlbench)
# data(PimaIndiansDiabetes)
# # Set age on a 5-year scale
# PimaIndiansDiabetes$age <- PimaIndiansDiabetes$age/5
# # Recode diabetes as 0/1
# PimaIndiansDiabetes$diabetes <- ifelse(PimaIndiansDiabetes$diabetes=="pos" , 1 , 0)
# # Logistic model predicting diabetes over BMI, age and glucose
# par(mfrow=c(2,2))
# for(i in 3:6){
#   myformula <- diabetes ~ mass + age * rcs( glucose , i)
#   model <- glm(myformula , data = PimaIndiansDiabetes , family = "binomial")
#   plotINT(rcsOR( var2values = 20:150
#                  , model = model , data = PimaIndiansDiabetes , var1 ="age", var2="glucose"
#                  , ci=TRUE , conf = 0.95 , ci.method = "bootstrap" , ci.boot.method = "norm") , ylim = c(0,3) , main = i)
# }
# # Linear model predicting BMI over diabetes, age and glucose
# par(mfrow=c(2,2))
# for(i in 3:6){
#   myformula2 <- mass ~ diabetes + age * rcs( glucose , i )
#   model2 <- glm(myformula2 , data = PimaIndiansDiabetes , family = "gaussian")
#   plotINT(intEST( var2values = 20:150
#                   , model = model2 , data = PimaIndiansDiabetes , var1 ="age", var2="glucose"
#                   , ci=TRUE , conf = 0.95 , ci.method = "bootstrap") , ylim = c(-1,3) , main = i)
# }
# # Cox model
# library(survival)
# library(rms)
# data(cancer)
# cancer$sex <- ifelse(cancer$sex==2 , 0 , cancer$sex)
# par(mfrow=c(2,2))
# for(i in 3:6){
#   myformula <- Surv(time, status) ~ ph.karno + ph.ecog + rcs(age,i)*sex
#   model <- cph(myformula , data = lung )
#   plotINT(intEST( var2values = 40:80
#                   , model = model , data = lung , var1 ="sex", var2="age"
#                   , ci=TRUE , conf = 0.95 , ci.method = "bootstrap" , R=100), ylim = c(-1,3) , main = i)
# }

# library(survival)
# library(rms)
# data(cancer)
# #----------------------------#
# # rcs continuos on continuos #
# #----------------------------#
# lung2 <- lung
# lung2$ph.karno <- as.vector(scale(lung2$ph.karno))
# lung2$ph.ecog <- as.vector(scale(lung2$ph.ecog))
# lung2$sex <- ifelse(lung2$sex==2 , 0 , 1)
# myformula <- Surv(time, status) ~ ph.ecog + age + rcs(ph.karno, 3)*sex
# model <- cph(myformula , data = lung2 , x = TRUE , y=TRUE)
# # Here sex is recognized as a binary and no warning is triggered
# myHR <- rcsHR( var2values =  seq(-3 , 2 , 0.5)
#              , model = model , data = lung2 , var1 ="sex", var2="ph.karno" , units=1
#              , center = 0, ci=TRUE , conf = 0.95 , ci.method = "delta"
#              , ci.boot.method = "norm" , R = 100 , parallel = "multicore")
# myHR2 <- rcsHR( var2values =  seq(-3 , 2 , 0.5)
#                      , model = model , data = lung2 , var1 ="sex", var2="ph.karno" , units=1
#                      , center = 0, ci=TRUE , conf = 0.95 , ci.method = "bootstrap"
#                      , ci.boot.method = "perc" , R = 100 , parallel = "multicore")
# par(mfrow=c(1,2))
# plotHR(myHR , xlab = "ph.karno" , main = "Delta" , log=TRUE)
# plotHR(myHR2 , xlab = "ph.karno" , main = "Bootstrap" , log=TRUE)
# #-------------------#
# # Interaction Model #
# #-------------------#
# myformula <- Surv(time, status) ~ sex*age + ph.karno + ph.ecog
# # Simple interaction can accept coxph or cph models. The result is the same
# modInt <- coxph(myformula , data = lung2 , x = TRUE)
# modInt <- cph(myformula , data = lung2 , x = TRUE)
# # # Exact
# # loglinHR( var2values =  c(30,50,60,70,80)
# #           , model = modInt , data = lung2 , var1 ="sex", var2="age" , units=1
# #           , center = 0, ci=TRUE , conf = 0.95 , ci.method = "normal"
# #           , ci.boot.method = "norm" , R = 100 , parallel = "multicore")
# # Delta
# loglinHR( var2values =  c(30,50,60,70,80)
#              , model = modInt , data = lung2 , var1 ="sex", var2="age" , units=1
#              , center = 0, ci=TRUE , conf = 0.95 , ci.method = "delta"
#              , ci.boot.method = "norm" , R = 100 , parallel = "multicore")
# # Boot
# loglinHR( var2values =  c(30,50,60,70,80)
#           , model = modInt , data = lung2 , var1 ="sex", var2="age" , units=1
#           , center = 0, ci=TRUE , conf = 0.95 , ci.method = "bootstrap"
#           , ci.boot.method = "perc" , R = 100 , parallel = "multicore")
# # Compare to cubic spline model
# myformula <- Surv(time, status) ~ sex*rcs(age,3) + ph.karno + ph.ecog
# modCub <- cph(myformula , data = lung2 , x = TRUE)
# rcsHR( var2values =  c(30,50,60,70,80)
#           , model = modCub , data = lung2 , var1 ="sex", var2="age" , units=1
#           , center = 0, ci=TRUE , conf = 0.95 , ci.method = "delta"
#           , ci.boot.method = "norm" , R = 100 , parallel = "multicore")
#
# # Example 3 on UKBB data
# # data_for_fig1 <- readRDS("C:/Users/AI880/Dropbox/TIMI/data_for_fig1.rds") # For Andrea
# data_for_fig1 <- readRDS("~/Dropbox/bwh/science/CAD/CAD_age/data_for_fig1.rds")
# data_for_fig1$score.sd <- as.numeric(data_for_fig1$score.sd)
# mod1 <- cph(Surv(days2mi,mifu) ~ score.sd*rcs(age,3)+sex+pc1+pc2+pc3+pc4+pc5
#             , data=data_for_fig1, x=TRUE, y=TRUE)
# # Takes a bit of time for bootstrap during load_all()
# myHR <- rcsHR( var2values =  40:55
#              , model = mod1 , data = data_for_fig1 , var1 ="score.sd", var2="age" , units=1
#              , center = 0, ci=TRUE , conf = 0.95 , ci.method = "delta"
#              , ci.boot.method = "norm" , R = 10 , parallel = "multicore")
# myHR2 <- rcsHR( var2values =  40:55
#                      , model = mod1 , data = data_for_fig1 , var1 ="score.sd", var2="age" , units=1
#                      , center = 0, ci=TRUE , conf = 0.95 , ci.method = "bootstrap"
#                      , ci.boot.method = "norm" , R = 3 , parallel = "multicore")
# # Try plotHR
# par(mfrow=c(1,2))
# plotHR(myHR , xlab = "Age" , main = "delta")
# plotHR(myHR2 , xlab = "Age" , main = "bootstrap")
#
# # Params for troubleshooting
# var2values =  30:60
# model = mod1
# data = data_for_fig1
# var1 ="score.sd"
# var2="age"
# units=1
# center=0
# ci=TRUE
# conf=0.95
# ci.method="delta"
# ci.boot.method="perc"
# parallel = "multicore"
# R = 100
#
# # Test parameters
# var2values =  c(30,50,60,70,80)
# model=model
# data=lung2
# var1="sex"
# var2="age"
# units=1
# center=0
# ci=TRUE
# conf=0.95
# ci.method="bootstrap"
# ci.boot.method="perc"
# parallel = "multicore"
# R = 100

#------------#
# OR example #
#------------#
# library(mlbench)
# data(PimaIndiansDiabetes)
# myformula <- diabetes ~ mass + age * rcs(glucose, 3)
# model <- lrm(myformula , data = PimaIndiansDiabetes )
# var2values = 20:50
# model = model
# data = PimaIndiansDiabetes
# var1 ="age"
# var2="glucose"
# ci=TRUE
# conf = 0.95
# ci.method = "delta"
