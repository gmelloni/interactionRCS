library(survival)
library(rms)
library(pkgload)
data(cancer)

lung2 <- lung
lung2$ph.karno <- as.vector(scale(lung2$ph.karno))
lung2$ph.ecog <- as.vector(scale(lung2$ph.ecog))
lung2$sex <- ifelse(lung2$sex==2 , 0 , 1)
myformula <- Surv(time, status) ~ ph.karno + sex + rcs(age, 3)*sex
model <- cph(myformula , data = lung2 , x = TRUE , y=TRUE)
# Bootstrap Method returns wrong results!!!
# Here sex is recognized as a binary and no warning is triggered
myHR <- HRcubSpline( x = 40:80
             , model = model , data = lung2 , var1 ="sex", var2="age" , units=1
             , center = 0, ci=TRUE , conf = 0.95 , ci.method = "delta"
             , ci.boot.method = "norm" , R = 100 , parallel = "multicore")
myHR2 <- HRcubSpline( x = 40:80
                     , model = model , data = lung2 , var1 ="sex", var2="age" , units=1
                     , center = 0, ci=TRUE , conf = 0.95 , ci.method = "bootstrap"
                     , ci.boot.method = "perc" , R = 100 , parallel = "multicore")
par(mfrow=c(1,2))
plot.HRSpline(myHR , xlab = "Age" , main = "Delta")
plot.HRSpline(myHR2 , xlab = "Age" , main = "Bootstrap")

# Test parameters
x = c(30,50,60,70,80)
model=model
data=lung2
var1="sex"
var2="age"
units=1
center=0
ci=TRUE
conf=0.95
ci.method="bootstrap"
ci.boot.method="perc"
parallel = "multicore"
R = 100

# Interaction Model
myformula <- Surv(time, status) ~ sex*age + ph.karno + ph.ecog
modInt <- coxph(myformula , data = lung2 , x = TRUE)
modInt <- cph(myformula , data = lung2 , x = TRUE)
# # Exact
# HRSpline( x = c(30,50,60,70,80)
#           , model = modInt , data = lung2 , var1 ="sex", var2="age" , units=1
#           , center = 0, ci=TRUE , conf = 0.95 , ci.method = "exact"
#           , ci.boot.method = "norm" , R = 100 , parallel = "multicore")
# Delta
HRSpline( x = c(30,50,60,70,80)
             , model = modInt , data = lung2 , var1 ="sex", var2="age" , units=1
             , center = 0, ci=TRUE , conf = 0.95 , ci.method = "delta"
             , ci.boot.method = "norm" , R = 100 , parallel = "multicore")
# Boot
HRSpline( x = c(30,50,60,70,80)
          , model = modInt , data = lung2 , var1 ="sex", var2="age" , units=1
          , center = 0, ci=TRUE , conf = 0.95 , ci.method = "bootstrap"
          , ci.boot.method = "perc" , R = 100 , parallel = "multicore")
# Compare to cubic spline model
myformula <- Surv(time, status) ~ sex*rcs(age,3) + ph.karno + ph.ecog
modCub <- cph(myformula , data = lung2 , x = TRUE)
HRcubSpline( x = c(30,50,60,70,80)
          , model = modCub , data = lung2 , var1 ="sex", var2="age" , units=1
          , center = 0, ci=TRUE , conf = 0.95 , ci.method = "delta"
          , ci.boot.method = "norm" , R = 100 , parallel = "multicore")

# Here Delta method works
# data_for_fig1 <- readRDS("C:/Users/AI880/Dropbox/TIMI/data_for_fig1.rds") # For Andrea
data_for_fig1 <- readRDS("~/Dropbox/bwh/science/CAD/CAD_age/data_for_fig1.rds")
data_for_fig1$score.sd <- as.numeric(data_for_fig1$score.sd)
mod1 <- cph(Surv(days2mi,mifu) ~ score.sd*rcs(age,3)+sex+pc1+pc2+pc3+pc4+pc5
            , data=data_for_fig1, x=TRUE, y=TRUE)
# Takes too long for bootstrap during load_all()
myHR <- HRcubSpline( x = 40:55
             , model = mod1 , data = data_for_fig1 , var1 ="score.sd", var2="age" , units=1
             , center = 0, ci=TRUE , conf = 0.95 , ci.method = "delta"
             , ci.boot.method = "norm" , R = 10 , parallel = "multicore")
myHR2 <- HRcubSpline( x = 40:55
                     , model = mod1 , data = data_for_fig1 , var1 ="score.sd", var2="age" , units=1
                     , center = 0, ci=TRUE , conf = 0.95 , ci.method = "bootstrap"
                     , ci.boot.method = "norm" , R = 3 , parallel = "multicore")
# Params for troubleshooting
x = 30:60
model = mod1
data = data_for_fig1
var1 ="score.sd"
var2="age"
units=1
center=0
ci=TRUE
conf=0.95
ci.method="delta"
ci.boot.method="perc"
parallel = "multicore"
R = 100

# Try plot
par(mfrow=c(1,2))
plot.HRSpline(myHR , xlab = "Age" , main = "delta")
plot.HRSpline(myHR2 , xlab = "Age" , main = "bootstrap")
