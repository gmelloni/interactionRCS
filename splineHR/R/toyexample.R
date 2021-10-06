library(survival)
library(rms)
data(cancer)
x <- 55:60

myformula <- Surv(time, status) ~ ph.karno*rcs(age, 3) + sex
model <- cph(myformula , data = lung , x = TRUE , y=TRUE)
HRcubSpline( x = 90:100
             , model = model , data = lung , var1 ="ph.karno", var2="age" , units=1
             , center = 0, ci=TRUE , conf = 0.95 , ci.method = "bootstrap"
             , ci.boot.method = "perc" , R = 100 , parallel = "multicore")
# Test parameters
x = 90:100
model = model
data = lung
var1 ="ph.karno"
var2="age"
units=1
center=0
ci=TRUE
conf=0.95
ci.method="bootstrap"
ci.boot.method="perc"
parallel = "multicore"
R = 100
# var1 <- "sex"

# Interaction Model
myformula <- Surv(time, status) ~ age*ph.karno + sex
model <- coxph(myformula , data = lung , x = TRUE)
