# Modelling Tuberculosis Incidence in Finland

---
title: "Thesis Analysis"
author: "Charles Ocran (MSc.)"
co-author: "Leena Kalliovirta (PhD)"
supervisor: "Matti J. Pirinen (Prof.)"
date: "2025-09-30"
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
library(readxl)
inf_data = read_excel("C:/Users/charl/OneDrive - Laurea-ammattikorkeakoulu/Desktop/Research/infectious_diseases.xlsx")

# Data description
#head(inf_data)
inf_data
``` 


```{r}
data = inf_data[c(2,3,5,6), c(-32,-33)]
dim(data)
```
##############
Extract DATA 
#############
```{r}
#x1 = t(tot.dat[1, ])
x2 = t(data[1, ])
x3 = t(data[2, ])

x2 = as.numeric(x2[-1])
x3 = as.numeric(x3[-1])
```



```{r, fig.height=5, fig.width=5}
X1 = ts(cbind(x2, x3), start = c(1995, 1), 
        frequency = 1)
ts.plot(X1, ylim = c(0,450), 
        col = 2:3, lty = 1, 
        lwd = 2, xlab = "Time(years)", 
        ylab = "Cases", 
        main = "Tuberculosis Subtypes")

legend("topright", legend = c(expression(paste("P-Tb", ~("y"[1]))),
expression(paste("O-Tb", ~("y"[2])))),
       col = 2:3, lty = 1, lwd = 2)

#legend("topright", legend = c("P-Tb (y1)", "O-Tb (y2)"),
       #col = 2:3, lty = 1, lwd = 2)
```


```{r, fig.height=4, fig.width=7}
# Plot the individual Series
par(mfrow = c(1,2))

x2 = ts(x2, start = c(1995, 1), frequency = 1)
x3 = ts(x3, start = c(1995, 1), frequency = 1)

plot.ts(x2, ylim = c(0,450), xlab = "Time (years)", ylab = "Cases",
        main = expression(paste("Subtype P-Tb", ~(italic("y")[1]))),
        col = 2)
plot.ts(x3, ylim = c(0,230), xlab = "Time (years)", ylab = "Cases",
        main = expression(paste("Subtype O-Tb", ~(italic("y")[2]))),
        col = 7)
```

################################
      Preliminary analysis
################################

```{r}
# Checking if subtype data relate
ols.m = lm(x2 ~ x3) #for Tb cases
summary(ols.m)
```

```{r}
# Report Confidence Intervals of Estimates
confint(ols.m)
```

```{r, fig.height=5, fig.width=6}
acf(x2, main = "ACF of P-Tb") # P-Tb
pacf(x2, main = "PACF of P-Tb") 
```


```{r}
acf(x3, main = "ACF of O-Tb")
pacf(x3, main = "PACF of O-Tb")
```


#######################
# Co-integration Test #
#######################
```{r}
library(urca)
ols.m.resid = ols.m$residuals
co.y = ur.df(ols.m.resid, type = "none", selectlags = "AIC")
co.y@teststat
co.y@cval
```

```{r}
#############
  Solution
#############

#Since the critical value in absolute terms is greater than the absolute critical value at #$5\%$ significance level, we can reject the null hypothesis of non-stationarity of the #series.

#That implies there is a long run relationship between the two variables $y1$ and $y2$ in our investigation.
#Thus the series are co-integrated!
```


#######################################
      PART 2 - UNIVARIATE ANALYSIS
#######################################

```{r}
# Fitting AR(1) + trend + constant model for Tb data
# Fit lm for Tb subtype 1
#t3 = 3:30
#y3 = x2 # start point
#y23 = x2[-c(29,30)] # lag 2

# y1 data

t = 2:30
y = x2[-c(1)] # start point
yl1 = x2[-c(30)] # lag 1
cons = rep(1, 29)  # constant term
trend = 1:29  # trend term
Y = cbind(y, yl1, trend)   #  y23,
cnames = c("y", "yl1", "trend") # "yl2",
y.dat = data.frame(Y)
colnames(y.dat) = cnames
head(y.dat)
```


Let's fit the model next using lm():

```{r}
# AR(1) Fit
fit.lm1 = lm(y ~ ., data = y.dat)
summary(fit.lm1)
```

```{r}
# Report Confidence Intervals of Estimates
confint(fit.lm1)
```

```{r, fig.height=4, fig.width=9}
res_lm1 = residuals(fit.lm1)

# Visualize acf & pacf of residuals
par(mfrow = c(1, 2))
acf(res_lm1, main = "Residuals ACF")
pacf(res_lm1, main = "Residuals PACF")
```

```{r}
library(forecast)
checkresiduals(res_lm1, ylab = "Residuals of y1-Model M2")
```

######################


```{r}
qqnorm(rstudent(fit.lm1))
abline(h=0, qqnorm(rstudent(fit.lm1)))
```


#########################
  MAKE THE PREDICTIONS
#########################

```{r, fig.height=5, fig.width=6}
# Predict with data points
x.vals = seq(min(x2), max(x2), length = 29)
y.vals = predict(fit.lm1, newdata = data.frame(x2 = x.vals))

Pred.P_Tb = ts(y.vals, start = c(1995, 1), frequency = 1)
y.b = cbind(y, Pred.P_Tb)
ts.plot(y.b, ylim = c(0,450), xlab = "Time (years)", 
        ylab = "Cases of P-Tb",
        main = expression(paste("Adjusted AR(1) Predictions of Target", ~italic("y"[1]))), 
        col = c(5,7))
legend("topright", legend = c("M2 P-Tb", "M2 Pred P-Tb"), 
       lty = 1, lwd = 2, col = c(5,7))
```


```{r}
# Print the predictions
print(Pred.P_Tb)
```

###################
Check the errors ##
###################

```{r}
# Mean Absolute Error
abs.e = abs(y - Pred.P_Tb)
err.1 = mean(abs.e)
paste("The mean absolute errors of predictions = ", err.1)
```

```{r}
paste("The mean square error of the predictions = ", mean((Pred.P_Tb - y)^2))
```


```{r}
# Extrapolate with the pred. model
# let's see how much the cases were in 1994

#n1994 = (389.74639 - 225.4457 + 5.9146) / 0.4049

n1994 = (395.66103 - 225.4457 + 5.9146) / 0.4049
paste("The value in the subtype cases in 1994 = ", n1994)
```

########################################
PREDICT NEW DATA Using Predicted Values 
########################################

```{r}
# What about 2024/ 2025/ 2026
#n2024 = 225.4457 + (0.4049*92.95008) - 5.9146*31
#n2025 = 225.4457 + (0.4049*79.72859) - 5.9146*32 
#n2026 = 225.4457 + (0.4049*68.46061) - 5.9146*33

n2024 = 225.4457 + (0.4049*98.86472) - 5.9146*30
n2025 = 225.4457 + (0.4049*n2024) - 5.9146*31 
n2026 = 225.4457 + (0.4049*n2025) - 5.9146*32
n2027 = 225.4457 + (0.4049*n2026) - 5.9146*33

C1 = cbind(n1994, n2024, n2025, n2026, n2027) 
rownames(C1) = c("x2 Cases")
print(C1)

```

#############################
  Upgrade to AR2 Model - y1
#############################

```{r}
# AR2 + const + trend model

t = 3:30
y.y = x2[-c(1, 2)] # start point
yl.1 = x2[-c(1, 30)] # lag 1
yl.2 = x2[-c(29, 30)]
cons = rep(1, 28)  # constant term
tr.end = 1:28  # trend term
Y.2 = cbind(y.y, yl.1, yl.2, tr.end)   #  y23,
cnames2 = c("ys", "yL1", "yL2", "trend") 
y.Dat = data.frame(Y.2)
colnames(y.Dat) = cnames2
head(y.Dat)
```

#####################
Fit Model
```{r}
AR2.ct.model = lm(ys ~ ., data = y.Dat)
summary(AR2.ct.model)
```

```{r}
# Report Confidence Intervals of x2-Estimates
confint(AR2.ct.model)
```

###################
See Residuals next
###################

```{r, fig.height=4, fig.width=8}
ar2.ct.residuals = residuals(AR2.ct.model)

# plots
par(mfrow = c(1, 2))
acf(ar2.ct.residuals, main = "ACF of y1-Model M4")
pacf(ar2.ct.residuals, main = "PACF of y1-Model M4")
```

```{r}
checkresiduals(ar2.ct.residuals, ylab = "Residuals of y1-Model M4")
```


##########################
   MAKE THE PREDICTIONS
##########################

```{r, fig.height=5, fig.width=6}

# Predict with data points
x.vals.2 = seq(min(x2), max(x2), length = 28)
y.vals.2 = predict(AR2.ct.model, data = data.frame(x2 = x.vals.2))

Pred.P_Tb2 = ts(y.vals.2, start = c(1995, 1), frequency = 1)
y.b2 = cbind(y.y, Pred.P_Tb2)
ts.plot(y.b2, ylim = c(0,450), 
        xlab = "Time (years)", ylab = "Cases of P-Tb",
        main = expression(paste("AR(2) Predictions of Target", ~italic("y"[1]))),
        col = c(5,7))
legend("topright", legend = c("M4 P-Tb", "M4 Pred P-Tb"), 
       lty = 1, lwd = 2, col = c(5,7))
```

#####################
Check the errors ####
#####################

```{r}
# Mean Absolute Error
abs.e2 = abs(y.y - Pred.P_Tb2)
err.2 = mean(abs.e2)
paste("The mean absolute errors of predictions = ", err.2)
```


```{r}
paste("The mean square error of the predictions = ", mean((Pred.P_Tb2 - y.y)^2))
```



####################
DO NEW PREDICTIONS
####################
```{r}
# z = 145.75 + 0.1784*z1 + 0.3783*z2 - 3.7669*tr
# Predict 2023/ 2024/ 2025/ 2026/ 2027/ 2028/ 2029

# 2023
z1 = Pred.P_Tb2[28]
z2 = Pred.P_Tb2[27]
tr = 29
C23 = 145.75 + 0.1784*z1 + 0.3783*z2 - 3.7669*tr

# 2024
z29 = C23
z28 = Pred.P_Tb2[28]
tr2 = 30
C24 = 145.75 + 0.1784*z29 + 0.3783*z28 - 3.7669*tr2


# 2025
z30 = C24
z29 = z29
tr31 = 31
C25 = 145.75 + 0.1784*z30 + 0.3783*z29 - 3.7669*tr31

# 2026
z31 = C25
z30 = z30
tr32 = 32
C26 = 145.75 + 0.1784*z31 + 0.3783*z30 - 3.7669*tr32 #C26 = z32

# 2027
tr33 = 33
C27 = 145.75 + 0.1784*C26 + 0.3783*C25 - 3.7669*tr33

# 2028
tr34 = 34
C28 = 145.75 + 0.1784*C27 + 0.3783*C26 - 3.7669*tr34

# 2029
tr35 = 35
C29 = 145.75 + 0.1784*C28 + 0.3783*C27 - 3.7669*tr35


Z = cbind(C23, C24, C25, C26)
rownames(Z) = c("P-Tb Cases")
print(Z)
```

################################
  Begin Analysis for Target x3 
################################
```{r}
t = 2:30
y2y = x3[-c(1)]  # start
y3.l1 = x3[-c(30)] # l1
trenD = 1:29  # trend term
Y.33 = cbind(y2y, y3.l1, trenD)
conames3 = c("y2y", "yyl1", "trnd")
y.data3 = data.frame(Y.33)
colnames(y.data3) = conames3
head(y.data3)
```

```{r}
# Fit the 1st Case
ar1.lm = lm(y2y ~ ., data = y.data3)
summary(ar1.lm)
```

```{r}
confint(ar1.lm)
```

###############
check residuals
```{r}
res_ar1.lm = resid(ar1.lm)
checkresiduals(ar1.lm, xlab = "Time (years)", ylab = "Residuals of y2-Model M3")
```

############
PREDICTIONS
############

```{r, fig.height=5, fig.width=6}

# Predict with data points
x.vals3 = seq(min(x3), max(x3), length = 29)
y.vals3 = predict(ar1.lm, data = data.frame(x3 = x.vals3))

Pred.O_Tb1 = ts(y.vals3, start = c(1995, 1), frequency = 1)
y.b1 = cbind(Pred.O_Tb1, ts(y2y, start = c(1995, 1), frequency = 1))
ts.plot(y.b1, ylim = c(0,230), xlab = "Time (years)", ylab = "Cases of O-Tb",
        main = expression(paste("Adjusted AR(1) Predictions of Target", ~italic("y"[2]))), col = c(3,7))
legend("topright", legend = c("M3 O-Tb", "M3 Pred O-Tb"), 
       lty = 1, lwd = 2, col = c(3,7))
```

#############
Check Errors
#############

```{r}
err.pred.ar1 = abs(Pred.O_Tb1 - y2y)
mae.ar1.x3 = mean(err.pred.ar1) 
paste("The prediction error quantity (MAE) is computed to be = ", mae.ar1.x3)
```

```{r}
paste("The prediction error quantity (MSE) is computed to be = ", mean((Pred.O_Tb1 - y2y)^2))
```


#####################
  PREDICT NEW CASES
#####################
```{r}
# Look back in 1994
# yt = 44.4892 + 0.7231yl1 -1.2448tr

c1994 = (Pred.O_Tb1[1] - 44.892 + 1.2448) / 0.7231
paste("The predicted cases back in 1994 = ", c1994)
```

```{r}
# print the prediction values
print(Pred.O_Tb1)
```

####################
  More Predictions
####################

```{r}
c2024 = 44.4892 + 0.7231*Pred.O_Tb1[29] - 1.2448*30
c2025 = 44.4892 + 0.7231*c2024 - 1.2448*31
c2026 = 44.4892 + 0.7231*c2025 - 1.2448*32
c2027 = 44.4892 + 0.7231*c2026 - 1.2448*33
c2028 = 44.4892 + 0.7231*c2027 - 1.2448*34
c2029 = 44.4892 + 0.7231*c2028 - 1.2448*35

CD = cbind(c1994, c2024, c2025, c2026)
rownames(CD) = c("Cases of O-Tb")
print(CD)
```


############################
  AR2 Model-X3 Begins here
############################

Similarly, we shall do the same analysis for data x3:

```{r}
t = 3:30
yy = x3[-c(1,2)]  # start
y3l1 = x3[-c(1,30)] # l1
y3l2 = x3[-c(29,30)] # l2
cons = rep(1, 28)  # const term
tre = 1:28  # trend term
Y33 = cbind(yy, y3l1, y3l2, tre)
conames = c("yy", "yyl1", "yyl2", "trd")
y.data = data.frame(Y33)
colnames(y.data) = conames
head(y.data)
```

############
Fit Model
```{r}
f.lm2 = lm(yy ~ ., data = y.data)
summary(f.lm2)
```

```{r}
# Report Confidence Intervals of x3-Estimates
confint(f.lm2)
```

###########

```{r, fig.height=4, fig.width=7}
# Model lm2 residuals
res_lm2 = residuals(f.lm2)

par(mfrow = c(1, 2))
acf(res_lm2, main = "Residuals ACF")
pacf(res_lm2, main = "Residuals PACF")
```

##################
Check ts Residuals
```{r}
checkresiduals(res_lm2, xlab = "Time (years)", ylab = "Residuals of y2-Model M5")
```

######################
PREDICTION NEXT HERE
######################
```{r, fig.height=5, fig.width=6}
# Predict with data points
x.val2 = round(seq(max(x3), min(x3), length = 28))
y.val2 = predict(f.lm2, data = data.frame(x3 = x.val2))

Pred.O_Tb = ts(y.val2, start = c(1995, 1), frequency = 1)
y.b2 = cbind(Pred.O_Tb, yy)
ts.plot(y.b2, ylim = c(0,230), 
        xlab = "Time (years)", ylab = "Cases of P-Tb", lwd = 2,
        main = expression(paste("AR(2) Model Predictions of Target", ~italic("y"[2]))), 
        col = c(5,7))
legend("topright", legend = c("M5 O-Tb", "M5 Pred O-Tb"), 
       lty = 1, lwd = 2, col = c(5,7))
```

```{r}
# Check MAE
mae = abs(Pred.O_Tb - yy)
er2 = mean(mae)
paste("The mean absolute error of the predictions = ", er2)
```

```{r}
paste("The mean square error of the predictions = ", mean((Pred.O_Tb - yy)^2))
```




```{r}
# PREDICTING NEW CASES
n.2023 = 55.589 + 0.8533*Pred.O_Tb[28] - 0.18*Pred.O_Tb[27] - 1.6431*29

n.2024 = 55.589 + 0.8533*n.2023 - 0.18*Pred.O_Tb[28] - 1.6431*30 
n.2025 = 55.589 + 0.8533*n.2024 - 0.18*n.2023 - 1.6431*31
n.2026 = 55.589 + 0.8533*n.2025 - 0.18*n.2024 - 1.6431*32

n.2027 = 55.589 + 0.8533*n.2026 - 0.18*n.2025 - 1.6431*33
n.2028 = 55.589 + 0.8533*n.2027 - 0.18*n.2026 - 1.6431*34
n.2029 = 55.589 + 0.8533*n.2028 - 0.18*n.2027 - 1.6431*35

cas = cbind(n.2023, n.2024, n.2025, n.2026)
rownames(cas) = c("O-Tb Cases")
print(cas)
```


##############
Load Libraries

```{r}
library(tseries)
library(forecast)
library(vars)
```

########################
   PART 2 of ANALYSIS   
########################

Here, we construct a data generating process (DGP) for the tuberculosis subtypes data, and once we identify the model to be adequate, we further use it to predict future observations of the Tb data.


```{r}
# Model selection order 
d.tb = cbind(x2, x3)
colnames(d.tb) = c("P-Tb", "O-Tb")

info.tb = VARselect(d.tb,lag.max = 5, type = "const")
info.tb$selection
```

##################
MODEL FORMULATION
##################
```{r}
# Model with constant/ trend terms
y1 = x2; y2 = x3
up_x.trtb = cbind(y1, y2)
modnew_m23 = VAR(up_x.trtb, p = 2, type = c("both"))
summary(modnew_m23)
```

#####################################
   Compute Crude CI95 of Estimates
#####################################
```{r}
### CI95 for VAR x2 Estimates
Estimate = c(0.006123, 0.938731, 0.405163,
             -0.459183, 127.236144, -2.937498)
Std.Error = c(0.182133, 0.294428, 0.170114,
              0.352624, 72.904573, 1.911785)
low.ci95 = c(Estimate + c(-1)*1.96*Std.Error)
up.ci95 = c(Estimate + c(1)*1.96*Std.Error)
x2.ci95 = cbind(Estimate, Std.Error, low.ci95, up.ci95)
print(x2.ci95)
```

```{r}
### CI95 for VAR x3 Estimates
Estimates = c(-0.05503, 0.92442, 0.29670,
             -0.43870, -5.97039, -0.12357)
Std.Errors = c(0.11566, 0.18698, 0.10803,
              0.22393, 46.29832, 1.21409)

lower.ci95 = c(Estimates + c(-1)*1.96*Std.Errors)
upp.ci95 = c(Estimates + c(1)*1.96*Std.Errors)
x3.ci95 = cbind(Estimates, Std.Errors, lower.ci95, upp.ci95)
print(x3.ci95)
```


```{r}
# Check the residual values
res_mnew = resid(modnew_m23)
head(res_mnew)
```

```{r}
#### Visualize residuals
checkresiduals(res_mnew[, 1], ylab = "Residuals of y1-Model M6")
```

```{r}
checkresiduals(res_mnew[, 2], ylab = "Residuals of y2-Model M7")
```

```{r, fig.height=5, fig.width=5}
# Cross correlation function of Residuals

par(mfrow = c(2,2))
ccf(res_mnew[, 1], res_mnew[, 1], main = "r11")
ccf(res_mnew[, 1], res_mnew[, 2], main = "r12")
ccf(res_mnew[, 2], res_mnew[, 1], main = "r21")
ccf(res_mnew[, 2], res_mnew[, 2], main = "r22")
```

```{r, fig.height=5, fig.width=5}
# Check the acf&pacf
acf(res_mnew)
pacf(res_mnew)
```

#######################
    More Diagnostics
#######################
```{r, fig.height=7, fig.width=5}
pl.stab = stability(modnew_m23, type = "OLS-CUSUM")
plot(pl.stab)
```

```{r}
# Check the roots
roots(modnew_m23)
```

```{r}
# Normality test
normality.test(modnew_m23)
```

```{r}
# Autocorrelation test
serial.test(modnew_m23, type = "BG")
```

```{r}
# AR conditional heteroscedasticity
arch.test(modnew_m23)
```

```{r}
library(vars)
# AR conditional heteroscedasticity
#arch.test(modnew_m23)
```


##############################
    Forecast with VAR MODEL
##############################
PART A:
```{r}
var.pred = predict(modnew_m23, start = c(1995, 1), 
                   n.ahead = 5, ci = 0.95,
                   dumvar = NULL, frequency = 1, 
                   interval = "prediction")
var.pred
```

```{r, fig.height=9, fig.width=5}
#par(mfrow =c(1,2))
plot(var.pred, xlab = "Time (years)", ylab = "Cases")
```


```{r, fig.height=5, fig.width=6}
# Refined graphs
var.ts = var.pred$endog
f = var.pred$fcst
fc = cbind(f$y1[,1], f$y2[,1])
fc.ts = ts(fc, start = c(2025,1))

ci.1 = rbind(c(129), f$y1[,2:3])
ci.2 = rbind(c(56), f$y2[,2:3])

tci1 = ts(ci.1, start = c(2024,1)) # Confidence interval 
tci2 = ts(ci.2, start = c(2024,1))


x.data = rbind(var.ts, fc.ts)
xp = seq(1995,2029,1)

# Future values of y1
plot(NULL, xlim = c(1995,2029), 
     ylim = c(0,450), 
     xlab = "Time (years)",
     ylab = "Cases",
     main = expression(paste("Forecast of Series", ~italic("y"[1]))))
lines(1995:2024,  x.data[1:30,1], type = "l")
lines(2024:2029, x.data[30:35,1], lty = 2, 
      col = "magenta", cex = 2)
lines(2024:2029, tci1[,1], lty = 2, col = "blue")
lines(2024:2029, tci1[,2], lty = 2, col = "blue")
abline(v = 2024, lty = 2, col = "red")
```

```{r, fig.height=5, fig.width=6}
# Future values of y2
plot(NULL, xlim = c(1995,2029), 
     ylim = c(-30,230), 
     xlab = "Time (years)",
     ylab = "Cases",
     main = expression(paste("Forecast of Series", ~italic("y"[2]))))
lines(1995:2024, x.data[1:30, 2], type = "l")
lines(2024:2029, x.data[30:35, 2], lty = 2, 
      col = "magenta", cex = 2)
lines(2024:2029, tci2[,1], lty = 2, col = "blue")
lines(2024:2029, tci2[,2], lty = 2, col = "blue")
abline(v = 2024, lty = 2, col = "red")
```


```{r, fig.height=7,fig.width=4}
# Axis
#plot(var.pred,  
     #xaxt = "n",
     #xlab = "Time (years)",
     #ylab = "Cases")
#axis(1, at = seq(1,35, by = 5)  ,
     #labels = c(1995,2000,2010,2015,2020,2025,2030))
```


##################################
   UPDATED VAR MODEL EVALUATION
##################################
Now, let us evaluate the performance of the updated VAR model in making predictions.
Here, we shall divide the data into training and testing sets. 80\% train data and 20\% 
test data.
```{r}
# Settings 
index = 1:21
tr.data = window(up_x.trtb, end = c(2018, 1))
te.data = window(up_x.trtb, start = c(2019, 1))
tail(tr.data)
```

##############
  Forecasts
##############
```{r}
# Forecasting 
var.fit = VAR(tr.data, p = 2, type = "both")
var.forc = predict(var.fit, n.ahead = 6, ci = 0.95,
                   dumvar = NULL, interval = "prediction") 
var.forc
```

```{r}
te.data[,1]
```

```{r}
te.data[,2]
```

###################
  Forecast for x3
###################
```{r, fig.height=5, fig.width=6}
plot(var.forc, names = "y2", 
     xlab = "Time (years)", ylab = "Cases")
points(c(25:30), te.data[, 2])
```

##################
   Alternative
##################
```{r}
library(forecast)

# Using forecast function
h = 6
var.for = forecast(var.fit, h) 
var.for
```

#####################
Check Model Accuracy
#####################

```{r}
# Forecast Accuracy
acc_for = accuracy(var.for, te.data, D = NULL, d = NULL)
print(acc_for)
```

###################
  Check More Here
###################

```{r}
# forecast for x2
df.x2 = as.data.frame(var.forc$fcst$y1)
df.x2 = df.x2[, 1]  
df.x2
```

```{r}
# forecast for x3
df.x3 = as.data.frame(var.forc$fcst$y2)
df.x3 = df.x3[, 1] 
df.x3
```

```{r}
# Forecast error x2
for_error2 = abs(df.x2 - te.data[, 1])
mae.x2 = mean(for_error2)
mae.x2
```

```{r}
# Error forecast x3
for_error3 = abs(df.x3 - te.data[, 2])
mae.x3 = mean(for_error3)
mae.x3
```

```{r, fig.height=9, fig.width=5}
# In sample prediction
library(tseries)
plot(var.forc, xlab = c("Time(years)"), 
     ylab = c("Cases")
     )
#points(te.data) 
```


```{r, fig.height=5, fig.width=6}
# Refine the x-axis
f1 = var.forc$fcst$y1 
f2 = var.forc$fcst$y2

f1.ts = ts(rbind(c(155), f1), start = c(2018,1))
f2.ts = ts(rbind(c(70), f2), start = c(2018,1))

xd.new = rbind(tr.data, cbind(f1.ts[-1,1], f2.ts[-1,1]))
xd.new2 = ts(xd.new, start = c(1995,1))
ci.f1 = f1.ts[,2:3] #CI 95
ci.f2 = f2.ts[,2:3]


t1 = ts(c(155,70), start = c(2018,1))
te.data1 = ts(rbind(t1, te.data), start = c(2018,1))
# Visualize y1
plot(NULL, xlim = c(1995,2025), 
     ylim = c(0,450),
     xlab = "Time (years)",
     ylab = "Cases",
     main = expression(paste("Out-Sample Forecast of Series", ~italic("y"[1]))))
points(2018:2024, te.data1[, 1], type = "l", col = "gold")
lines(1995:2018, xd.new2[1:24, 1], type = "l")

lines(2018:2024, xd.new2[24:30, 1], lty = 2, col = "magenta")
lines(2018:2024, ci.f1[,1], lty = 2, col = "blue")
lines(2018:2024, ci.f1[,2], lty = 2, col = "blue")
abline(v = 2018, lty = 2, col = "red")
```


```{r, fig.height=5, fig.width=6}
# In-sample F-y2

# Visualize y1
plot(NULL, xlim = c(1995,2025), ylim = c(-30,230),
     xlab = "Time (years)",
     ylab = "Cases",
     main = expression(paste("Out-Sample Forecast of Series",
                             ~italic("y"[2]))))
lines(1995:2018, xd.new2[1:24, 2], type = "l")
points(2018:2024, te.data1[, 2], 
       type = "l", col = "gold")

#xd = seq(1995,2024,1)
#plot.ts(xd.new2[, 2], type = "l", 
        #ylim = c(-30,230),
        #xlab = "Time (years)",
        #ylab = "Cases",
        #main = expression(paste("Out-Sample Forecast of Series", ~italic("y"[2]))))
#points(2019:2024, xd.new2[25:30, 2], col = 2)
lines(2018:2024, xd.new2[24:30, 2], lty = 2, col = "magenta")
lines(2018:2024, ci.f2[,1], lty = 2, col = "blue")
lines(2018:2024, ci.f2[,2], lty = 2, col = "blue")
abline(v = 2018, lty = 2, col = "red")

```



######################
#  Further Analysis  #
######################

First, let's get the out-sample forecasts augmented 
to each sub-type disease data.

```{r}
# Forecast 5yrs for X2-AR2
pred.ar2x2 = rbind(C25, C26, C27, C28, C29)
#rownames(pred.ar2x2) = c("2025", "2026", "2027",
                         #"2028", "2029")

colnames(pred.ar2x2) = c("fcst.X2")
print(as.numeric(pred.ar2x2))
```

```{r}
# Forecast 5yrs for X3-AR2
n.2029 = 55.589 + 0.8533*n.2028 - 0.18*n.2027 - 1.6431*35
pred.ar2x3 = rbind(n.2025, n.2026, n.2027, n.2028, n.2029)
#rownames(pred.ar2x3) = c("2025", "2026", "2027",
                         #"2028", "2029")

#colnames(pred.ar2x3) = c("fcst.X3")
print(as.numeric(pred.ar2x3))
```

```{r}
# Univariate forecast ts y1
f.AR2 = cbind(pred.ar2x2, pred.ar2x3)
ar2.f = ts(f.AR2, start = c(2025, 1), frequency = 1)
colnames(ar2.f) = NULL
ar2.x2 = ar2.f[,1]
ar2.x2
```

```{r}
# Univariate forecast ts y2
ar2.x3 = ar2.f[,2]
ar2.x3
```

```{r}
# Forecasts ts for x2/y1 VAR
fvarX2 = ts(var.pred$fcst$y1, start = c(2025, 1))
fx2 = fvarX2[,1]
colnames(fx2) = NULL
fx2
```

```{r}
# Forecasts ts for x3/y2 VAR
fvarX3 = ts(var.pred$fcst$y2, start = c(2025, 1))
fx3 = fvarX3[,1]
colnames(fx3) = NULL
fx3
```

```{r}
X.f = cbind(ar2.f, fx2, fx3)
colnames(X.f) = NULL
X.f
```

```{r, fig.width=7}
# Visualize the prognosis
par(mar = c(4,4,4,9), xpd = T)
ts.plot(X.f, 
        xlab = "Time (years)",
        ylab = "Forecasts Values",
        main = "Forecasts Overview",
        col = 2:5, lwd = 1.5, lty = 1)
legend("topright", inset = c(-0.4, 0.2), title = "Models",
       legend = c(expression(paste("y"[1],AR2)), 
                  expression(paste("y"[2],AR2)), 
                  expression(paste("y"[1],VAR)), 
                  expression(paste("y"[2],VAR))),
       col = 2:5, lwd = 1.5)
grid()
```

```{r}
# Alternative
# Combining the Forecasts methods
```


```{r, fig.height=5, fig.width=4}
library(forecast)
#autoplot(fx.VAR2, xlab("Time (years)")) 
```



```{r}
Cx2 = (X.f[,1] + X.f[,3]) / 2
Cx3 = (X.f[,2] + X.f[,4]) / 2
Cd = cbind(Cx2, Cx3)
colnames(Cd) = NULL
print(Cd)
```

```{r, fig.width=7}
# Visuals
par(mar = c(4,4,4,9), xpd = T)
X.f2 = cbind(X.f, Cd)
ts.plot(X.f2, 
        xlab = "Time (years)",
        ylab = "Forecasts Values",
        main = "Forecasts Overview",
        col = 2:7, lwd = 1.5, type = "l")
legend("topright", inset = c(-0.4, 0.2), title = "Models",
       #legend = c("y1.AR2", "y2.AR2", "y1.VAR", "y2.VAR",
                  #"Cy1", "Cy2"),
       legend = c(expression(paste("y"[1],AR2)), 
                  expression(paste("y"[2],AR2)), 
                  expression(paste("y"[1],VAR)), 
                  expression(paste("y"[2],VAR)),
                  expression(paste(italic(C),"y"[1])), 
                  expression(paste(italic(C),"y"[2]))),
       col = 2:7, lwd = 1.5)
grid()
```



```{r}
library(ggplot2)

# Visualize the combined forecast
autoplot(X1[,1], ylab = "P-Tb Cases") + 
  autolayer(Cx2, series = "Combined") + 
  autolayer(X.f2[,1], series = "X2 AR2") +
  autolayer(X.f2[,3], series = "X2 VAR2") +
  ggtitle("P-Tb target disease forecasts")
  
```


########################
# VAR CROSS VALIDATION #
########################

```{r}
# CV Implementation
dat.Y = up_x.trtb

k = 24
Tt = 30
test.D = dat.Y[(Tt-(Tt-k-1)):Tt, ]

# Model 1 forecast
i = 1
test_Model1 = VAR(dat.Y[1:(k + i - 1), ], 
                   p = 2,
                   type = "both")
f1 = predict(test_Model1, n.ahead = 1, 
               ci = 0.95, dumvar = NULL, 
               interval = "prediction")

y1f1 = f1$fcst$y1[1:4]
y2f1 = f1$fcst$y2[1:4]
```



```{r}
# Model 2 forecast
i = 2
test_Model2 = VAR(dat.Y[1:(k + i - 1), ], 
                   p = 2,
                   type = "both")
f2 = predict(test_Model2, n.ahead = 1, 
               ci = 0.95, dumvar = NULL, 
               interval = "prediction")

y1f2 = f2$fcst$y1[1:4]
y2f2 = f2$fcst$y2[1:4]
```



```{r}
# Model 3 forecast
i = 3
test_Model3 = VAR(dat.Y[1:(k + i - 1), ], 
                   p = 2,
                   type = "both")
f3 = predict(test_Model3, n.ahead = 1, 
               ci = 0.95, dumvar = NULL, 
               interval = "prediction")
y1f3 = f3$fcst$y1[1:4]
y2f3 = f3$fcst$y2[1:4]
```



```{r}
# Model 4 forecast
i = 4
test_Model4 = VAR(dat.Y[1:(k + i - 1), ], 
                   p = 2,
                   type = "both")
f4 = predict(test_Model4, n.ahead = 1, 
               ci = 0.95, dumvar = NULL, 
               interval = "prediction")

y1f4 = f4$fcst$y1[1:4]
y2f4 = f4$fcst$y2[1:4]
```



```{r}
# Model 5 forecast
i = 5
test_Model5 = VAR(dat.Y[1:(k + i - 1), ], 
                   p = 2,
                   type = "both")
f5 = predict(test_Model5, n.ahead = 1, 
               ci = 0.95, dumvar = NULL, 
               interval = "prediction")

y1f5 = f5$fcst$y1[1:4]
y2f5 = f5$fcst$y2[1:4]
```



```{r}
# Model 6 forecast
i = 6
test_Model6 = VAR(dat.Y[1:(k + i - 1), ], 
                   p = 2,
                   type = "both")
f6 = predict(test_Model6, n.ahead = 1, 
               ci = 0.95, dumvar = NULL, 
               interval = "prediction")

y1f6 = f6$fcst$y1[1:4]
y2f6 = f6$fcst$y2[1:4]
```



```{r}
# Combine Results in TS 
y1fcas.vals = ts(rbind(y1f1, y1f2, y1f3, y1f4, y1f5, y1f6), 
                 start = c(2019,1))
cnames_cv1 = c("fcst", "lower", "upper", "CI")
colnames(y1fcas.vals) = cnames_cv1


y2fcas.vals = ts(rbind(y2f1, y2f2, y2f3, y2f4, y2f5, y2f6),
               start = c(2019,1))
```


```{r}
# Compute Error Measure
y1.forc.err = y1fcas.vals[, 1] - test.D[, 1]
mae.fy1 = mean(abs(y1.forc.err))
print(paste("The y1 insample forecast error in terms of the MAE quantity = ", mae.fy1))
```

```{r}
y2.forc.err = y2fcas.vals[, 1] - test.D[, 2]
mae.fy2 = mean(abs(y2.forc.err))
print(paste("The y2 insample forecast error in terms of the MAE quantity = ", mae.fy2))
```

```{r, fig.height=5, fig.width=6}
fore1 = ts(rbind(rep(c(155),4), y1fcas.vals))

# Plotting the results
plot(NULL, xlim = c(1995,2024), ylim = c(0, 450),
     xlab = "Time (years)", ylab = "Cases",
     main = expression(paste("Holdout-sample Forecast of ", "y"[1])))
lines(ts(dat.Y[1:24, 1], start = c(1995,1)), lty = 1)
points(2018:2024, te.data1[, 1], 
       type = "l", col = "gold")
lines(ts(fore1[,1], start = c(2018,1)), 
      lty = 2, col = "magenta")
lines(ts(fore1[,2], start = c(2018,1)), 
      lty = 2, col = "blue")
lines(ts(fore1[,3], start = c(2018,1)), 
      lty = 2, col = "blue")
abline(v = c(2018), lty = 2, col = "red")
```



#######################
      NEXT IS y2
#######################

```{r, fig.height=5, fig.width=6}
fore2 = ts(rbind(rep(c(70),4), y2fcas.vals))

# Plotting the results
plot(NULL, xlim = c(1995,2024), ylim = c(0, 230),
     xlab = "Time (years)", ylab = "Cases",
     main = expression(paste("Holdout-sample Forecast of ", "y"[2])))
lines(ts(dat.Y[1:24, 2], start = c(1995,1)), lty = 1)
points(2018:2024, te.data1[,2],
       type = "l", col = "gold")
lines(ts(fore2[,1], start = c(2018,1)), 
      lty = 2, col = "magenta")
lines(ts(fore2[,2], start = c(2018,1)), 
      lty = 2, col = "blue")
lines(ts(fore2[,3], start = c(2018,1)), 
      lty = 2, col = "blue")
abline(v = c(2018), lty = 2, col = "red")

```


```{r}
# Computing the Error Measures for y1
# 1. ME
e.y1 = y1fcas.vals[,1] - test.D[,1]
me.y1 = mean(e.y1)

# 2. MAE
mae.1y = mean(abs(e.y1))

# 3. RMSE 
y1.rmse = sqrt(mean(e.y1^2))

# 4. MASE
N = 24
j = 2:N
Y = dat.Y[1:N, 1]
Q.scale = (1/(N-1))*sum(abs(Y[j]-Y[j-1]))

y1.mase = mae.1y / Q.scale

# 5. MPE
n = 6
y1.mpe = (1/n)*100*(sum((e.y1/test.D[,1])))

# 6. MAPE
y1.mape = (1/n)*100*(sum(abs(e.y1/test.D[,1])))

y1.Errors = cbind(me.y1, y1.rmse, mae.1y, y1.mpe, y1.mape, y1.mase)
co.names = c("ME", "RMSE", "MAE", "MPE", "MAPE", "MASE")
colnames(y1.Errors) = co.names
rownames(y1.Errors) = "y1 Test Sets"
y1.Errors
```


```{r}

# 1. ME
e.y2 = y2fcas.vals[,1] - test.D[,2]
me.y2 = mean(e.y2)

# 2. MAE
mae.2y = mean(abs(e.y2))

# 3. RMSE 
y2.rmse = sqrt(mean(e.y2^2))

# 4. MASE
N = 24
j = 2:N
yY = dat.Y[1:N, 2]
Q.scale2 = (1/(N-1))*sum(abs(yY[j]-yY[j-1]))

y2.mase = mae.2y / Q.scale2

# 5. MPE
n = 6
y2.mpe = (1/n)*100*(sum((e.y2/test.D[,2])))

# 6. MAPE
y2.mape = (1/n)*100*(sum(abs(e.y2/test.D[,2])))

y2.Errors = cbind(me.y2, y2.rmse, mae.2y, y2.mpe, y2.mape, y2.mase)
co.names = c("ME", "RMSE", "MAE", "MPE", "MAPE", "MASE")
colnames(y2.Errors) = co.names
rownames(y2.Errors) = "y2 Test Sets"
y2.Errors

```

```{r}
# Combine results
Errors.T = rbind(y1.Errors, y2.Errors)
print(Errors.T)
```

###############
CV METHOD 2
###############

```{r}
#h = 1
#k = 21
#Tt = 30
#i = 1:(Tt-k-h+1)
#i = 1
#test.D = dat.Y[(Tt-(Tt - k - h - i + 1)):Tt, ]

# Model 1 forecast
#test_Model = c()
#fcs = c()

#teMod.1 = VAR(dat.Y[1:(k + i - 1), ], 
                   #p = 2, type = "both")

#fcs.1 = predict(teMod.1, n.ahead = 3, 
                #ci = 0.95, dumvar = NULL, 
               #interval = "prediction")
```



```{r}
# 
#test_Mod = c()
#for (i in 1:(Tt-k-h+1)) {
  #test_Mod[[i]] = VAR(dat.Y[1:(k + i - 1), ], 
  #p = 2, type = "both")
  #return(test_Model)
#}
#test_Mod[[1]]
```







```{r}
# By Prof. Hyndman Approach
library(fpp3)
y.1 = dat.Y[,1]
y.2 = dat.Y[,2]
dat.y = cbind(y.1,y.2) |>
  as_tsibble(pivot_longer = FALSE)
dat.y_stretch <- dat.y |>
  stretch_tsibble(.init = 21, .step = 1)
#dat.y_stretch %>%
  #stretch_tsibble()

fit.VAR = dat.y_stretch %>%
  model(VAR(vars(y.1,y.2) ~ AR(2)))

fit.VAR |>
  forecast(h = 1) |>
  accuracy(dat.y)
```

######################
STLL ON VAR EVALUATION
######################

```{r}

```


```{r}
#est.dat = c()
#for (j in 1:9) {
  #est.dat[[j]] = ts(matrix(0, nrow = 21+j-1, ncol = 2), start = c(1995, 1))
  #return(matrix(unlist(est.dat[[j]]), ncol = 3))
#}

```

# Expanding Window Idea ####

```{r}
# 1. 
y1 = x2
y2 = x3
up_x.trtb2 = up_x.trtb

D1 = window(up_x.trtb2, end = c(2015, 1))
D2 = window(up_x.trtb2, end = c(2016, 1))
D3 = window(up_x.trtb2, end = c(2017, 1))
D4 = window(up_x.trtb2, end = c(2018, 1))
D5 = window(up_x.trtb2, end = c(2019, 1))

D6 = window(up_x.trtb2, end = c(2020, 1))
D7 = window(up_x.trtb2, end = c(2021, 1))
D8 = window(up_x.trtb2, end = c(2022, 1))
D9 = window(up_x.trtb2, end = c(2023, 1))

```



```{r}
# Models
library(vars)
M1 = vars::VAR(D1, p = 2, type = "both")
M2 = vars::VAR(D2, p = 2, type = "both")
M3 = vars::VAR(D3, p = 2, type = "both")
M4 = vars::VAR(D4, p = 2, type = "both")
M5 = vars::VAR(D5, p = 2, type = "both")
M6 = vars::VAR(D6, p = 2, type = "both")
M7 = vars::VAR(D7, p = 2, type = "both")
M8 = vars::VAR(D8, p = 2, type = "both")
M9 = vars::VAR(D9, p = 2, type = "both")
```



```{r}
library(vars)
# forecasts
fcs.E1 = predict(M1, n.ahead = 9,
                 ci = 0.95, 
              dumvar = NULL)
fcs.E2 = predict(M2, n.ahead = 8,
                 ci = 0.95, 
                dumvar = NULL)
fcs.E3 = predict(M3, n.ahead = 7,
                 ci = 0.95, 
                dumvar = NULL)
fcs.E4 = predict(M4, n.ahead = 6,
                 ci = 0.95, 
                dumvar = NULL)
fcs.E5 = predict(M5, n.ahead = 5,
                 ci = 0.95, 
                dumvar = NULL)
fcs.E6 = predict(M6, n.ahead = 4,
                 ci = 0.95, 
                dumvar = NULL)
fcs.E7 = predict(M7, n.ahead = 3,
                 ci = 0.95, 
                dumvar = NULL)
fcs.E8 = predict(M8, n.ahead = 2,
                 ci = 0.95, 
                dumvar = NULL)
fcs.E9 = predict(M9, n.ahead = 1,
                 ci = 0.95, 
                dumvar = NULL)
```



```{r}
fcs.E1$fcst$y1[9,1]
fcs.E1$fcst$y2[9,1]
fcs.E2$fcst$y1[8,1]
fcs.E2$fcst$y2[8,1]
fcs.E3$fcst$y1[7,1]
fcs.E3$fcst$y2[7,1]
fcs.E4$fcst$y1[6,1]
fcs.E4$fcst$y2[6,1]
fcs.E5$fcst$y1[5,1]
fcs.E5$fcst$y2[5,1]
fcs.E6$fcst$y1[4,1]
fcs.E6$fcst$y2[4,1]
fcs.E7$fcst$y1[3,1]
fcs.E7$fcst$y2[3,1]
fcs.E8$fcst$y1[2,1]
fcs.E8$fcst$y2[2,1]
fcs.E9$fcst$y1[1,1]
fcs.E9$fcst$y2[1,1]
```









```{r}
# Expanding Window
library(zoo)
library(forecast)
```






















######################
  ADF Unit Root Test
######################

```{r}
# For x2 target disease
#adf_t.x2 = adf.test(x2, alternative = "stationary")
#print(adf_t.x2)
```

```{r}
# For x3 target disease
#adf_t.x3 = adf.test(x3, alternative = "stationary")
#print(adf_t.x3)
```

###################
  KPSS Trend Test
###################

```{r}
# Test for x2
#kpss.x2 = kpss.test(x2, null = "Trend")
#print(kpss.x2)
```

```{r}
# Test for x3
#kpss.x3 = kpss.test(x3, null = "Trend")
#print(kpss.x3)
```

#####################
    MORE ANALYSIS
#####################

```{r, fig.height=7,fig.width=5}
library(urca)
y1.pp.ts = ur.pp(x2, type = "Z-tau", 
                 model = "trend", lags = "short")
plot(y1.pp.ts)
```

```{r}
summary(y1.pp.ts)
```

```{r}

```


