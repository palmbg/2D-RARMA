########################################################################
source("simu_2drarma.r")
source("fit_2drarma.r")


alpha =  0.5052 # intercept (beta)
phi = c(-0.1708,0.2035,0.5689) # AR parameters
theta = c(0.4804,-0.0151,0.4270) # MA parameters

y = simu.2drarma(alpha, phi, theta, p = 1, q = 1, n = 80 , k = 80) # ARMA simulation


library(raster)
plot(raster(y), col = gray.colors(255)) 
image(y)
graphics::image(t(apply(y, 2, rev)),axes=FALSE, col=gray(1:100/100))

fit1 = fit.2drarma(y, 1, 1); fit1$model # adjust
fit2 = fit.2drarma(y, 2, 2); fit2$model # adjust
fit3 = fit.2drarma(y, 1, 2); fit3$model # adjust

y = simu.2drarma(alpha, NA, theta, p = NA, q = 1, n = 10 , k = 10) # MA simu
fit4 = fit.2drarma(y, 1, NA); fit4$model
fit5 = fit.2drarma(y, 2, NA); fit5$model

y = simu.2drarma(alpha, phi, NA, p = 1, q = NA, n = 10 , k = 10) # AR simu
fit6 = fit.2drarma(y, NA, 1); fit6$model
fit7 = fit.2drarma(y, NA, 2); fit7$model

