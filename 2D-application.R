########################################################################

source("fit_2drarma.r")

library(R.matlab)
library(raster)

im_1 = readMat("m_1_p_1.mat") # read the image
im_1 = im_1$im
im = raster(im_1)
plot(im, col = gray.colors(255))

n1 = k1 = 80

y = im_1[1346:(1346+n1), 920:(920+k1)] # selected region
plot(raster(y), col = gray.colors(255))

fit1 = fit.2drarma(y, 1, 1)
fit1$model
fit2 = fit.2drarma(y, 1, 2)
fit2$model
fit3 = fit.2drarma(y, 1, NA)
fit3$model

