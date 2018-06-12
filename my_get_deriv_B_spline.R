library(fda)
library(R.matlab)

setwd('/home/cjm/Desktop/Summer/Data')
data = readMat("theta_test_v.mat")
#x = data$theta[1, ]*4
x = data$theta.test.v
# can get breaks from knt2brk function in matlab
breaks = c(0,1/2,1)*pi

bS = bsplineS(x, breaks, nderiv=1)
writeMat('deriv_B_spline_test_v2.mat', bS=bS)
