
# cointegration test

coint.test(y = Y, X = X[,2],d=4)

fit1<-arima(Y,xreg = scale(X[,-1]),include.mean = F)
fit1
library(aTSA)
adf.test(fit1$residuals)

library(aTSA)
test1 <- 
  for (i in 3:6) {
    adf.test(diff(log(unlist(data_raw[,i]))))
    print(Box.test(diff(log(unlist(data_raw[,i])))))
  }
