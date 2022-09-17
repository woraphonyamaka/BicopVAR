# Hello, world!# This is an  function named 'Bivariate copula based VAR model'
# Written by Dr.Worphon Yamaka,
# Center of excellence in Econometrics, Faculty of Economics,
# Chaing Mai University
# Impulse response function

IRF=function(Y,A,Sigma,h.step) {
  Binverse=t((chol(Sigma)))
  m=ncol(Y)
  IRFfit=matrix(c(NA),h.step,4)
  IRFfit[1,]=c(Binverse)
  for( i in 2:h.step){
    IRFfit[i,]=(A^i)%*%Binverse
  }

  par(mfrow=c(2,2))
  plot(ts(IRFfit[,1]), ylab="response of y1", main="response of y1 to the shock of y1", lwd=2, lty=2)
  abline(h=0, col="red",lwd=2)
  plot(ts(IRFfit[,2]), ylab="response of y2",main="response of y2 to the shock of y1", lwd=2, lty=2)
  abline(h=0, col="red",lwd=2)
  plot(ts(IRFfit[,3]), ylab="response of y1 ",main="response of y1 to the shock of y2", lwd=2, lty=2)
  abline(h=0, col="red",lwd=2)
  plot(ts(IRFfit[,4]), ylab="response of y2 ", main="response of y2 to the shock of y2",lwd=2, lty=2)
  abline(h=0, col="red",lwd=2)
  IRFfit
}
