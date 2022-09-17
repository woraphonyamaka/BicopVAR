# Hello, world!# This is an  function named 'Bivariate copula based VAR model'
# Written by Dr.Worphon Yamaka,
# Center of excellence in Econometrics, Faculty of Economics,
# Chaing Mai University
# Main function


Bicopmix=function(Y, p,dist1,dist2,copula){
  ##------MixCOP VAR------------------
  cop2=function(param,family){

    if (family=="normal"){
      rho1=param
      c = normalCopula(c(rho1), dim = 2)
    }
    if (family=="studentt"){
      rho1=param[1]
      df=param[2]
      c = tCopula(c(rho1),df=df, dim=2)
    }

    if (family=="frank"){
      c = archmCopula("frank", param= param, dim = 2)
    }


    if (family=="joe"){
      c = archmCopula("joe", param= param, dim = 2)
    }
    if (family=="gumbel"){
      c = archmCopula("gumbel", param= param, dim = 2)
    }
    if (family=="clayton"){
      c = archmCopula("clayton", param= param, dim = 2)
    }

    c
  }

  loglikmixCopula2 <- function(param, u, family1, family2) {

    K=length(param)
    w11=param[K]
    w12=1-param[K]

    if (family1=="studentt"|family2=="studentt"){

      mC <- mixCopula(list(cop(param=param[1], family=family1),
                           cop(param=c(param[2],param[3]), family=family2)),
                      c(w11,w12))
    }else{

      mC <- mixCopula(list(cop(param=param[1], family=family1),
                           cop(param=param[2], family=family2)),
                      c(w11,w12))
    }
    cop.param <- getTheta(mC , freeOnly = TRUE, attr = TRUE)
    lower <- attr(cop.param, "param.lowbnd")
    upper <- attr(cop.param, "param.upbnd")

    ## FIXME-JY: param range check is only a *part* of validity check
    d.uM=  (dCopula(u, copula=mC , log=FALSE, checkPar=FALSE))+0.0000000001

    LL=sum(log(d.uM))
    #  cat("Sum of log Likelihood for MS-mixcop Model ->",sprintf("%4.4f",c(LL)),"\n")
    return(LL)
  }

  ###================================================
  copvar2 =function(param,Y, p,dist1,dist2,copula){
    n<-nrow(Y);	 	# # of observations in data set
    m<-ncol(Y);	 	# # of variables in data set

    if (dist1=="norm"){
      coef1=param[1:(1+(m*p))]
    }
    if (dist1=="std"){
      coef1=param[1:(1+1+(m*p))]
    }

    if (dist2=="norm"){
      coef2=param[((length(coef1)+1)):(length(coef1)+1+(m*p))]
    }
    if (dist2=="std"){
      coef2=param[(length(coef1)+1):(length(coef1)+1+1+(m*p))]
    }



    KK=length(coef1)+length(coef2)
    sigma11=param[KK+1]
    sigma12=param[KK+2]


    var.beta0 <- c(coef1[1],coef2[1]) # intercepts

    var.beta1 <- matrix(c(coef1[2:(1+m*p)],
                          coef2[2:(1+m*p)]), m, byrow=TRUE)


    y=Y
    season = NULL
    exogen = NULL
    lag.max = NULL
    y <- as.matrix(y)
    if (any(is.na(y)))
      stop("\nNAs in y.\n")
    if (ncol(y) < 2)
      stop("The matrix 'y' should contain at least two variables. For univariate analysis consider ar() and arima() in package stats.\n")
    if (is.null(colnames(y))) {
      colnames(y) <- paste("y", 1:ncol(y), sep = "")
      warning(paste("No column names supplied in y, using:",
                    paste(colnames(y), collapse = ", "), ", instead.\n"))
    }
    colnames(y) <- make.names(colnames(y))
    y.orig <- y
    obs <- dim(y)[1]
    K <- dim(y)[2]

    sample <- obs - p
    ylags <- embed(y, dimension = p + 1)[, -(1:K)]
    temp1 <- NULL
    for (i in 1:p) {
      temp <- paste(colnames(y), ".l", i, sep = "")
      temp1 <- c(temp1, temp)
    }
    colnames(ylags) <- temp1
    yend <- y[-c(1:p), ]

    rhs <- cbind(ylags, rep(1, sample))
    colnames(rhs) <- c(colnames(ylags), "const")


    if (!(is.null(season))) {
      season <- abs(as.integer(season))
      dum <- (diag(season) - 1/season)[, -season]
      dums <- dum
      while (nrow(dums) < obs) {
        dums <- rbind(dums, dum)
      }
      dums <- dums[1:obs, ]
      colnames(dums) <- paste("sd", 1:ncol(dums), sep = "")
      rhs <- cbind(rhs, dums[-c(1:p), ])
    }
    if (!(is.null(exogen))) {
      exogen <- as.matrix(exogen)
      if (!identical(nrow(exogen), nrow(y))) {
        stop("\nDifferent row size of y and exogen.\n")
      }
      if (is.null(colnames(exogen))) {
        colnames(exogen) <- paste("exo", 1:ncol(exogen),
                                  sep = "")
        warning(paste("No column names supplied in exogen, using:",
                      paste(colnames(exogen), collapse = ", "), ", instead.\n"))
      }
      colnames(exogen) <- make.names(colnames(exogen))
      tmp <- colnames(rhs)
      rhs <- cbind(rhs, exogen[-c(1:p), ])
      colnames(rhs) <- c(tmp, colnames(exogen))
    }
    datamat <- as.data.frame(rhs)
    colnames(datamat) <- colnames(rhs)

    Xlag=as.matrix(datamat[,c(1:(K*p))])

    cbind(yend,Xlag)
    e1=t(t(yend)-var.beta0 - tcrossprod(var.beta1, Xlag))

    if (dist1=="norm"){
      like1=dnorm(e1[,1], mean=0,sd = sigma11,log=FALSE) +0.00000001
      mar1=pnorm(e1[,1], mean=0,sd = sigma11)
    }
    if (dist2=="norm"){
      like2=dnorm(e1[,2], mean=0,sd = sigma12,log=FALSE) +0.00000001
      mar2=pnorm(e1[,2], mean=0,sd = sigma12)
    }


    if (dist1=="std"){
      df=coef1[length(coef1)]
      like1=dstd(e1[,1], mean=0,sd = sigma11,nu = df,log=FALSE) +0.00000001
      mar1=pstd(e1[,1], mean=0,sd = sigma11,nu = df)
    }
    if (dist2=="std"){
      df=coef2[length(coef2)]
      like2=dstd(e1[,2], mean=0,sd = sigma12,nu = df,log=FALSE) +0.00000001
      mar2=pstd(e1[,2], mean=0,sd = sigma12,nu = df)
    }


    U=cbind(mar1,mar2)-0.000000001

    if (copula=="normal"){
      rho1=param[(KK+3)]
      cop1=normalCopula(rho1, dim=m,dispstr = "un")}
    if (copula=="frank"){
      rho1=param[(KK+3)]
      cop1 <- archmCopula("frank", rho1, dim = m)}
    if (copula=="joe"){
      rho1=param[(KK+3)]
      cop1 <- archmCopula("joe",rho1, dim = m)}
    if (copula=="clayton"){
      rho1=param[(KK+3)]
      cop1 <- archmCopula("clayton",rho1, dim = m)}
    if (copula=="gumbel"){
      rho1=param[(KK+3)]
      cop1 <- archmCopula("gumbel",rho1, dim = m)}
    if (copula=="student"){
      rho1=param[(KK+3):(KK+4)]
      cop1 <- tCopula(rho1[1], dim=m,df = rho1[2],dispstr = "un") }


    if (copula=="NS"){
      rho1=param[(KK+3):(KK+5)]
      cop1 <- normalCopula(rho1[1], dim=m,dispstr = "ex")
      cop2 <- tCopula(rho1[2], dim=m,df = 5,dispstr = "ex")
      w=rho1[3]
      copmix =(w*(dCopula(U,cop1,log=FALSE))+((1-w)*(dCopula(U,cop2,log=FALSE))))+0.00000000000001
    }
    if (copula=="NC"){
      rho1=param[(KK+3):(KK+5)]
      cop1 <- normalCopula(rho1[1], dim=m,dispstr = "ex")
      cop2 <- archmCopula("clayton",rho1[2], dim = m)
      w=rho1[3]
      copmix =(w*(dCopula(U,cop1,log=FALSE))+((1-w)*(dCopula(U,cop2,log=FALSE))))+0.00000000000001
    }
    if (copula=="NF"){
      rho1=param[(KK+3):(KK+5)]
      cop1 <- normalCopula(rho1[1], dim=m,dispstr = "ex")
      cop2 <- archmCopula("frank",rho1[2], dim = m)
      w=rho1[3]
      copmix =(w*(dCopula(U,cop1,log=FALSE))+((1-w)*(dCopula(U,cop2,log=FALSE))))+0.00000000000001
    }
    if (copula=="NG"){
      rho1=param[(KK+3):(KK+5)]
      cop1 <- normalCopula(rho1[1], dim=m,dispstr = "ex")
      cop2 <- archmCopula("gumbel",rho1[2], dim = m)
      w=rho1[3]
      copmix =(w*(dCopula(U,cop1,log=FALSE))+((1-w)*(dCopula(U,cop2,log=FALSE))))+0.00000000000001
    }
    if (copula=="NJ"){
      rho1=param[(KK+3):(KK+5)]
      cop1 <- normalCopula(rho1[1], dim=m,dispstr = "ex")
      cop2 <- archmCopula("joe",rho1[2], dim = m)
      w=rho1[3]
      copmix =(w*(dCopula(U,cop1,log=FALSE))+((1-w)*(dCopula(U,cop2,log=FALSE))))+0.00000000000001
    }
    if (copula=="SC"){
      rho1=param[(KK+3):(KK+5)]
      cop1 <- tCopula(rho1[1], dim=m,df = 5,dispstr = "ex")
      cop2 <- archmCopula("clayton",rho1[2], dim = m)
      w=rho1[3]
      copmix =(w*(dCopula(U,cop1,log=FALSE))+((1-w)*(dCopula(U,cop2,log=FALSE))))+0.00000000000001
    }
    if (copula=="SF"){
      rho1=param[(KK+3):(KK+5)]
      cop1 <- tCopula(rho1[1], dim=m,df = 5,dispstr = "ex")
      cop2 <- archmCopula("frank",rho1[2], dim = m)
      w=rho1[3]
      copmix =(w*(dCopula(U,cop1,log=FALSE))+((1-w)*(dCopula(U,cop2,log=FALSE))))+0.00000000000001
    }
    if (copula=="SG"){
      rho1=param[(KK+3):(KK+5)]
      cop1 <- tCopula(rho1[1], dim=m,df = 5,dispstr = "ex")
      cop2 <- archmCopula("gumbel",rho1[2], dim = m)
      w=rho1[3]
      copmix =(w*(dCopula(U,cop1,log=FALSE))+((1-w)*(dCopula(U,cop2,log=FALSE))))+0.00000000000001
    }
    if (copula=="SJ"){
      rho1=param[(KK+3):(KK+5)]
      cop1 <- tCopula(rho1[1], dim=m,df = 5,dispstr = "ex")
      cop2 <- archmCopula("joe",rho1[2], dim = m)
      w=rho1[3]
      copmix =(w*(dCopula(U,cop1,log=FALSE))+((1-w)*(dCopula(U,cop2,log=FALSE))))+0.00000000000001
    }
    if (copula=="CF"){
      rho1=param[(KK+4):(KK+5)]
      cop1 <- archmCopula("frank",rho1[1], dim = m)
      cop2 <- archmCopula("clayton",rho1[2], dim = m)
      w=rho1[3]
      copmix =(w*(dCopula(U,cop1,log=FALSE))+((1-w)*(dCopula(U,cop2,log=FALSE))))+0.00000000000001
    }
    if (copula=="CG"){
      rho1=param[(KK+3):(KK+5)]
      cop1 <- archmCopula("gumbel",rho1[1], dim = m)
      cop2 <- archmCopula("clayton",rho1[2], dim = m)
      w=rho1[3]
      copmix =(w*(dCopula(U,cop1,log=FALSE))+((1-w)*(dCopula(U,cop2,log=FALSE))))+0.00000000000001
    }
    if (copula=="CJ"){
      rho1=param[(KK+3):(KK+5)]
      cop1 <- archmCopula("joe",rho1[1], dim = m)
      cop2 <- archmCopula("clayton",rho1[2], dim = m)
      w=rho1[3]
      copmix =(w*(dCopula(U,cop1,log=FALSE))+((1-w)*(dCopula(U,cop2,log=FALSE))))+0.00000000000001
    }
    if (copula=="FG"){
      rho1=param[(KK+3):(KK+5)]
      cop1 <- archmCopula("frank", rho1[1], dim = m)
      cop2 <- archmCopula("gumbel",rho1[2], dim = m)
      w=rho1[3]
      copmix =(w*(dCopula(U,cop1,log=FALSE))+((1-w)*(dCopula(U,cop2,log=FALSE))))+0.00000000000001
    }
    if (copula=="FJ"){
      rho1=param[(KK+3):(KK+5)]
      cop1 <- archmCopula("frank", rho1[1], dim = m)
      cop2 <- archmCopula("joe",rho1[2], dim = m)
      w=rho1[3]
      copmix =(w*(dCopula(U,cop1,log=FALSE))+((1-w)*(dCopula(U,cop2,log=FALSE))))+0.00000000000001
    }
    if (copula=="GJ"){
      rho1=param[(KK+3):(KK+5)]
      cop1 <- archmCopula("gumbel",rho1[1], dim = m)
      cop2 <- archmCopula("joe",rho1[2], dim = m)
      w=rho1[3]
      copmix =(w*(dCopula(U,cop1,log=FALSE))+((1-w)*(dCopula(U,cop2,log=FALSE))))+0.00000000000001
    }


    if (copula=="FG"|copula=="GJ"|copula=="SJ"|copula=="SG"|copula=="SF"|
        copula=="SC"|copula=="NJ"|copula=="NG"|copula=="NF"|copula=="NC"|
        copula=="CG"|copula=="NS"|copula=="FJ"|copula=="CJ"|copula=="CF"){
      likelihood=log(like1)+log(like2)+log(copmix)
    }else{
      copp1= dCopula(U,cop1,log=FALSE)+0.00000001
      likelihood=log(like1)+log(like2)+log(copp1)
    }

    LL=sum(likelihood)
    if (is.infinite(LL))  # control for optimization
      LL<--1000*100
    if (is.nan(LL))  # control for optimization
      LL<--10000*100
    if (is.na(LL))  # control for optimization
      LL<--10000*100
    cat("Sum of log Likelihood for COP-VAR Model ->",sprintf("%4.4f",LL),"\n")
    return(LL)

  }



  ############ START HERE ##################
  m=ncol(Y)
  # two dimension copula model ==> 1 correlations
  if (copula=="normal"){
    coppar=rep(0.5,1)
    cop1=normalCopula(rep(0.5,1), dim=m,dispstr = "un")}
  if (copula=="frank"){
    coppar=2
    cop1 <- archmCopula("frank", coppar, dim = m)}
  if (copula=="joe"){
    coppar=2
    cop1 <- archmCopula("joe",coppar, dim = m)}
  if (copula=="clayton"){
    coppar=2
    cop1 <- archmCopula("clayton",coppar, dim = m)}

  if (copula=="gumbel"){
    coppar=2
    cop1 <- archmCopula("gumbel",coppar, dim = m)}
  if (copula=="student"){
    coppar=c(rep(0.5,1),5)
    cop1 <- tCopula(coppar[1], dim=m,df = coppar[2],dispstr = "un") }
  if (copula=="normal"|copula=="student"|copula=="gumbel"|copula=="clayton"|copula=="joe"|copula=="frank"){
    copup = cop1@param.upbnd
    coplow = cop1@param.lowbnd
  }
  if (copula=="NS"){
    coppar=c(0.5,0.5,0.5)
    coplow=c(-0.3,-0.3,0.001)
    copup=c(0.9,0.9,0.99)
  }
  if (copula=="NC"){
    coppar=c(0.5,2,0.5)
    coplow=c(-0.3,0.1,0.001)
    copup=c(0.9,Inf,0.99)
  }
  if (copula=="NF"){
    coppar=c(0.5,2,0.5)
    coplow=c(-0.3,0.1,0.001)
    copup=c(0.9,Inf,0.99)
  }
  if (copula=="NG"){
    coppar=c(0.5,2,0.5)
    coplow=c(-0.3,0.1,0.001)
    copup=c(0.9,Inf,0.99)
  }
  if (copula=="NJ"){
    coppar=c(0.5,2,0.5)
    coplow=c(-0.3,1.1,0.001)
    copup=c(0.9,Inf,0.99)
  }
  if (copula=="SC"){
    coppar=c(0.5,2,0.5)
    coplow=c(-0.3,0.1,0.001)
    copup=c(0.9,Inf,0.99)
  }
  if (copula=="SF"){
    coppar=c(0.5,2,0.5)
    coplow=c(-0.3,0.1,0.001)
    copup=c(0.9,Inf,0.99)
  }
  if (copula=="SG"){
    coppar=c(0.5,2,0.5)
    coplow=c(-0.3,0.1,0.001)
    copup=c(0.9,Inf,0.99)
  }
  if (copula=="SJ"){
    coppar=c(0.5,2,0.5)
    coplow=c(-0.3,1.1,0.001)
    copup=c(0.9,Inf,0.99)
  }
  if (copula=="CF"){
    coppar=c(2,2,0.5)
    coplow=c(0.1,0.1,0.001)
    copup=c(Inf,Inf,0.99)
  }
  if (copula=="CG"){
    coppar=c(2,2,0.5)
    coplow=c(0.1,1.1,0.01)
    copup=c(Inf,Inf,0.99)
  }
  if (copula=="CJ"){
    coppar=c(5,5,0.5)
    coplow=c(0.1,1.1,0.001)
    copup=c(Inf,Inf,0.99)
  }
  if (copula=="FG"){
    coppar=c(2,2,0.5)
    coplow=c(0.1,1.1,0.001)
    copup=c(Inf,Inf,0.99)
  }

  if (copula=="FJ"){
    coppar=c(2,2,0.5)
    coplow=c(0.1,1.1,0.001)
    copup=c(Inf,Inf,0.99)
  }
  if (copula=="GJ"){
    coppar=c(4,10,0.5)
    coplow=c(1.5,1.1,0.001)
    copup=c(Inf,Inf,0.99)
  }
  ##=========================



  ##MLE starting value
  colnames(Y)=c("y2","y1")
  varmodel <- VAR(Y , p = p, type = "const")
  Dist=coef(varmodel)
  AIC(varmodel)
  BIC(varmodel)

  Dist1=rev(Dist$y1[,1])
  Dist2=rev(Dist$y2[,1])


  sigma=rep(1,ncol(Y))
  lowcoef=rep(-Inf,length(Dist1))
  upcoef=rep(Inf,length(Dist1))

  param=c(Dist1,Dist2,sigma=sigma,cop_par=coppar)
  lower =c(lowcoef,lowcoef,rep(0.000001,2),coplow+0.01)
  upper =c(upcoef,upcoef,rep(Inf,2),copup-0.01)



  model1 <- optim(par=param,fn=copvar2,Y=Y,p=p,dist1=dist1,dist2=dist2,copula=copula,
                  control = list(maxit=100000,fnscale=-1),method="L-BFGS-B",
                  lower =lower,upper =upper, hessian=TRUE )

  #######################################

  se=sqrt(diag(solve(-model1$hessian)))
  coef=model1$par
  n=nrow(Y)
  BIC= -2*model1$value + length(coef)*log(n)
  AIC = -2*model1$value + 2*length(coef)
  model1$value
  stat=coef/se
  pvalue <- 2*(1 - pnorm(abs(stat)))
  result=cbind(coef,se,pvalue )

  output <- list(result= result,
                 AIC=AIC,
                 BIC=BIC,
                 loglike=model1$value
  )
  return(output)
}
