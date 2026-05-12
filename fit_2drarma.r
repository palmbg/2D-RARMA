fit.2drarma = function (y, ar = NA, ma = NA)
{
  source("functions.r")

  
  
  maxit1 = 1000
  
  ar = max(ar)
  ma = max(ma)
  p = max(ar) # AR order
  q = max(ma) # MA order
  n = dim(y)[1] # n is the number of rows
  k = dim(y)[2] # k is the number of columns 
  m = max(p,q,na.rm=T)
  
    ynew = round(log(y),4)

  # inicializacao dos parametros alpha e phi (beta)
  if(any(is.na(ar)==F)) # se nao tem componente ar
  {
    p1 = (p+1)^2-1
    
    XX = c()
    
    for (i in (m+1):n)
    {
      for (j in (m+1):k)
      {
        xx1 = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
        XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
      }
    }
    
    P = XX[,2:dim(XX)[2]]
    Y = as.matrix(XX[,1])

    Z = cbind(rep(1,(n-m)*(k-m)),P)
    
    x = as.matrix(Z)
    ajuste = lm.fit(x, Y)
    mqo = round(c(ajuste$coef),4)
    
  }else{
    ynew = log(y)
    Z = as.matrix(rep(1,(n-m)*(k-m)))
    
    Y = as.vector(ynew[(m+1):n,(m+1):k])
    
    x = as.matrix(Z)
    ajuste = lm.fit(x, Y)
    mqo = round(c(ajuste$coef),4)
  }
  
  
  ############
  
  ######### ARMA model
  if(any(is.na(ar)==F) && any(is.na(ma)==F))
  { 
    q1 = (q+1)^2-1
   
    reg = c(mqo, rep(0,q1)) # initializing the parameter values
    
    loglik = function(z) 
    {
      alpha = z[1]
      phi = z[2:(p1+1)] 
      theta = z[(p1+2):(p1+q1+1)]
      
      eta = error = matrix(0, ncol=k,nrow=n)
      
      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          xx = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
          y_new1 = as.matrix(xx[(length(xx)-1):1])
          
          yy = as.vector(t(error[(i-ma):i,(j-ma):j]))
          error_new = as.matrix(yy[(length(yy)-1):1])
          
          eta[i,j]  = round(alpha + phi%*%y_new1 + theta%*%error_new,4)
          error[i,j] = round(ynew[i,j]-eta[i,j],4)
        }
        
      }
      
      mu = round(as.vector(exp((t(eta[(m+1):n,(m+1):k])))),4)
      y1 = round(as.vector(t(y[(m+1):n,(m+1):k])),4)
      
      ll = round(suppressWarnings(-sum(log(pi/2)+log(y1)-log(2*mu^2)-(pi*y1^2)/(4*(mu^2)))),4)
      
      sum(ll)
  
    } 
    
    escore = function(z) 
    {
      alpha = z[1]
      phi = z[2:(p1+1)] 
      theta = z[(p1+2):(p1+q1+1)]
      
      eta = error = matrix(0, ncol=k,nrow=n)
      
      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          xx = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
          y_new1 = as.matrix(xx[(length(xx)-1):1])
          
          yy = as.vector(t(error[(i-ma):i,(j-ma):j]))
          error_new = as.matrix(yy[(length(yy)-1):1])
          
          eta[i,j]  = round(alpha + phi%*%y_new1 + theta%*%error_new,4)
          error[i,j] = round(ynew[i,j]-eta[i,j],4)
        }
        
      }
      
      mu = round(as.vector(exp(t(eta[(m+1):n,(m+1):k]))),4)
      y1 = round(as.vector(t(y[(m+1):n,(m+1):k])),4)
  
      dmu = round(as.vector(((pi*(y1^2))/(2*(mu^3))-(2)/(mu))),4)
      
      mT = round(diag(mu),4)
      
      ###################################
      
      # dalpha
      deta.dalpha = matrix(0, nrow = n, ncol=k)
      
      for(i in (m+1):n)
      {
        for(j in (m+1):k)
        {
          xx = as.vector(t(deta.dalpha[(i-ma):i,(j-ma):j]))
          deta.dalpha1 = as.matrix(xx[(length(xx)-1):1])
          deta.dalpha[i,j] = 1 - theta%*%deta.dalpha1
        }
        
      }
      
      a = as.vector(deta.dalpha[(m+1):n,(m+1):k])
      
      # dphi 
      P1 = rbind(rep(0,p1),P)
      deta.dphi = matrix(0, ncol=p1,nrow=dim(P1)[1])
      dsum.phi = matrix(0, ncol=p1,nrow=q1)
    
      for(i in m:dim(P1)[1])
      {
        if( i == m)
        {
          deta.dphi[i,] = P1[i,] 
        }else{
          
          for(j in 1:q1)
          {
            dsum.phi[j,] = theta[j] * deta.dphi[i-ma,]
          }
          deta.dphi[i,] = P1[i,] - apply(dsum.phi,2,sum)
        }
      }
      
      rP = deta.dphi[-m,]
      
      #dtheta
      XX = c()
      
      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          xx1 = as.vector(t(error[(i-ma):i,(j-ma):j]))
          XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
        }
      }
      
      R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
      deta.dtheta = matrix(0, ncol=q1,nrow=dim(R)[1])
      dsum.theta = matrix(0, ncol=q1,nrow=q1)
      
        for(i in m:dim(R)[1])
        {
          if( i == m)
          {
            deta.dtheta[i,] = R[i,] 
          }else{
            for(j in 1:q1)
            {
              dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
            }
            deta.dtheta[i,] = R[i,] - apply(dsum.theta,2,sum)
          }
        }
      rR = deta.dtheta[-m,]
      
      ###################################
      
      Ualpha =  - a %*% mT %*% dmu
      Uphi =    - t(rP) %*% mT %*% dmu
      Utheta =  - t(rR) %*% mT %*% dmu

      rval = round(c(Ualpha,Uphi,Utheta),4)
    }
    
    if(any(is.na(ar))==F) names_phi = c(paste("phi",1:p1,sep=""))
    
    if(any(is.na(ma))==F) names_theta = (paste("theta",1:q1,sep=""))
    
    names_par = c("alpha",names_phi,names_theta)
    
    opt = optim(reg, loglik, 
                 escore,
                 method = "BFGS", 
                 hessian = TRUE,
                 control = list(maxit = 400, abstol = 1e-6,
                                factr = 1e20))
    
    
    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")
    
    z = c()
    z$conv = opt$conv
    coef = opt$par
    names(coef) = names_par
    z$coeff = round(coef,4)

    alpha = coef[1]
    phi = coef[2:(p1+1)]
    theta = coef[(p1+2):(p1+q1+1)]


    z$alpha = alpha
    z$phi = phi
    z$theta = theta
    
    etahat = errorhat = matrix(0, ncol=k,nrow=n)
    
    for (i in (m+1):n)
    {
      for (j in (m+1):k)
      {
        xx = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
        y_new1 = as.matrix(xx[(length(xx)-1):1])
        
        yy = as.vector(t(errorhat[(i-ma):i,(j-ma):j]))
        errorhat_new = as.matrix(yy[(length(yy)-1):1])
        
        etahat[i,j]  = round(alpha + phi%*%y_new1 + theta%*%errorhat_new,4)
        errorhat[i,j] = round(ynew[i,j]-etahat[i,j],4)
      }
      
    }

    z$fitted = round(exp(etahat[(m+1):n,(m+1):k]),4)
    z$etahat = round(etahat[(m+1):n,(m+1):k],4)
    z$errorhat = round(errorhat[(m+1):n,(m+1):k],4)
  
    y1 = round(as.vector(t(y[(m+1):n,(m+1):k])),4)
    
    ###################################
    
    # dalpha
    deta.dalpha = matrix(0, nrow = n, ncol=k)
    
    for(i in (m+1):n)
    {
      for(j in (m+1):k)
      {
        xx = as.vector(t(deta.dalpha[(i-ma):i,(j-ma):j]))
        deta.dalpha1 = as.matrix(xx[(length(xx)-1):1])
        deta.dalpha[i,j] = 1 - theta%*%deta.dalpha1
      }
      
    }
    
    a = as.vector(deta.dalpha[(m+1):n,(m+1):k])
    
    # dphi 
    P1 = rbind(rep(0,p1),P)
    deta.dphi = matrix(0, ncol=p1,nrow=dim(P1)[1])
    dsum.phi = matrix(0, ncol=p1,nrow=q1)
    
    for(i in m:dim(P1)[1])
    {
      if( i == m)
      {
        deta.dphi[i,] = P1[i,] 
      }else{
        
        for(j in 1:q1)
        {
          dsum.phi[j,] = theta[j] * deta.dphi[i-ma,]
        }
        deta.dphi[i,] = P1[i,] - apply(dsum.phi,2,sum)
      }
    }
    
    rP = deta.dphi[-m,]
    
    ##dtheta
    XX = c()
    
    for (i in (m+1):n)
    {
      for (j in (m+1):k)
      {
        xx1 = as.vector(t(errorhat[(i-ma):i,(j-ma):j]))
        XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
      }
    }
    
    R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
    deta.dtheta = matrix(0, ncol=q1,nrow=dim(R)[1])
    dsum.theta = matrix(0, ncol=q1,nrow=q1)
    
    for(i in m:dim(R)[1])
    {
      if( i == m)
      {
        deta.dtheta[i,] = R[i,] 
      }else{
        for(j in 1:q1)
        {
          dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
        }
        deta.dtheta[i,] = R[i,] - apply(dsum.theta,2,sum)
      }
    }
    rR = deta.dtheta[-m,]
    
    ###################################
    
    muhat2 = round(as.vector(exp(t(etahat[(m+1):n,(m+1):k]))),4)
    W = round(diag(((4)/(muhat2^2))*(muhat2^2)),4)

    Kaa = t(a) %*% W %*% a
    Kpa = t(rP) %*% W %*% a
    Kap = t(Kpa)
    Kta = t(rR) %*% W %*% a
    Kat = t(Kta)
    Kpp = t(rP) %*% W %*% rP
    Kpt = t(rP) %*% W %*% rR
    Ktp = t(Kpt)
    Ktt = t(rR) %*% W %*% rR

    K = rbind(
      cbind(Kaa,Kap,Kat),
      cbind(Kpa,Kpp,Kpt),
      cbind(Kta,Ktp,Ktt)
    )
    
    K = round(K,4)
    

  }  
  
  
  ##### only AR model
  if(any(is.na(ar)==F) && any(is.na(ma)==T))
  {
    
    p1 = (p+1)^2-1
    q1 = 0
    reg = c(mqo) # initializing the parameter values
 
    
    loglik = function(z) 
    {
      alpha = z[1]
      phi = z[2:(p1+1)] 
      
      eta = matrix(0, ncol=k,nrow=n)
      
      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          xx = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
          y_new1 = as.matrix(xx[(length(xx)-1):1])
          eta[i,j]  = round(alpha + phi%*%y_new1,4) 
        }
        
      }
      
      mu = round(as.vector(exp((t(eta[(m+1):n,(m+1):k])))),4)
      y1 = round(as.vector(t(y[(m+1):n,(m+1):k])),4)
      
      ll = round(suppressWarnings(-sum(log(pi/2)+log(y1)-log(2*mu^2)-(pi*y1^2)/(4*(mu^2)))),4)
      
      sum(ll)
      
    } 
    
    escore = function(z) 
    {
      alpha = z[1]
      phi = z[2:(p1+1)] 
      
      eta = matrix(0, ncol=k,nrow=n)
      
      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          xx = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
          y_new1 = as.matrix(xx[(length(xx)-1):1])
          eta[i,j]  = round(alpha + phi%*%y_new1,4)
        }
        
      }
      
      mu = round(as.vector(exp(t(eta[(m+1):n,(m+1):k]))),4)
      y1 = round(as.vector(t(y[(m+1):n,(m+1):k])),4)
      
      dmu = round(as.vector(((pi*(y1^2))/(2*(mu^3))-(2)/(mu))),4)
      
      mT = round(diag(mu),4)
      a = as.vector(rep(1,(n-m)*(k-m)))
      
      Ualpha =  - t(a) %*% mT %*% dmu
      Uphi =    - t(P) %*% mT %*% dmu
      
      rval = round(c(Ualpha,Uphi),4)
    }
    
    if(any(is.na(ar))==F) names_phi = c(paste("phi",1:p1,sep=""))
    
    names_par = c("alpha",names_phi)
    
    opt = optim(reg, loglik, 
                escore,
                method = "BFGS", 
                hessian = TRUE,
                control = list(maxit = 400, abstol = 1e-6,
                               factr = 1e20))
    
    
    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")
    
    z = c()
    z$conv = opt$conv
    coef = (opt$par)
    names(coef) = names_par
    z$coeff = round(coef,4)
    
    alpha = coef[1]
    phi = coef[2:(p1+1)]

    z$alpha = alpha
    z$phi = phi
    
    etahat = matrix(0, ncol=k,nrow=n)
    
    for (i in (m+1):n)
    {
      for (j in (m+1):k)
      {
        xx = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
        y_new1 = as.matrix(xx[(length(xx)-1):1])
        etahat[i,j]  = round(alpha + phi%*%y_new1,4)
  
      }
      
    }
    
    
    z$fitted = round(exp(etahat[(m+1):n,(m+1):k]),4)
    z$etahat = round(etahat[(m+1):n,(m+1):k],4)
    
    y1 = round(as.vector(t(y[(m+1):n,(m+1):k])),4)
    
    muhat2 = round(as.vector(exp(t(etahat[(m+1):n,(m+1):k]))),4)
    W = round(diag(((4)/(muhat2^2))*(muhat2^2)),4)
    
    a = as.vector(rep(1,(n-m)*(k-m)))
    
    Kaa = t(a) %*% W %*% a
    Kpa = t(P) %*% W %*% a
    Kap = t(Kpa)
    Kpp = t(P) %*% W %*% P
    
    K = rbind(
      cbind(Kaa,Kap),
      cbind(Kpa,Kpp)
    )
    
    K = round(K,4)

    
  }
  
  ######### MA model
  if(any(is.na(ar)==T) && any(is.na(ma)==F))
  { 
    p1 = 0
    q1 = (q+1)^2-1
    reg = c(mqo[1],rep(0,q1)) # initializing the parameter values
    
    loglik = function(z) 
    {
      alpha = z[1]
      theta = z[2:(q1+1)]
      
      eta = error = matrix(0, ncol=k,nrow=n)
      
      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          yy = as.vector(t(error[(i-ma):i,(j-ma):j]))
          error_new = as.matrix(yy[(length(yy)-1):1])
          
          eta[i,j]  = round(alpha + theta%*%error_new,4)
          error[i,j] = round(ynew[i,j]-eta[i,j],4)
        }
        
      }
      
      mu = round(as.vector(exp((t(eta[(m+1):n,(m+1):k])))),4)
      y1 = round(as.vector(t(y[(m+1):n,(m+1):k])),4)
      
      ll = round(suppressWarnings(-sum(log(pi/2)+log(y1)-log(2*mu^2)-(pi*y1^2)/(4*(mu^2)))),4)
      
      sum(ll)
      
    } 
    
    escore = function(z) 
    {
      alpha = z[1]
      theta = z[2:(q1+1)]
      
      eta = error = matrix(0, ncol=k,nrow=n)
      
      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          
          yy = as.vector(t(error[(i-ma):i,(j-ma):j]))
          error_new = as.matrix(yy[(length(yy)-1):1])
          
          eta[i,j]  = round(alpha + theta%*%error_new,4)
          error[i,j] = round(ynew[i,j]-eta[i,j],4)
        }
        
      }
      
      mu = round(as.vector(exp(t(eta[(m+1):n,(m+1):k]))),4)
      y1 = round(as.vector(t(y[(m+1):n,(m+1):k])),4)
      
      dmu = round(as.vector(((pi*(y1^2))/(2*(mu^3))-(2)/(mu))),4)
      
      mT = round(diag(mu),4)
      
      ###################################
      
      # dalpha
      deta.dalpha = matrix(0, nrow = n, ncol=k)
      
      for(i in (m+1):n)
      {
        for(j in (m+1):k)
        {
          xx = as.vector(t(deta.dalpha[(i-ma):i,(j-ma):j]))
          deta.dalpha1 = as.matrix(xx[(length(xx)-1):1])
          deta.dalpha[i,j] = 1 - theta%*%deta.dalpha1
        }
        
      }
      
      a = as.vector(deta.dalpha[(m+1):n,(m+1):k])

      #dtheta
      XX = c()
      
      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          xx1 = as.vector(t(error[(i-ma):i,(j-ma):j]))
          XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
        }
      }
      
      R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
      deta.dtheta = matrix(0, ncol=q1,nrow=dim(R)[1])
      dsum.theta = matrix(0, ncol=q1,nrow=q1)
      
      for(i in m:dim(R)[1])
      {
        if( i == m)
        {
          deta.dtheta[i,] = R[i,] 
        }else{
          for(j in 1:q1)
          {
            dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
          }
          
          deta.dtheta[i,] = R[i,] - apply(dsum.theta,2,sum)
        }
      }
      rR = deta.dtheta[-m,]
      
      ###################################
      
      Ualpha =  - a %*% mT %*% dmu
      Utheta =  - t(rR) %*% mT %*% dmu
      
      rval = round(c(Ualpha,Utheta),4)
    }
    
    if(any(is.na(ma))==F) names_theta = (paste("theta",1:q1,sep=""))
    
    names_par = c("alpha",names_theta)
    
    opt = optim(reg, loglik, 
                escore,
                method = "BFGS", 
                hessian = TRUE,
                control = list(maxit = 400, abstol = 1e-6,
                               factr = 1e20))
    
    
    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")
    
    z = c()
    z$conv = opt$conv
    coef = opt$par
    names(coef) = names_par
    z$coeff = round(coef,4)
    
    alpha = coef[1]
    theta = coef[2:(q1+1)]
    
    z$alpha = alpha
    z$theta = theta
    
    etahat = errorhat = matrix(0, ncol=k,nrow=n)
    
    for (i in (m+1):n)
    {
      for (j in (m+1):k)
      {
      
        yy = as.vector(t(errorhat[(i-ma):i,(j-ma):j]))
        errorhat_new = as.matrix(yy[(length(yy)-1):1])
        
        etahat[i,j]  = round(alpha + theta%*%errorhat_new,4)
        errorhat[i,j] = round(ynew[i,j]-etahat[i,j],4)
      }
      
    }
    
    
    z$fitted = round(exp(etahat[(m+1):n,(m+1):k]),4)
    z$etahat = round(etahat[(m+1):n,(m+1):k],4)
    z$errorhat = round(errorhat[(m+1):n,(m+1):k],4)
    
    y1 = round(as.vector(t(y[(m+1):n,(m+1):k])),4)
    
    ###################################
    
    # dalpha
    deta.dalpha = matrix(0, nrow = n, ncol=k)
    
    for(i in (m+1):n)
    {
      for(j in (m+1):k)
      {
        xx = as.vector(t(deta.dalpha[(i-ma):i,(j-ma):j]))
        deta.dalpha1 = as.matrix(xx[(length(xx)-1):1])
        deta.dalpha[i,j] = 1 - theta%*%deta.dalpha1
      }
      
    }
    
    a = as.vector(deta.dalpha[(m+1):n,(m+1):k])
    
    #dtheta
    XX = c()
    
    for (i in (m+1):n)
    {
      for (j in (m+1):k)
      {
        xx1 = as.vector(t(errorhat[(i-ma):i,(j-ma):j]))
        XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
      }
    }
    
    R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
    deta.dtheta = matrix(0, ncol=q1,nrow=dim(R)[1])
    dsum.theta = matrix(0, ncol=q1,nrow=q1)
    
    for(i in m:dim(R)[1])
    {
      if( i == m)
      {
        deta.dtheta[i,] = R[i,] 
      }else{
        for(j in 1:q1)
        {
          dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
        }
        deta.dtheta[i,] = R[i,] - apply(dsum.theta,2,sum)
      }
    }
    rR = deta.dtheta[-m,]

    
    ###################################
    
    muhat2 = round(as.vector(exp(t(etahat[(m+1):n,(m+1):k]))),4)
    W = round(diag(((4)/(muhat2^2))*(muhat2^2)),4)
    
    Kaa = t(a) %*% W %*% a
    Kta = t(rR) %*% W %*% a
    Kat = t(Kta)
    Ktt = t(rR) %*% W %*% rR
    
    K = rbind(
      cbind(Kaa,Kat),
      cbind(Kta,Ktt)
    )
    
    K = round(K,4)
    
  }
  
  z$image = y
  y1 = y[(m+1):n,(m+1):k]
  z$fitted = exp(z$etahat)

  ##############################################

  # residuals

  z$resid = round(qnorm(pr(y1,z$fitted)),4) # quantile residuals
  # 
  
  ############################################
  
  m1 = (dim(y)[1]*dim(y)[2]) - (dim(z$fitted)[1]*dim(z$fitted)[2])
  z$loglik1 = opt$value*((n*k)/((n*k)-m1))
  z$aic1 = -2*z$loglik1+2*length(z$coeff)
  z$aic2 = -2*opt$value+2*length(z$coeff)

  vcov = chol2inv(chol(K))
  z$vcov = round(vcov,4)

  stderror = sqrt(diag(vcov))
  z$stderror = round(stderror,4)

  z$zstat = round(abs(z$coef/stderror),4)
  z$pvalues = round(2*(1 - pnorm(z$zstat) ),4)

  model_presentation = cbind(round(z$coef,4),round(z$stderror,4),
                             round(z$zstat,4),round(z$pvalues,4))
  colnames(model_presentation)=c("Estimate","Std. Error","z value","Pr(>|z|)")

  z$model = model_presentation

  # 
   return(z)
  
  

}

